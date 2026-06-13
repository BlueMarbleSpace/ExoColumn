#!/usr/bin/env python3
"""
check_co2_lbl.py — validation harness for the dense-CO2 LBL physics in
lbl_co2_benchmark.py.  Three checks on a single representative dense-CO2 layer:

  1. chi OFF (pure Lorentz) vs RADIS Voigt abscoeff  -> our per-line sum reproduces
     RADIS (expected agreement to ~1%, since Doppler << Lorentz here).
  2. chi ON vs OFF -> the PH89 sub-Lorentzian wing suppression is visible.
  3. CO2-CO2 CIA magnitude at this layer is physically sized.

Run: python tools/check_co2_lbl.py
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lbl_co2_benchmark import _lorentz_chi_sum, CO2CIA, KB
from radis import SpectrumFactory

# representative dense layer near the middle of the max-greenhouse column
T, p_mb, x = 290.0, 4000.0, 0.90          # ~4 bar, 90% CO2
wmin, wmax, wstep, cut = 580.0, 720.0, 0.01, 500.0
n_co2 = (p_mb * 100.0) / (KB * T) * 1e-6 * x      # molec/cm3

sf = SpectrumFactory(wavenum_min=wmin - cut, wavenum_max=wmax + cut, wstep=wstep,
                     molecule='CO2', isotope='1,2,3', pressure=1.0, mole_fraction=x,
                     truncation=cut, verbose=0,
                     warnings={'AccuracyError': 'ignore', 'AccuracyWarning': 'ignore'})
sf.fetch_databank('hitran')
s = sf.eq_spectrum(Tgas=T, pressure=p_mb / 1013.25, mole_fraction=x)
wn_full, abscoeff = s.get('abscoeff', wunit='cm-1')     # cm-1 (= sigma * n_co2)

# our grid (subset) and per-line params from df1
m = (wn_full >= wmin) & (wn_full <= wmax)
wn = wn_full[m]
radis_sigma = abscoeff[m] / n_co2                       # cm2/molec (RADIS Voigt)

df = sf.df1
nu0 = df['shiftwav'].values.astype(np.float64)
S = df['S'].values.astype(np.float64)
gL = df['hwhm_lorentz'].values.astype(np.float64)
keep = S > 1e-6 * S.max()
nu0, S, gL = nu0[keep], S[keep], gL[keep]
B1 = 0.0888 - 0.16 * np.exp(-0.0041 * T)
B2 = 0.0526 * np.exp(-0.00152 * T)

sig_off = _lorentz_chi_sum(wn, nu0, S, gL, np.float64(B1), np.float64(B2),
                           np.float64(cut), np.int64(0))   # chi=1, pure Lorentz
sig_on = _lorentz_chi_sum(wn, nu0, S, gL, np.float64(B1), np.float64(B2),
                          np.float64(cut), np.int64(1))    # chi PH89

# check 1: chi-off vs RADIS Voigt
rel = np.abs(sig_off - radis_sigma) / (radis_sigma + 1e-30)
band = radis_sigma > radis_sigma.max() * 1e-4            # where there is signal
print("=== check 1: Lorentz(chi=1) vs RADIS Voigt, %d lines kept ===" % keep.sum())
print(f"  median rel.diff where signal: {np.median(rel[band])*100:.2f}%")
print(f"  mean   rel.diff where signal: {np.mean(rel[band])*100:.2f}%")
print(f"  band-integrated sigma: ours={np.trapezoid(sig_off,wn):.4e} "
      f"RADIS={np.trapezoid(radis_sigma,wn):.4e}  "
      f"ratio={np.trapezoid(sig_off,wn)/np.trapezoid(radis_sigma,wn):.4f}")

# check 2: chi effect (sample a window point 619 cm-1, between band lines)
print("\n=== check 2: PH89 chi suppression ===")
for w0 in (605.0, 619.0, 700.0):
    j = np.argmin(np.abs(wn - w0))
    print(f"  nu={wn[j]:.1f} cm-1: sigma chi-off={sig_off[j]:.3e}  "
          f"chi-on={sig_on[j]:.3e}  ratio={sig_on[j]/(sig_off[j]+1e-40):.3f}")
print(f"  band-integrated sigma: chi-on/chi-off = "
      f"{np.trapezoid(sig_on,wn)/np.trapezoid(sig_off,wn):.4f}")

# check 3: CIA magnitude
print("\n=== check 3: CO2-CO2 CIA ===")
cia = CO2CIA()
wn_cia = np.arange(10.0, 2000.0, 1.0)
absco = cia.absco(wn_cia, T)
jpk = np.argmax(absco)
print(f"  windows loaded: {[(w['w1'],w['w2'],len(w['T'])) for w in cia.win]}")
print(f"  peak CIA cross-section {absco[jpk]:.3e} cm^5/molec^2 at "
      f"{wn_cia[jpk]:.0f} cm-1 (T={T:.0f}K)")
# CIA absorption coefficient at peak for this layer: alpha = absco * n_co2^2
print(f"  -> alpha_CIA at peak = {absco[jpk]*n_co2**2:.3e} cm-1 "
      f"(tau over 10 km = {absco[jpk]*n_co2**2*1e6:.2f})")
