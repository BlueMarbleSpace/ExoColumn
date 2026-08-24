#!/usr/bin/env python3
"""
lbl_co2_benchmark.py — dense-CO2 (outer-HZ) line-by-line clear-sky OLR benchmark,
the maximum-greenhouse twin of tools/lbl_olr_benchmark.py.

Context: this is the OHZ analogue of Kopparapu et al. (2013) Figure 1 (their
"early Mars" 2-bar CO2 LBL-vs-SMART check), applied instead to *our own*
maximum-greenhouse limit column — pCO2 = 8.869 bar over 1 bar N2 at Ts = 273 K,
the Seff-minimum of reference/max_greenhouse (Seff = 0.395).  It produces an
independent line-by-line OLR for the SAME column ExoRT n68 sees, so the dense-CO2
F_IR offset can be split into "ExoRT n68 vs LBL" and (optionally) "CLIMA vs LBL".

Dense CO2 needs two opacity sources the IHZ tool (lbl_olr_benchmark.py) omits and
that RADIS does not supply; this module adds them, on top of CO2 Voigt lines:

  * CO2 sub-Lorentzian chi-factor (Perrin & Hartmann 1989) on the line wings,
    truncated at 500 cm-1 from each line centre — Kopparapu's exact prescription
    ("Lorentzian overestimates absorption in the far wings"; their §2.1).
    Coefficients via the ClearSky.jl implementation (Baum/Wordsworth lineage):
        B1 = 0.0888 - 0.16*exp(-0.0041 T);  B2 = 0.0526*exp(-0.00152 T)
        chi = 1                                        (dnu < 3 cm-1)
            = exp(-B1*(dnu-3))                         (3  <= dnu < 30)
            = exp(-B1*27 - B2*(dnu-30))                (30 <= dnu < 120)
            = exp(-B1*27 - B2*90 - 0.0232*(dnu-120))   (dnu >= 120)
  * CO2-CO2 collision-induced absorption (HITRAN CIA 2024, Karman et al.), the
    far-IR (1-750 cm-1) rototranslational band that sits under the 273 K Planck
    peak plus the weaker 1200-1800 cm-1 induced band.

Line data: RADIS/HITRAN2020 provides the per-line T/P-scaled centres, intensities
and widths (sf.df1); we replace only the lineshape convolution so the PH89 chi can
be applied per line.  In the dense layers the Lorentz HWHM (~0.7 cm-1) exceeds the
Doppler width (~5e-4 cm-1) by ~1300x, so a pure-Lorentz x chi sum is used (the
saturated band core emits at the local T regardless of its sub-cm-1 shape, so the
negligible Doppler-core error does not affect OLR; verified vs RADIS Voigt at
chi=1, see tools/check_co2_lbl.py).

H2O (saturated at 273 K, a trace at ~10 bar) lines + MT_CKD are included for
completeness via the shared helpers; their contribution is small.

Result (10-2000 cm-1, the pCO2 = 8.87 bar Seff-minimum column):
  LBL pure-Lorentz (opaque bound)  = 44.1 W/m2
  ExoRT n68 (this work)            = 75.5 W/m2
  LBL PH89 chi  (transparent bound)= 78.5 W/m2
ExoRT n68 sits well INSIDE the CO2 wing-treatment envelope.  The entire spread
lives in the 875-1200 cm-1 (8-12 um) window: it is NOT a CIA effect — HITRAN-2024
CO2-CO2 CIA (which carries the Gruszka-Borysow far-IR + Baranov 7um bands) AND
ExoRT's own GB/Baranov-derived CIA are both ~0 across 875-1150 cm-1 — it is the
CO2 far-line-wing (sub-Lorentzian chi) treatment, which swings the column OLR by
~47 W/m2 between pure-Lorentzian and PH89 wings.  This is the documented dense-CO2
window-continuum uncertainty (Yang et al. 2016); ExoRT n68's window opacity is
consistent with an intermediate wing treatment, between the two LBL bounds.
Both wing bounds are computed and saved (olr_nu_full = PH89 chi, olr_nu_full_nochi
= pure-Lorentz); the 2-panel figure uses the PH89-chi case as the reference LBL.

Usage:
  python tools/lbl_co2_benchmark.py [exocol_out.nc] [--nlay 40] [--wstep 0.02]
                                    [--out PREFIX]

Shares planck_nu, schwarzschild_olr, MTCKD and load_column with
lbl_olr_benchmark.py (imported); only the CO2 line/CIA physics is new.
"""
import os, sys, argparse
import numpy as np
import netCDF4 as nc
from numba import njit, prange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lbl_olr_benchmark import (planck_nu, schwarzschild_olr, MTCKD, load_column,
                               H, C, KB, NA, G, MW_H2O, MW_CO2, DIFFUSIVITY)

CIA_NC = os.path.expanduser('~/.radisdb/cia/CO2-CO2_2024.cia')
CIA_URL = 'https://hitran.org/data/CIA/main/CO2-CO2_2024.cia'
# CIA spectral windows used for OLR (cm-1); non-overlapping to avoid double counting.
# Far-IR rototranslational band (dominant under the 273 K Planck peak) + the
# 1200-1800 induced band.  Higher HITRAN CIA windows (>2500 cm-1) are off the
# 0-2000 cm-1 OLR range and are skipped.
CIA_WINDOWS = [(1.0, 750.0), (1150.0, 1850.0)]


# ----------------------------------------------------------------------------
# CO2-CO2 CIA (HITRAN CIA 2024, Karman et al.)
# ----------------------------------------------------------------------------
class CO2CIA:
    """HITRAN-format CO2-CO2 CIA reader -> absco(wn, T) [cm^5 / molecule^2].
    tau = absco * n_CO2^2 [molec^2/cm^6] * L [cm]."""

    def __init__(self, path=CIA_NC, windows=CIA_WINDOWS):
        if not (os.path.isfile(path) and os.path.getsize(path) > 1e4):
            import urllib.request
            os.makedirs(os.path.dirname(path), exist_ok=True)
            print(f"downloading HITRAN CO2-CO2 CIA -> {path}")
            urllib.request.urlretrieve(CIA_URL, path)
        blocks = self._parse(path)
        # group blocks into the chosen windows; build (T, wn) grids per window
        self.win = []
        for w1, w2 in windows:
            sel = [b for b in blocks if abs(b['w1'] - w1) < 1 and abs(b['w2'] - w2) < 1]
            if not sel:
                continue
            sel.sort(key=lambda b: b['T'])
            Ts = np.array([b['T'] for b in sel])
            wn = sel[0]['wn']                      # common grid within a window
            K = np.vstack([np.interp(wn, b['wn'], b['k']) for b in sel])  # [nT,nwn]
            self.win.append(dict(w1=w1, w2=w2, T=Ts, wn=wn, K=K))

    @staticmethod
    def _parse(path):
        blocks = []
        with open(path) as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            l = lines[i]
            if 'CO2-CO2' in l:
                p = l.split()
                w1, w2, npts, T = float(p[1]), float(p[2]), int(p[3]), float(p[4])
                data = np.array([lines[i + 1 + j].split() for j in range(npts)], float)
                blocks.append(dict(w1=w1, w2=w2, T=T,
                                   wn=data[:, 0], k=data[:, 1]))
                i += npts + 1
            else:
                i += 1
        return blocks

    def absco(self, wn_out, T):
        """CIA cross-section [cm^5/molec^2] on wn_out at temperature T."""
        out = np.zeros_like(wn_out, float)
        for w in self.win:
            Ts = w['T']
            Tc = min(max(T, Ts[0]), Ts[-1])        # clamp to table T-range
            j = np.searchsorted(Ts, Tc).clip(1, len(Ts) - 1)
            f = (Tc - Ts[j - 1]) / (Ts[j] - Ts[j - 1])
            kT = (1 - f) * w['K'][j - 1] + f * w['K'][j]     # interp in T
            m = (wn_out >= w['wn'][0]) & (wn_out <= w['wn'][-1])
            out[m] += np.interp(wn_out[m], w['wn'], kT)      # interp in wn
        return np.clip(out, 0.0, None)


# ----------------------------------------------------------------------------
# CO2 line absorption with Perrin & Hartmann (1989) sub-Lorentzian chi-factor
# ----------------------------------------------------------------------------
@njit(parallel=True, fastmath=True, cache=True)
def _lorentz_chi_sum(wn, nu0, S, gL, B1, B2, dnu_cut, chi_on):
    """k(wn) [cm2/molec] = sum_i S_i * Lorentz(wn-nu0_i; gL_i) * chi_PH89(|dnu|).
    chi_on=0 sets chi=1 everywhere (pure Lorentz, for validation vs RADIS Voigt).
    Loop over output grid (prange) so each k[j] is written by one thread."""
    nwn = wn.size
    nline = nu0.size
    k = np.zeros(nwn)
    inv_pi = 1.0 / np.pi
    for j in prange(nwn):
        wj = wn[j]
        acc = 0.0
        for i in range(nline):
            d = wj - nu0[i]
            if d < 0.0:
                d = -d
            if d > dnu_cut:
                continue
            chi = 1.0
            if chi_on != 0:
                if d < 3.0:
                    chi = 1.0
                elif d < 30.0:
                    chi = np.exp(-B1 * (d - 3.0))
                elif d < 120.0:
                    chi = np.exp(-B1 * 27.0 - B2 * (d - 30.0))
                else:
                    chi = np.exp(-B1 * 27.0 - B2 * 90.0 - 0.0232 * (d - 120.0))
            g = gL[i]
            acc += S[i] * g * inv_pi / (d * d + g * g) * chi
        k[j] = acc
    return k


def co2_line_tau(lay, wn, wstep, chi=True, dnu_cut=500.0, tau_min=1e-3):
    """tau(nlay, nwn) for CO2 lines with the PH89 chi-factor.  Uses RADIS to get
    the T/P-scaled HITRAN2020 line parameters per layer, then a numba Lorentz x chi
    sum.

    Weak-line pruning is done on PEAK COLUMN OPTICAL DEPTH, not on line strength
    relative to the band maximum.  A line is kept when

        tau_peak = S / (pi * gamma_L) * N_col  >  tau_min          (default 1e-3)

    i.e. when it could contribute more than tau_min of opacity over the whole CO2
    column.  The previous criterion (S > 1e-6 * S.max()) was wrong for dense
    columns in two ways.  (1) In absolute terms it scaled with whatever the
    strongest line in the *requested spectral range* happened to be, so widening
    the range from 80-1220 to 1-2500 cm-1 silently made it 12x harsher (S.max()
    jumps from the 15 um band head to the nu3 4.3 um band head).  (2) On the
    8.87-bar max-greenhouse column (N_col ~ 1.3e26 molec/cm2) a line sitting
    exactly at that threshold still had a peak optical depth of ~215 through the
    column, so the "negligible" lines being discarded were optically thick by two
    orders of magnitude.  The damage was concentrated in the 800-1200 cm-1
    window, where 17540 of 17610 in-window lines were being thrown away and no
    strong band exists to carry the opacity in their place: the window OLR was
    inflated 15.6 -> 28.2 W/m2 and the 10-2000 cm-1 total 78.5 -> 91.7 W/m2.
    Convergence in tau_min is checked in tools/check_co2_lbl.py."""
    from radis import SpectrumFactory
    nlay = len(lay['T'])
    wmin, wmax = wn[0], wn[-1]
    # RADIS is used ONLY to populate df1 (per-line T/P-scaled centre, intensity,
    # Lorentz HWHM) — all grid-independent — so its own output grid is set coarse
    # to skip the (discarded) Voigt convolution.  The fine lineshape sum with the
    # PH89 chi is done on `wn` by _lorentz_chi_sum.
    sf = SpectrumFactory(wavenum_min=max(1.0, wmin - dnu_cut),
                         wavenum_max=wmax + dnu_cut, wstep=2.0, molecule='CO2',
                         isotope='1,2,3', pressure=1.0, mole_fraction=0.9,
                         truncation=5.0, verbose=0,
                         warnings={'AccuracyError': 'ignore',
                                   'AccuracyWarning': 'ignore'})
    sf.fetch_databank('hitran')
    tau = np.zeros((nlay, len(wn)))
    N_col = float(np.sum(lay['N_co2']))      # whole-column CO2, for the tau_min test
    for k in range(nlay):
        T = lay['T'][k]; p_mb = lay['p_mb'][k]; x = lay['x_co2'][k]
        # df1 carries T/P-scaled centre (shiftwav), intensity S, Lorentz HWHM
        sf.eq_spectrum(Tgas=T, pressure=p_mb / 1013.25, mole_fraction=x)
        df = sf.df1
        nu0 = df['shiftwav'].values.astype(np.float64)
        S = df['S'].values.astype(np.float64)
        gL = df['hwhm_lorentz'].values.astype(np.float64)
        # peak column optical depth this line could contribute, on its own
        tau_pk = S / (np.pi * np.maximum(gL, 1e-12)) * N_col
        keep = tau_pk > tau_min
        nu0, S, gL = nu0[keep], S[keep], gL[keep]
        B1 = 0.0888 - 0.16 * np.exp(-0.0041 * T)
        B2 = 0.0526 * np.exp(-0.00152 * T)
        kabs = _lorentz_chi_sum(wn, nu0, S, gL, np.float64(B1), np.float64(B2),
                                np.float64(dnu_cut), np.int64(1 if chi else 0))
        # kabs is cross-section [cm2/molec]; tau = kabs * N_co2 [molec/cm2]
        tau[k] = kabs * lay['N_co2'][k]
        print(f"  [CO2] layer {k+1:2d}/{nlay}  T={T:6.1f}  p={p_mb:9.2f} mb"
              f"  nlines={keep.sum():7d}  max(tau)={tau[k].max():.2e}", flush=True)
    return tau


def h2o_tau(lay, wn, wstep):
    """tau(nlay, nwn) for the trace H2O: RADIS/HITRAN2020 lines (Voigt, no chi —
    H2O is dilute here) + the AER MT_CKD continuum, interpolated onto wn.  The
    max-greenhouse column carries ~2 cm precipitable water, which leaves a small
    rotation-band/window signature even though CO2 dominates."""
    from radis import SpectrumFactory
    sf = SpectrumFactory(wavenum_min=max(1.0, wn[0] - 25.0), wavenum_max=wn[-1] + 25.0,
                         wstep=wstep, molecule='H2O', isotope='1,2,3', pressure=1.0,
                         mole_fraction=1e-3, truncation=25.0, verbose=0,
                         warnings={'AccuracyError': 'ignore', 'AccuracyWarning': 'ignore'})
    sf.fetch_databank('hitran')
    ck = MTCKD()
    tau_l = np.zeros((len(lay['T']), len(wn)))
    tau_c = np.zeros_like(tau_l)
    for k in range(len(lay['T'])):
        T = lay['T'][k]; p_mb = lay['p_mb'][k]
        x = max(lay['x_h2o'][k], 1e-12)
        s = sf.eq_spectrum(Tgas=T, pressure=p_mb / 1013.25, mole_fraction=x)
        w, kabs = s.get('abscoeff', wunit='cm-1')          # cm-1 (= sigma * n_h2o)
        n_h2o = (p_mb * 100.0) / (KB * T) * 1e-6 * x        # molec/cm3
        sigma = np.interp(wn, w, kabs) / max(n_h2o, 1e-30)  # cm2/molec
        tau_l[k] = sigma * lay['N_h2o'][k]
        tau_c[k] = ck.absco(wn, p_mb, T, lay['x_h2o'][k]) * lay['N_h2o'][k]
    return tau_l, tau_c


def co2_cia_tau(lay, wn):
    """tau(nlay, nwn) for CO2-CO2 CIA."""
    cia = CO2CIA()
    tau = np.zeros((len(lay['T']), len(wn)))
    for k in range(len(lay['T'])):
        T = lay['T'][k]; p_mb = lay['p_mb'][k]; x = lay['x_co2'][k]
        n_co2 = (p_mb * 100.0) / (KB * T) * 1e-6 * x      # molec/cm3
        path = lay['N_co2'][k] / max(n_co2, 1e-30)        # cm  (N = n*L)
        tau[k] = cia.absco(wn, T) * n_co2 * n_co2 * path
    return tau


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ncfile', nargs='?',
                    default=os.path.join(os.path.dirname(HERE), 'iofiles', 'exocol_out.nc'))
    ap.add_argument('--nlay', type=int, default=40)
    ap.add_argument('--wmin', type=float, default=10.0)
    ap.add_argument('--wmax', type=float, default=2000.0)
    ap.add_argument('--wstep', type=float, default=0.02)
    ap.add_argument('--g', type=float, default=G,
                    help='surface gravity [m/s2] (default: ExoRT Earth value)')
    ap.add_argument('--out', default=os.path.join(
        os.path.dirname(HERE), 'reference', 'max_greenhouse', 'lbl_olr_co2_maxgh'))
    args = ap.parse_args()

    lay = load_column(args.ncfile, args.nlay, g=args.g)
    print(f"column: Ts={lay['ts']:.1f} K, {len(lay['T'])} layers, g={args.g} m/s2, "
          f"x_co2 sfc={lay['x_co2'][-1]:.4f}, x_h2o sfc={lay['x_h2o'][-1]:.3e}")

    wn = np.arange(args.wmin, args.wmax + args.wstep, args.wstep)
    print("CO2 lines (PH89 chi, sub-Lorentzian) ...", flush=True)
    tau_co2 = co2_line_tau(lay, wn, args.wstep, chi=True)
    print("CO2 lines (pure Lorentz, chi off) ...", flush=True)
    tau_co2_noX = co2_line_tau(lay, wn, args.wstep, chi=False)
    print("CO2-CO2 CIA (HITRAN 2024 = GB far-IR + Baranov 7um) ...", flush=True)
    tau_cia = co2_cia_tau(lay, wn)
    print("H2O (trace) lines + MT_CKD ...", flush=True)
    tau_h2o_l, tau_h2o_c = h2o_tau(lay, wn, args.wstep)
    tau_h2o = tau_h2o_l + tau_h2o_c

    # The dense-CO2 8-12 um window OLR is controlled by the CO2 far-wing treatment,
    # NOT by CIA (HITRAN 2024 + ExoRT's own GB/Baranov CIA are both ~0 in 875-1150).
    # The two wing treatments bracket the LBL: pure-Lorentz (opaque) <-> PH89 (open).
    cases = {
        'PH89 chi  + CIA + H2O (transparent bound)': tau_co2 + tau_cia + tau_h2o,
        'pure-Lorentz + CIA + H2O (opaque bound)': tau_co2_noX + tau_cia + tau_h2o,
    }
    print(f"\n=== dense-CO2 LBL OLR, Ts={lay['ts']:.0f} K column ===")
    results = {}
    for name, tau in cases.items():
        olr_nu = schwarzschild_olr(wn, tau, lay['T'], lay['ts'])
        olr = np.trapezoid(olr_nu, wn)
        results[name] = (olr, olr_nu)
        print(f"  {name:28s}: OLR = {olr:7.2f} W/m2")

    band = lay['band']
    if band:
        edges = band['edges']
        olr_nu = results['PH89 chi  + CIA + H2O (transparent bound)'][1]
        olr_nu_op = results['pure-Lorentz + CIA + H2O (opaque bound)'][1]
        rows = []
        for i in range(len(band['olr'])):
            n1, n2 = edges[i], edges[i + 1]
            if n2 < args.wmin or n1 > args.wmax:
                lbl_b = np.nan
            else:
                m = (wn >= n1) & (wn < n2)
                lbl_b = np.trapezoid(olr_nu[m], wn[m]) if m.sum() > 2 else 0.0
            rows.append((n1, n2, lbl_b, band['olr'][i]))
        print(f"\n  per-n68-band (W/m2): LBL vs ExoRT band_lwup_toa")
        for n1, n2, lbl_b, exo_b in rows:
            if np.isfinite(lbl_b) and (lbl_b > 0.3 or exo_b > 0.3):
                print(f"    {n1:7.0f}-{n2:7.0f}  LBL={lbl_b:7.2f}  ExoRT={exo_b:7.2f}"
                      f"  d={exo_b-lbl_b:+6.2f}")
        lbl_PH = sum(l for n1, n2, l, e in rows if np.isfinite(l))
        lbl_op = float(np.nansum([np.trapezoid(olr_nu_op[(wn >= n1) & (wn < n2)],
                       wn[(wn >= n1) & (wn < n2)]) if ((wn >= n1) & (wn < n2)).sum() > 2
                       else np.nan for n1, n2, l, e in rows]))
        exo_t = sum(e for n1, n2, l, e in rows if np.isfinite(l))
        print(f"\n  common-range totals ({args.wmin:.0f}-{args.wmax:.0f} cm-1):")
        print(f"    LBL pure-Lorentz (opaque)  = {lbl_op:.2f} W/m2")
        print(f"    ExoRT n68 (this work)      = {exo_t:.2f} W/m2")
        print(f"    LBL PH89 chi (transparent) = {lbl_PH:.2f} W/m2")
        print(f"    => ExoRT sits {'INSIDE' if lbl_op<=exo_t<=lbl_PH else 'OUTSIDE'} "
              f"the CO2 wing-treatment envelope [{lbl_op:.1f}, {lbl_PH:.1f}]")
        np.savez_compressed(
            args.out + '.npz', wn=wn.astype(np.float32),
            olr_nu_full=olr_nu.astype(np.float32),          # PH89 chi (transparent)
            olr_nu_full_nochi=olr_nu_op.astype(np.float32),  # pure Lorentz (opaque)
            band_edges=edges, band_exo=band['olr'],
            T=lay['T'], p_mb=lay['p_mb'], x_co2=lay['x_co2'],
            x_h2o=lay['x_h2o'], ts=lay['ts'])
        print(f"\nsaved {args.out}.npz")


if __name__ == '__main__':
    main()
