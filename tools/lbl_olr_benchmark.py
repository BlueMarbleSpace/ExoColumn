#!/usr/bin/env python3
"""
lbl_olr_benchmark.py — line-by-line clear-sky OLR benchmark for an ExoColumn
column, to attribute the IHZ F_IR offset vs Kopparapu et al. (2013).

Context (see reference/moist_runaway/README.md, "LW offset" bullet): with the
atmospheric state verified to match CLIMA (T(P) ±0.07 K, H2O profiles 0-2%),
ExoColumn's OLR exceeds Kopparapu's tabulated FTIR by a bell that peaks at
~+20 W/m² at Ts = 300-320 K and vanishes at the steam plateau.  Yang et al.
(2016, ApJ 826:222) documented 10-25 W/m² OLR spreads among band models on
exactly this configuration and attributed them to window-region opacity.
This script computes an independent line-by-line OLR for the SAME column so
the offset can be split into "ExoRT n68 vs LBL" and "CLIMA vs LBL" parts.

Ingredients:
  * Lines: HITRAN for H2O and CO2 (via RADIS 0.17's fetch_databank('hitran'),
    cached in ~/.radisdb), Voigt profiles, RADIS default lineshape truncation;
    air-broadened widths (the N2 background is approximated as air — same
    approximation as the k-tables).  Line-list edition (methods note): RADIS
    serves HITRAN's *current* consolidated edition, HITRAN2020 — hitran.org no
    longer exposes the 2016 line-by-line edition through any of RADIS, astroquery,
    HAPI or the LBL web interface, and no 2016 .par files are cached locally.
    The ExoRT n68 k-tables this is benchmarked against are built from HITRAN2016,
    but the 2016->2020 updates to the H2O rotation/nu2 bands that dominate the
    10-3000 cm-1 OLR are sub-1% (and CO2 at 330 ppm contributes only the
    well-established 15 um band).  So the line-list edition is a sub-leading term,
    smaller than the n68-LBL residual (272.4 vs 269.7 W/m2, ~1%): that residual is
    the correlated-k (k-distribution) approximation, not the line list.
  * H2O continuum: faithful Python port of the AER MT_CKD water-vapour
    continuum (mt_ckd_h2o_module.f90; self + foreign, T-power-law on the self
    coefficient, (p/p_ref)(T_ref/T) density scaling, FASCOD radiation term).
    The coefficient file absco-ref_wv-mt-ckd.nc is downloaded from the AER-RC
    GitHub on first use and cached next to this script (© AER, redistributed
    per their research-use license with acknowledgment).
  * RT: no-scattering Schwarzschild integration with the standard diffusivity
    factor D = 1.66; isothermal layers; surface emissivity 1 at Ts.

Limitations (documented, all small for the 300 K target): no N2-N2 CIA
(matters < 1 W/m² at 1 bar under the saturated H2O rotation band); no CO2
chi-factor sub-Lorentzian wings (matters for dense CO2, not 330 ppm); RADIS
wstep under-resolves Doppler cores above ~50 Pa (the stratosphere is
isothermal at t_strato, so saturated-core emission temperature is unaffected).

Usage:
  python tools/lbl_olr_benchmark.py [exocol_out.nc] [--nlay 50] [--wstep 0.01]
                                    [--out PREFIX]

The input must be an ExoColumn output file (needs tmid, pint, h2ommr, co2mmr,
mw, ts, plus band_lwup_toa/wavenum_edge for the n68 comparison).
"""
import os, sys, argparse
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
MTCKD_NC = os.path.join(HERE, 'absco-ref_wv-mt-ckd.nc')
MTCKD_URL = ('https://raw.githubusercontent.com/AER-RC/MT_CKD_H2O/master/'
             'data/absco-ref_wv-mt-ckd.nc')

# physical constants (SI)
H = 6.62607015e-34; C = 2.99792458e8; KB = 1.380649e-23
NA = 6.02214076e23
G = 9.80616            # ExoRT exo_g
MW_H2O = 18.015268e-3  # kg/mol
MW_CO2 = 44.0095e-3
RADCN2 = H * C / KB * 100.0   # second radiation constant [cm K] = 1.43877
DIFFUSIVITY = 1.66

def planck_nu(wn_cm, T):
    """pi*B_nu per wavenumber [W/m2/(cm-1)] at wavenumber wn_cm [cm-1]."""
    nu = wn_cm * 100.0                       # m-1
    B = 2 * H * C**2 * nu**3 / np.expm1(H * C * nu / (KB * T))  # W/m2/sr/m-1
    return np.pi * B * 100.0                 # per cm-1


# ----------------------------------------------------------------------------
# MT_CKD H2O continuum (port of AER mt_ckd_h2o_module.f90)
# ----------------------------------------------------------------------------
class MTCKD:
    def __init__(self, path=MTCKD_NC):
        if not os.path.isfile(path):
            import urllib.request
            print(f"downloading MT_CKD coefficients -> {path}")
            urllib.request.urlretrieve(MTCKD_URL, path)
        with nc.Dataset(path) as ds:
            self.wvn = np.array(ds['wavenumbers'][:], float)
            self.self_ref = np.array(ds['self_absco_ref'][:], float)
            self.for_ref = np.array(ds['for_absco_ref'][:], float)
            self.self_texp = np.array(ds['self_texp'][:], float)
            self.p_ref = float(ds['ref_press'][:])   # 1013 mb
            self.t_ref = float(ds['ref_temp'][:])    # 296 K

    @staticmethod
    def _radfn(wn, T):
        """FASCOD radiation term [cm-1]."""
        xkt = T / RADCN2
        x = wn / xkt
        rad = np.where(x <= 0.01, 0.5 * x * wn,
                       wn * (1.0 - np.exp(-x)) / (1.0 + np.exp(-x)))
        return rad

    def absco(self, wn_out, p_mb, T, x_h2o):
        """Self+foreign continuum absorption coefficient [cm2/molec of H2O]
        on the wn_out grid, radiation term applied (radflag=.true.)."""
        rho_rat = (p_mb / self.p_ref) * (self.t_ref / T)
        sh2o = self.self_ref * (self.t_ref / T)**self.self_texp * x_h2o * rho_rat
        fh2o = self.for_ref * (1.0 - x_h2o) * rho_rat
        coeff = (sh2o + fh2o) * self._radfn(self.wvn, T)
        return np.interp(wn_out, self.wvn, coeff)


# ----------------------------------------------------------------------------
# Column handling
# ----------------------------------------------------------------------------
def load_column(path, nlay):
    """Read an ExoColumn output and coarsen to ~nlay layers (TOA->surface).
    Returns dict with per-layer T, p_mb, x_h2o, x_co2, N_h2o, N_co2 [molec/cm2]
    plus ts, and the n68 band data if present."""
    with nc.Dataset(path) as ds:
        tmid = np.array(ds['tmid'][:]); pint = np.array(ds['pint'][:])
        q = np.array(ds['h2ommr'][:]); qco2 = np.array(ds['co2mmr'][:])
        mwdry = float(ds['mw'][:]) * 1e-3            # kg/mol
        ts = float(ds['ts'][:])
        band = {}
        if 'band_lwup_toa' in ds.variables:
            band['olr'] = np.array(ds['band_lwup_toa'][:])
            band['edges'] = np.array(ds['wavenum_edge'][:])
    pdel = np.diff(pint)                              # [Pa], TOA->sfc
    Mtot = pdel / G                                   # kg/m2
    m_h2o = q * Mtot; M_dry = (1.0 - q) * Mtot
    m_co2 = qco2 * M_dry                              # qco2 is a DRY mmr
    N_h2o = m_h2o / MW_H2O * NA * 1e-4                # molec/cm2
    N_co2 = m_co2 / MW_CO2 * NA * 1e-4
    N_dry = M_dry / mwdry * NA * 1e-4
    N_tot = N_h2o + N_dry
    pmid = 0.5 * (pint[:-1] + pint[1:])

    # coarsen: contiguous chunks of equal count
    n = len(tmid)
    idx = np.array_split(np.arange(n), nlay)
    lay = dict(T=[], p_mb=[], x_h2o=[], x_co2=[], N_h2o=[], N_co2=[], N_dry=[])
    for ii in idx:
        w = pdel[ii] / pdel[ii].sum()
        lay['T'].append(float(np.sum(w * tmid[ii])))
        lay['p_mb'].append(float(np.sum(w * pmid[ii])) / 100.0)
        nh, nco, ntot = N_h2o[ii].sum(), N_co2[ii].sum(), N_tot[ii].sum()
        lay['N_h2o'].append(nh); lay['N_co2'].append(nco)
        lay['N_dry'].append(N_dry[ii].sum())
        lay['x_h2o'].append(nh / ntot); lay['x_co2'].append(nco / ntot)
    lay = {k: np.array(v) for k, v in lay.items()}
    lay['ts'] = ts; lay['band'] = band
    return lay


# ----------------------------------------------------------------------------
# RADIS optical depths
# ----------------------------------------------------------------------------
def layer_tau(lay, wmin, wmax, wstep):
    """tau(nlay, nwn) for H2O lines, CO2 lines (RADIS/HITRAN), and the MT_CKD
    H2O continuum, all on a common wavenumber grid."""
    from radis import SpectrumFactory
    nlay = len(lay['T'])
    wn = None
    tau_l = {}
    for mol, xkey, Nkey in (('H2O', 'x_h2o', 'N_h2o'), ('CO2', 'x_co2', 'N_co2')):
        sf = SpectrumFactory(wavenum_min=wmin, wavenum_max=wmax, wstep=wstep,
                             molecule=mol, isotope='1,2,3', pressure=1.0,
                             mole_fraction=1e-3, path_length=1.0,
                             verbose=0, warnings={'AccuracyError': 'ignore',
                                                  'AccuracyWarning': 'ignore'})
        sf.fetch_databank('hitran')
        tau = None
        for k in range(nlay):
            T = lay['T'][k]; p_mb = lay['p_mb'][k]
            x = max(lay[xkey][k], 1e-12); Nsp = lay[Nkey][k]
            s = sf.eq_spectrum(Tgas=T, pressure=p_mb / 1013.25, mole_fraction=x)
            w, kabs = s.get('abscoeff', wunit='cm-1')   # [cm-1] base-e
            if wn is None:
                wn = w;
            if tau is None:
                tau = np.zeros((nlay, len(w)))
            # path length s.t. x * n_tot * L = Nsp  ->  tau = kabs * L
            n_tot = (p_mb * 100.0) / (KB * T) * 1e-6    # molec/cm3
            L = Nsp / (x * n_tot)                       # cm
            tau[k] = kabs * L
            print(f"  [{mol}] layer {k+1:2d}/{nlay}  T={T:6.1f} K  p={p_mb:9.3f} mb"
                  f"  x={x:.3e}  max(tau)={tau[k].max():.2e}", flush=True)
        tau_l[mol] = tau
    # MT_CKD continuum
    ck = MTCKD()
    tau_c = np.zeros_like(tau_l['H2O'])
    for k in range(nlay):
        absco = ck.absco(wn, lay['p_mb'][k], lay['T'][k], lay['x_h2o'][k])
        tau_c[k] = absco * lay['N_h2o'][k]
    return wn, tau_l['H2O'], tau_l['CO2'], tau_c


def schwarzschild_olr(wn, tau_layers, T_layers, ts):
    """No-scattering OLR_nu with diffusivity D; layers TOA->surface."""
    D = DIFFUSIVITY
    # transmittance from TOA to the top of layer k
    tcum = np.exp(-D * np.cumsum(tau_layers, axis=0))
    t_above = np.vstack([np.ones_like(wn), tcum[:-1]])
    olr = planck_nu(wn, ts) * tcum[-1]               # surface term
    for k in range(len(T_layers)):
        olr += planck_nu(wn, T_layers[k]) * (t_above[k] - tcum[k])
    return olr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ncfile', nargs='?',
                    default=os.path.join(os.path.dirname(HERE), 'iofiles', 'exocol_out.nc'))
    ap.add_argument('--nlay', type=int, default=50)
    ap.add_argument('--wmin', type=float, default=10.0)
    ap.add_argument('--wmax', type=float, default=3000.0)
    ap.add_argument('--wstep', type=float, default=0.01)
    ap.add_argument('--out', default=os.path.join(
        os.path.dirname(HERE), 'reference', 'moist_runaway',
        'lbl_olr_benchmark_ts300'))
    args = ap.parse_args()

    lay = load_column(args.ncfile, args.nlay)
    print(f"column: Ts={lay['ts']:.1f} K, {len(lay['T'])} layers, "
          f"x_h2o sfc={lay['x_h2o'][-1]:.4f}, x_co2 sfc={lay['x_co2'][-1]:.3e}")

    wn, tau_h2o, tau_co2, tau_cont = layer_tau(lay, args.wmin, args.wmax, args.wstep)

    cases = {
        'lines only (H2O+CO2)': tau_h2o + tau_co2,
        'lines + MT_CKD continuum': tau_h2o + tau_co2 + tau_cont,
    }
    print(f"\n=== LBL OLR, Ts={lay['ts']:.0f} K column ===")
    results = {}
    for name, tau in cases.items():
        olr_nu = schwarzschild_olr(wn, tau, lay['T'], lay['ts'])
        olr = np.trapezoid(olr_nu, wn)
        # add the (tiny) Planck tail beyond [wmin, wmax]
        results[name] = (olr, olr_nu)
        print(f"  {name:28s}: OLR = {olr:7.2f} W/m2")

    # n68 band comparison if available
    band = lay['band']
    if band:
        edges = band['edges']
        print(f"\n  per-n68-band (W/m2): LBL(lines+cont) vs ExoRT band_lwup_toa")
        olr_nu = results['lines + MT_CKD continuum'][1]
        rows = []
        for i in range(len(band['olr'])):
            n1, n2 = edges[i], edges[i+1]
            if n2 < args.wmin or n1 > args.wmax:
                lbl_b = np.nan
            else:
                m = (wn >= n1) & (wn < n2)
                lbl_b = np.trapezoid(olr_nu[m], wn[m]) if m.sum() > 2 else 0.0
            rows.append((n1, n2, lbl_b, band['olr'][i]))
        for n1, n2, lbl_b, exo_b in rows:
            if np.isfinite(lbl_b) and (lbl_b > 0.5 or exo_b > 0.5):
                print(f"    {n1:7.0f}-{n2:7.0f}  LBL={lbl_b:7.2f}  ExoRT={exo_b:7.2f}"
                      f"  d={exo_b-lbl_b:+6.2f}")
        # fair totals over the common range only (bands fully inside [wmin,wmax])
        inr = [(l, e) for n1, n2, l, e in rows
               if np.isfinite(l) and n1 >= args.wmin and n2 <= args.wmax]
        lbl_t = sum(l for l, e in inr); exo_t = sum(e for l, e in inr)
        print(f"\n  common-range totals ({args.wmin:.0f}-{args.wmax:.0f} cm-1): "
              f"LBL = {lbl_t:.2f}  ExoRT n68 = {exo_t:.2f}  "
              f"(ExoRT - LBL = {exo_t-lbl_t:+.2f} W/m2)")
        np.savez_compressed(
            args.out + '.npz', wn=wn.astype(np.float32),
            olr_nu_lines=results['lines only (H2O+CO2)'][1].astype(np.float32),
            olr_nu_cont=results['lines + MT_CKD continuum'][1].astype(np.float32),
            band_edges=edges, band_exo=band['olr'],
            T=lay['T'], p_mb=lay['p_mb'], x_h2o=lay['x_h2o'], ts=lay['ts'])
        print(f"\nsaved {args.out}.npz")


if __name__ == '__main__':
    main()
