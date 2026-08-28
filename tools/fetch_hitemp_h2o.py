#!/usr/bin/env python3
"""
fetch_hitemp_h2o.py — compute HITEMP-2010 (and HITRAN) H2O shortwave cross-sections
with RADIS and band-average them onto the ExoRT n68 grid, for the n68-vs-HITEMP
absorption comparison (plot_h2o_sw_absorption.py).

WHY this is a separate, user-run step
-------------------------------------
HITEMP downloads are gated behind a free HITRAN.org account.  RADIS reads the
credentials from the environment variables HITRAN_EMAIL and HITRAN_PASSWORD (or
prompts interactively and caches them encrypted in ~/radis.json on first use).
Run this ONCE, in your own terminal, with your credentials available, e.g.:

    HITRAN_EMAIL='you@example.com' HITRAN_PASSWORD='...' python3 tools/fetch_hitemp_h2o.py

(or just `python3 tools/fetch_hitemp_h2o.py` and answer the getpass prompt).
The first run downloads the HITEMP H2O line list for the SW range and may take a
few minutes + a few hundred MB; subsequent runs use the RADIS cache.

Output: tools/h2o_radis_bands.npz  (band-averaged sigma on the n68 grid), which
plot_h2o_sw_absorption.py overlays.  No credentials are needed to plot afterwards.
"""
import os
import sys
import numpy as np
import netCDF4 as nc

# ExoRT tree.  Override with the EXORT_ROOT environment variable (the same
# name the build uses in config.mk) if ExoRT lives elsewhere.
EXORT = os.environ.get('EXORT_ROOT', '/models/ExoRT')
CFILE = os.path.join(EXORT, 'data/continuum/KH2O_MTCKD3.3_SELF.FRGN_n68_ngauss.nc')
OUT   = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'h2o_radis_bands.npz')

# Shortwave range to cover the near-IR H2O bands (1-5 um: 2.7, 1.87, 1.38, 1.13,
# 0.94 um) where the solar near-IR energy and the HITEMP-vs-HITRAN difference live.
WN_MIN, WN_MAX = 2000.0, 12000.0      # cm-1  (0.83 - 5.0 um)
WSTEP = 0.01                           # cm-1 (resolve ~0.02 cm-1 lines at 0.1 bar)
PRESSURE_BAR = 0.1                     # 100 hPa
XMOL = 0.1                             # H2O mole fraction (broadening; sigma per molec)
TEMPS = [300.0, 500.0]                 # K (500 = ExoRT n68 k-table ceiling)


def _patch_radis_register():
    """Work around a RADIS 0.17 bug (hitempapi.py:1372): for the on-demand H2O/CO2
    HITEMP databases, HITEMPDatabaseManager.wmin/wmax are None, and register() formats
    them into an info string -> TypeError.  The line list still downloads fine; only
    this cosmetic metadata crashes.  Fill the Nones before register (runtime patch in
    OUR code; the installed package source is left untouched)."""
    try:
        from radis.api.hitempapi import HITEMPDatabaseManager as M
    except Exception:
        return
    if getattr(M.register, '_wminmax_guard', False):
        return
    _orig = M.register
    def register(self, download):
        if getattr(self, 'wmin', None) is None:
            self.wmin = 0.0
        if getattr(self, 'wmax', None) is None:
            self.wmax = 30000.0          # HITEMP-2010 H2O spans ~0-30000 cm-1
        return _orig(self, download)
    register._wminmax_guard = True
    M.register = register


def n68_band_edges():
    cds = nc.Dataset(CFILE)
    return np.array(cds['NU_LOW'][:]), np.array(cds['NU_HIGH'][:])


def band_average(wn, sig, nu_lo, nu_hi):
    """Linear spectral mean of sigma(wn) over each n68 band [nu_lo,nu_hi].
    This matches the correlated-k band-mean (sum_g g_w k_g = <k> over the band)."""
    out = np.full(nu_lo.shape, np.nan)
    for i, (lo, hi) in enumerate(zip(nu_lo, nu_hi)):
        m = (wn >= lo) & (wn < hi)
        if m.any():
            out[i] = np.nanmean(sig[m])
    return out


def radis_sigma(databank, T):
    from radis import calc_spectrum
    s = calc_spectrum(WN_MIN, WN_MAX, molecule='H2O', isotope='1,2,3',
                      Tgas=T, pressure=PRESSURE_BAR, mole_fraction=XMOL,
                      databank=databank, wstep=WSTEP, verbose=False,
                      warnings={'AccuracyError': 'ignore',
                                'GaussianBroadeningWarning': 'ignore'})
    try:
        wn, sig = s.get('xsection', wunit='cm-1')        # cm2/molecule
    except Exception:
        wn, ab = s.get('abscoeff', wunit='cm-1')          # cm-1
        kB = 1.380649e-23
        n = XMOL * (PRESSURE_BAR * 1e5) / (kB * T)         # m-3
        sig = ab / (n * 1e-6)                              # / cm-3  -> cm2
    return np.asarray(wn), np.asarray(sig)


def main():
    _patch_radis_register()
    nu_lo, nu_hi = n68_band_edges()
    res = {'nu_lo': nu_lo, 'nu_hi': nu_hi}
    for db in ('hitemp', 'hitran'):
        for T in TEMPS:
            print(f"  computing {db} H2O sigma at {T:.0f} K "
                  f"({WN_MIN:.0f}-{WN_MAX:.0f} cm-1) ...", flush=True)
            try:
                wn, sig = radis_sigma(db, T)
            except Exception as e:
                print(f"  *** {db}@{T:.0f}K failed: {e}", file=sys.stderr)
                if db == 'hitemp':
                    print("  (HITEMP needs HITRAN.org credentials — see header.)",
                          file=sys.stderr)
                continue
            res[f'sig_{db}_{int(T)}'] = band_average(wn, sig, nu_lo, nu_hi)
            print(f"    done ({len(wn)} pts).", flush=True)
    np.savez(OUT, **res)
    print(f"Saved: {OUT}\n  keys: {sorted(res.keys())}")


if __name__ == '__main__':
    main()
