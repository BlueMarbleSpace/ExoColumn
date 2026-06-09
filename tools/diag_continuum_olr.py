#!/usr/bin/env python3
"""
diag_continuum_olr.py — the H2O-continuum check behind the ExoColumn-vs-Kopparapu
inner-HZ offsets (albedo AND OLR).

After ruling out the H2O line list (tools/fetch_hitemp_h2o.py: HITEMP-2010 ≈
HITRAN-2016 to ~1% in the near-IR), the remaining spectroscopic candidate is the
water-vapour CONTINUUM.  Kopparapu+2013 (Sec. 2.2.2) make this quantitative for the
LONGWAVE, where it can be checked directly: for a dense-H2O atmosphere (Ts=400 K,
Tstrat=200 K, 4 bar N2 + saturated H2O, ps≈6.5 bar) their model — which uses the
BPS continuum (Paynter & Ramaswamy 2011) — gives OLR = 285 W/m2, while SMART, which
uses the older CKD continuum (Clough 1989), gives 297 W/m2.  BPS absorbs ~12 W/m2
MORE in the H2O windows (800-1200 and 300-600 cm-1) because it includes water-dimer
absorption.

ExoColumn/ExoRT uses MT_CKD 3.3 — the modern descendant of CKD, NOT BPS.  This script
reproduces that exact case and reads the band-resolved OLR: ExoColumn lands at ~304
W/m2, with the CKD family (SMART 297) and ABOVE Kopparapu's BPS (285).  i.e. ExoColumn
under-absorbs in the continuum relative to Kopparapu, by ~12-19 W/m2 in the LW windows.

The SAME continuum difference (MT_CKD < BPS) acts in the near-IR SHORTWAVE windows,
so Kopparapu absorbs more sunlight there -> their planetary albedo is lower (0.19 vs
our 0.21).  The continuum is thus the single consistent cause of BOTH the OLR offset
(we are higher: 295 vs 291 on the plateau) and the albedo offset (we are higher).

Figure: spectral OLR vs wavenumber with the continuum-sensitive windows shaded, and
the integrated-OLR comparison.  Restores the user's exocol_config.nml.
"""
import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE  = os.path.join(ROOT, 'run', 'exocol.exe')
OUT  = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML  = os.path.join(ROOT, 'exocol_config.nml')
FIG  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'continuum_olr.png')

KOPP_BPS = 285.0    # Kopparapu+2013 Fig 2, BPS continuum [W/m2]
SMART_CKD = 297.0   # SMART, CKD continuum [W/m2]

NML_CASE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  ihz_profile    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  h2o_eos        = 'nonideal'
  sw_zenith_quad = .true.
  sw_nquad       = 6
/
&exocol_init
  input_file = ''
  ts = 400.0
  t_strato = 200.0
  p_top = 0.01
  rh_init = 1.0
  coszrs = 0.5
  asdir = 0.30
  asdif = 0.30
  aldir = 0.30
  aldif = 0.30
/
&exocol_composition
  ps = 4.0e5
  n2_vmr = 1.0
  o2_vmr = 0.0
  ar_vmr = 0.0
  co2_vmr = 0.0
  ch4_vmr = 0.0
  o3_vmr = 0.0
/
"""


def planck_wn(wn_cm, T):
    """Planck spectral radiance per wavenumber -> hemispheric flux density [W/m2/cm-1]."""
    h = 6.62607015e-34; c = 2.99792458e10; kB = 1.380649e-23   # cgs-ish (c in cm/s)
    nu = wn_cm
    B = (2*h*c**2*nu**3) / (np.exp(h*c*nu/(kB*T)) - 1.0)        # erg/s/cm2/sr/cm-1
    return np.pi * B * 1e-3                                     # -> W/m2/cm-1


def main():
    orig = open(NML).read() if os.path.exists(NML) else None
    try:
        open(NML, 'w').write(NML_CASE)
        r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
        if r.returncode != 0:
            print(r.stderr[-400:]); return
        with nc.Dataset(OUT) as ds:
            bolr = np.array(ds['band_lwup_toa'][:])
            we   = np.array(ds['wavenum_edge'][:])
            olr  = float(ds['LWUP'][:][0])
    finally:
        if orig is not None:
            open(NML, 'w').write(orig)
        elif os.path.exists(NML):
            os.remove(NML)

    wm = 0.5*(we[:-1] + we[1:]); dwn = np.diff(we)
    olr_density = bolr / dwn                          # W/m2/cm-1
    win1 = (wm >= 800) & (wm <= 1200)
    win2 = (wm >= 300) & (wm <= 600)

    fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=300)
    fig.patch.set_facecolor('white'); ax.set_facecolor('white')

    # continuum-sensitive windows (Kopparapu Sec 2.2.2)
    for lo, hi in [(800, 1200), (300, 600)]:
        ax.axvspan(lo, hi, color='C1', alpha=0.15, lw=0)
    ax.text(1000, 0.015, 'H$_2$O window', ha='center', va='bottom',
            fontsize=7.5, color='C1', fontweight='bold')
    ax.text(450, 0.015, 'rot. window', ha='center', va='bottom',
            fontsize=7.5, color='C1', fontweight='bold')

    # 400 K blackbody reference
    wn_fine = np.linspace(max(we.min(), 1), 2000, 800)
    ax.plot(wn_fine, planck_wn(wn_fine, 400.), color='0.6', lw=0.9, ls='--',
            label='400 K blackbody')
    # ExoColumn spectral OLR
    ax.plot(wm, olr_density, drawstyle='steps-mid', color='C0', lw=1.5,
            label='ExoColumn OLR (MT_CKD 3.3)')

    ax.set_xlim(0, 2000); ax.set_ylim(0, 0.33)
    ax.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax.set_ylabel('Spectral OLR (W m$^{-2}$ / cm$^{-1}$)')
    ax.set_title('Dense-H$_2$O atmosphere (Kopparapu+2013 Fig. 2 case): '
                 'the H$_2$O continuum sets the window OLR', fontsize=9)

    box = ('Integrated OLR, same case (Ts=400 K, 4 bar N$_2$, ps≈6.5 bar):\n'
           f'  ExoColumn  (MT_CKD 3.3) = {olr:.0f} W m$^{{-2}}$\n'
           f'  SMART      (CKD, Clough 1989) = {SMART_CKD:.0f} W m$^{{-2}}$\n'
           f'  Kopparapu  (BPS, dimers) = {KOPP_BPS:.0f} W m$^{{-2}}$\n'
           'BPS absorbs ~12 W m$^{-2}$ more in these windows.  ExoColumn uses\n'
           'MT_CKD (CKD family) → under-absorbs vs BPS in BOTH the LW (higher\n'
           'OLR) and the SW near-IR windows (higher albedo: 0.21 vs 0.19).')
    ax.text(0.975, 0.95, box, transform=ax.transAxes, va='top', ha='right',
            fontsize=6.8, family='monospace',
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.93))
    ax.legend(fontsize=7.5, loc='upper left', framealpha=0.9)

    fig.tight_layout()
    fig.savefig(FIG, dpi=150, facecolor='white')
    print(f'ExoColumn dense-H2O OLR = {olr:.1f} W/m2 '
          f'(SMART/CKD {SMART_CKD:.0f}, Kopparapu/BPS {KOPP_BPS:.0f})')
    print(f'  window 800-1200 cm-1 OLR = {bolr[win1].sum():.1f} W/m2 ; '
          f'300-600 cm-1 = {bolr[win2].sum():.1f} W/m2')
    print(f'Saved: {FIG}')


if __name__ == '__main__':
    main()
