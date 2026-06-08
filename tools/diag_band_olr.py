#!/usr/bin/env python3
"""
diag_band_olr.py  —  Band-resolved TOA OLR diagnostic for the inner-HZ runaway
greenhouse, built to answer one question:

    As Ts climbs into the runaway regime, does the H2O continuum close the
    8-12 um atmospheric window (so OLR caps at the Simpson-Nakajima limit), or
    does the window leak (so OLR keeps climbing with Ts)?

See memory project_steam_runaway_olr: with the correct Kasting (1988) non-ideal
steam adiabat the nonideal sweep shows OLR climbing past the S-N limit instead
of capping.  Kasting (1988, p.474) attributes the runaway cap to the *absence
of IR windows* — the H2O continuum makes a dense steam atmosphere optically
thick at ALL infrared wavelengths.  ExoRT wires the MT_CKD 3.3 self+foreign
continuum into the n68 path (calc_opd_mod.F90:323-328), so this diagnostic
checks, band by band, whether that window actually closes.

ExoColumn now writes the band-resolved TOA OLR (band_lwup_toa, W/m2 per n68
band) and the band edges (wavenum_edge, cm-1) into iofiles/exocol_out.nc.

For a sweep of prescribed Ts (flux_only mode, same column construction as
tools/hz_inner.py), this script reads band_lwup_toa and produces:

  (a) Spectral OLR density (W m-2 / cm-1) vs wavenumber, one curve per Ts,
      with the 8-12 um window shaded.  Planck flux pi*B at a few T overlaid.
  (b) Brightness temperature T_b(nu) vs wavenumber, one curve per Ts.
      If the window closes, T_b in the window tracks the cold tropopause
      (~t_strato); if it leaks, T_b in the window tracks the (hot) surface.

A summary table reports total OLR, in-window OLR, and out-of-window OLR vs Ts
so the leak is quantified.

Usage:
    python tools/diag_band_olr.py
Env overrides:
    HZ_H2O_EOS   ideal | nonideal   (default nonideal)
    BAND_TS      comma-separated Ts list, e.g. "300,400,500,600"
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

# ---------------------------------------------------------------------------
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

H2O_EOS = os.environ.get('HZ_H2O_EOS', 'nonideal')
TAG     = '' if H2O_EOS == 'ideal' else '_' + H2O_EOS

_default_ts = "300,450,600,700,1000,1300,1500"   # spans sub- & super-critical (Tc=647 K)
TS_VALUES = np.array([float(x) for x in
                      os.environ.get('BAND_TS', _default_ts).split(',')])

ALBEDO   = 0.24229
T_STRATO = 200.0
SN_LIMIT = 282.0       # W/m2  Simpson-Nakajima OLR limit

# 8-12 um atmospheric window in wavenumber [cm-1]
WIN_LO, WIN_HI = 1.0e4 / 12.0, 1.0e4 / 8.0     # 833.3 .. 1250 cm-1

# Physical constants (SI)
_H = 6.62607015e-34    # J s
_C = 2.99792458e8      # m/s
_KB = 1.380649e-23     # J/K
_SIGMA = 5.670374419e-8

NML_TEMPLATE = """\
&exocol_nml
  flux_only   = .true.
  variable_ps = .true.
  ihz_profile = .true.
  o3_profile  = 'none'
  msdist      = 1.0
  h2o_eos     = '{h2o_eos}'
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = {t_strato:.1f}
  p_top      = 1.0
  rh_init    = 1.0
  coszrs     = 0.5
  asdir      = {albedo:.4f}
  asdif      = {albedo:.4f}
  aldir      = {albedo:.4f}
  aldif      = {albedo:.4f}
/
&exocol_composition
  n2_vmr   = 0.78
  o2_vmr   = 0.210
  ar_vmr   = 0.01
  co2_vmr  = 3.3e-4
  ch4_vmr  = 0.0
  o3_vmr   = 0.0
/
"""


# ---------------------------------------------------------------------------
def planck_wn(wn_cm, T):
    """Planck spectral radiance B(nu,T) in W/(m^2 sr cm^-1).
    wn_cm: wavenumber [cm^-1] (scalar or array).  T: temperature [K]."""
    wn_m = np.asarray(wn_cm, dtype=float) * 100.0      # cm^-1 -> m^-1
    x = _H * _C * wn_m / (_KB * T)
    # radiance per unit (m^-1): 2 h c^2 nu^3 / (exp(x)-1)
    B_per_m = 2.0 * _H * _C**2 * wn_m**3 / np.expm1(x)
    return B_per_m * 100.0                              # per m^-1 -> per cm^-1


def brightness_temp(F_band, dwn, wn_mid):
    """Invert band-integrated upward flux F_band [W/m2] to a brightness
    temperature, treating the band as monochromatic at wn_mid with width dwn.
    F_band = pi * B(wn_mid, T_b) * dwn  ->  solve for T_b."""
    B = F_band / (np.pi * dwn)                          # W/(m^2 sr cm^-1)
    wn_m = wn_mid * 100.0
    c1 = 2.0 * _H * _C**2 * wn_m**3 * 100.0             # same scaling as planck_wn
    c2 = _H * _C * wn_m / _KB
    Tb = np.full_like(B, np.nan)
    good = B > 0
    Tb[good] = c2[good] / np.log1p(c1[good] / B[good])
    return Tb


def run_one(ts):
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO,
                              h2o_eos=H2O_EOS)
    orig = None
    if os.path.exists(NML_PATH):
        with open(NML_PATH) as f:
            orig = f.read()
    try:
        with open(NML_PATH, 'w') as f:
            f.write(nml)
        result = subprocess.run([EXE], cwd=ROOT, capture_output=True,
                                text=True, timeout=180)
        if result.returncode != 0:
            print(f"  FAIL Ts={ts:.0f} K  rc={result.returncode}")
            print(result.stderr[-300:] if result.stderr else '')
            return None
        with nc.Dataset(OUT_NC) as ds:
            band_olr = ds['band_lwup_toa'][:]
            wn_edge  = ds['wavenum_edge'][:]
            lwup     = ds['LWUP'][:]
        return dict(band_olr=np.array(band_olr), wn_edge=np.array(wn_edge),
                    olr_tot=float(lwup[0]))
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"Executable not found: {EXE}")

    # Planck normalisation self-check: integrate pi*B over a fine wn grid at
    # 288 K and compare to sigma T^4 (should agree to <0.5%).
    wn_chk = np.linspace(1.0, 6000.0, 60000)
    _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    flux_chk = np.pi * _trapz(planck_wn(wn_chk, 288.0), wn_chk)
    print(f"Planck check: pi*int(B) @288K = {flux_chk:7.2f} W/m2  "
          f"(sigma T^4 = {_SIGMA*288.0**4:7.2f})")

    print(f"\nBand-resolved OLR diagnostic  —  H2O EOS = {H2O_EOS}")
    print(f"  Ts values: {', '.join(f'{t:.0f}' for t in TS_VALUES)} K\n")

    results = {}
    wn_edge = None
    for ts in TS_VALUES:
        r = run_one(ts)
        if r is None:
            print(f"  Ts={ts:5.1f} K  SKIPPED")
            continue
        results[ts] = r
        wn_edge = r['wn_edge']

    if not results:
        print("No successful runs.")
        return

    wn_mid = 0.5 * (wn_edge[:-1] + wn_edge[1:])
    dwn    = np.diff(wn_edge)
    win    = (wn_mid >= WIN_LO) & (wn_mid <= WIN_HI)   # 8-12 um window
    far_ir = wn_mid < WIN_LO                            # rotational far-IR
    near_ir = wn_mid > WIN_HI                           # near-IR + SW-IR

    # OLR split into far-IR / window / near-IR shows WHERE the runaway flux
    # escapes: if the window closes but OLR climbs, the leak is elsewhere.
    print(f"\n{'Ts':>6} {'OLR_tot':>8} | {'farIR':>7} {'window':>7} "
          f"{'nearIR':>7} | {'Tb_far':>7} {'Tb_win':>7} {'Tb_nir':>7}")
    print(f"{'(K)':>6} {'(W/m2)':>8} | {'<833':>7} {'8-12um':>7} "
          f"{'>1250':>7} | {'(K)':>7} {'(K)':>7} {'(K)':>7}")
    for ts in sorted(results):
        bo = results[ts]['band_olr']
        Tb = brightness_temp(bo, dwn, wn_mid)
        print(f"{ts:6.0f} {results[ts]['olr_tot']:8.2f} | "
              f"{bo[far_ir].sum():7.2f} {bo[win].sum():7.2f} "
              f"{bo[near_ir].sum():7.2f} | "
              f"{np.nanmean(Tb[far_ir]):7.1f} {np.nanmean(Tb[win]):7.1f} "
              f"{np.nanmean(Tb[near_ir]):7.1f}")

    _plot(results, wn_edge, wn_mid, dwn, win)


def _plot(results, wn_edge, wn_mid, dwn, win):
    ts_list = sorted(results)
    norm = Normalize(vmin=min(ts_list), vmax=max(ts_list))
    cmap = plt.colormaps['inferno']

    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(7.0, 7.5), dpi=300)
    fig.patch.set_facecolor('white')

    # --- (a) spectral OLR density ---
    for ts in ts_list:
        bo = results[ts]['band_olr']
        spec = bo / dwn                       # W/m2/cm-1
        ax_a.plot(wn_mid, spec, color=cmap(norm(ts)), lw=1.3)
    # Planck flux pi*B reference curves (dashed grey) at a low and high T
    for Tref, ls in [(T_STRATO, ':'), (max(ts_list), '--')]:
        ax_a.plot(wn_mid, np.pi * planck_wn(wn_mid, Tref), color='0.4',
                  lw=0.8, ls=ls, label=f'$\\pi B$ at {Tref:.0f} K')
    ax_a.axvspan(WIN_LO, WIN_HI, color='C0', alpha=0.10)
    ax_a.set_xlim(0, 2500)
    ax_a.set_ylabel('Spectral OLR (W m$^{-2}$ / cm$^{-1}$)')
    ax_a.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax_a.legend(fontsize=7, framealpha=0.9, loc='upper right')
    ax_a.set_facecolor('white')
    ax_a.text(0.02, 0.96, '(a)', transform=ax_a.transAxes, va='top',
              fontsize=9, fontweight='bold')
    ax_a.text(0.5 * (WIN_LO + WIN_HI), ax_a.get_ylim()[1] * 0.92,
              '8–12 µm\nwindow', ha='center', va='top', fontsize=7, color='C0')

    # --- (b) brightness temperature spectrum ---
    for ts in ts_list:
        bo = results[ts]['band_olr']
        Tb = brightness_temp(bo, dwn, wn_mid)
        ax_b.plot(wn_mid, Tb, color=cmap(norm(ts)), lw=1.3)
    ax_b.axhline(T_STRATO, color='k', lw=0.8, ls=':',
                 label=f'$t_{{strato}}$ = {T_STRATO:.0f} K')
    ax_b.axvspan(WIN_LO, WIN_HI, color='C0', alpha=0.10)
    ax_b.set_xlim(0, 2500)
    ax_b.set_ylabel('Brightness temperature $T_b$ (K)')
    ax_b.set_xlabel('Wavenumber (cm$^{-1}$)')
    ax_b.legend(fontsize=7, framealpha=0.9, loc='upper right')
    ax_b.set_facecolor('white')
    ax_b.text(0.02, 0.96, '(b)', transform=ax_b.transAxes, va='top',
              fontsize=9, fontweight='bold')

    # shared colorbar for Ts
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax_a, ax_b], fraction=0.04, pad=0.02)
    cbar.set_label('$T_s$ (K)')

    out_dir = os.path.dirname(os.path.abspath(__file__))
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(out_dir, f'diag_band_olr{TAG}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"\nSaved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
