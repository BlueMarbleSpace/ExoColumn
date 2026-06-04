#!/usr/bin/env python3
"""
hz_inner.py  —  ExoColumn version of Kopparapu+2013 Figure 3 (inner HZ, G2V star).

For each prescribed surface temperature Ts, ExoColumn builds a moist-adiabat
column (cold start: fully saturated at all levels including stratosphere,
isothermal stratosphere at t_strato = 200 K), calls ExoRT radiation once
(flux_only mode), and reads the output profiles and fluxes.

Four-panel figure matching Kopparapu+2013 Fig. 3:
  (a) OLR and absorbed SW vs Ts
  (b) Planetary albedo vs Ts
  (c) Effective stellar flux Seff = OLR / ASR vs Ts
  (d) H2O VMR vertical profiles for selected Ts values

Composition: N2=0.78 + O2=0.21 + Ar=0.01 + CO2=3.3e-4 + H2O (Kopparapu 2013 Earth tuning; no O3/CH4).
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

# ---------------------------------------------------------------------------
TS_VALUES = np.arange(200, 425, 5, dtype=float)   # K

ALBEDO   = 0.24229   # surface albedo calibrated to Ts=288 K (Kopparapu composition, fixed cold-trap)
T_STRATO = 200.0     # isothermal stratosphere cap [K]

MW_H2O   = 18.015    # g/mol

# Kopparapu+2013 Table 1 inner-edge Seff (G2V, 1 M_Earth)
KOPP_RUNAWAY = 1.0140
KOPP_MOIST   = 1.0129
SN_LIMIT     = 282.0    # W/m²  Simpson-Nakajima OLR limit
MOIST_GH_VMR = 3.0e-3   # moist greenhouse stratospheric H2O VMR threshold

# Ts values at which to save the full H2O vertical profile for panel (d)
PROFILE_TS = [250.0, 280.0, 310.0, 340.0]

NML_TEMPLATE = """\
&exocol_nml
  flux_only   = .true.
  variable_ps = .true.
  o3_profile  = 'none'
  msdist      = 1.0
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

def run_one(ts):
    """
    Run ExoColumn flux_only at surface temperature ts.
    Returns dict with scalar diagnostics and profiles, or None on failure.
    """
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO)
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
            lwup   = ds['LWUP'][:]
            swup   = ds['SWUP'][:]
            swdn   = ds['SWDN'][:]
            tmid   = ds['tmid'][:]
            pmid   = ds['pmid'][:]
            h2ommr = ds['h2ommr'][:]
            mwdry  = float(ds['mw'][:])

        olr   = float(lwup[0])
        asr   = float(swdn[0] - swup[0])
        alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
        seff  = olr / asr if asr > 1e-6 else np.nan

        # H2O volume mixing ratio profile: VMR ≈ h2ommr × mwdry / MW_H2O
        h2o_vmr = h2ommr * mwdry / MW_H2O

        return dict(olr=olr, asr=asr, alpha=alpha, seff=seff,
                    h2o_vmr=h2o_vmr, pmid=np.array(pmid),
                    tmid=np.array(tmid))
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"Executable not found: {EXE}")

    print("ExoColumn inner-HZ sweep  —  Figure 3 (Kopparapu+2013 style)")
    print(f"  Ts range    : {TS_VALUES[0]:.0f}–{TS_VALUES[-1]:.0f} K")
    print(f"  t_strato    : {T_STRATO:.0f} K  (isothermal, fully saturated)")
    print(f"  variable_ps : ps = p_N2(1 bar) + esat(Ts)")
    print(f"  composition : N2=0.78 + O2=0.21 + Ar=0.01 + CO2=3.3e-4 + H2O  [no O3/CH4]")
    print()

    rows = []
    profiles = {}   # ts → (pmid, h2o_vmr) for selected Ts
    for ts in TS_VALUES:
        r = run_one(ts)
        if r is None:
            print(f"  Ts={ts:5.1f} K  SKIPPED")
            continue
        print(f"  Ts={ts:5.1f} K  OLR={r['olr']:7.2f}  ASR={r['asr']:7.2f}"
              f"  α={r['alpha']:.3f}  Seff={r['seff']:.4f}")
        rows.append((ts, r['olr'], r['asr'], r['alpha'], r['seff']))
        if any(abs(ts - ts_p) < 0.5 for ts_p in PROFILE_TS):
            profiles[ts] = (r['pmid'], r['h2o_vmr'])

    if not rows:
        print("No successful runs.")
        return

    arr  = np.array(rows)
    ts_a, olr_a, asr_a, alp_a, seff_a = arr.T

    _plot(ts_a, olr_a, asr_a, alp_a, seff_a, profiles)


def _plot(ts, olr, asr, alpha, seff, profiles):
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), dpi=300)
    fig.patch.set_facecolor('white')
    (ax_a, ax_b), (ax_c, ax_d) = axes

    kw = dict(xlim=(ts[0], ts[-1]))

    # --- (a) OLR and ASR vs Ts ---
    ax_a.plot(ts, olr, color='C3', lw=1.5, label='OLR')
    ax_a.plot(ts, asr, color='C0', lw=1.5, label='Absorbed SW')
    ax_a.axhline(SN_LIMIT, color='k', lw=0.8, ls='--',
                 label=f'S-N limit ({SN_LIMIT:.0f} W m⁻²)')
    ax_a.set_ylabel('Flux (W m⁻²)')
    ax_a.set_xlabel('$T_s$ (K)')
    ax_a.legend(fontsize=7, framealpha=0.9)
    ax_a.set(**kw)
    ax_a.set_facecolor('white')
    ax_a.text(0.04, 0.96, '(a)', transform=ax_a.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (b) Planetary albedo vs Ts ---
    ax_b.plot(ts, alpha, color='C2', lw=1.5)
    ax_b.set_ylabel('Planetary albedo $\\alpha_p$')
    ax_b.set_xlabel('$T_s$ (K)')
    ax_b.set(**kw)
    ax_b.set_facecolor('white')
    ax_b.text(0.04, 0.96, '(b)', transform=ax_b.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (c) Seff vs Ts ---
    ax_c.plot(ts, seff, color='C1', lw=1.5, label='ExoColumn (G2V)')
    ax_c.axhline(KOPP_RUNAWAY, color='C3', lw=0.9, ls='--',
                 label=f'Kopparapu runaway ({KOPP_RUNAWAY:.4f})')
    ax_c.axhline(KOPP_MOIST, color='C0', lw=0.9, ls=':',
                 label=f'Kopparapu moist GH ({KOPP_MOIST:.4f})')
    ax_c.axhline(1.0, color='gray', lw=0.5, alpha=0.5)
    ax_c.set_ylabel('$S_{\\rm eff}$')
    ax_c.set_xlabel('$T_s$ (K)')
    ax_c.legend(fontsize=7, framealpha=0.9)
    ax_c.set(**kw)
    ax_c.set_facecolor('white')
    ax_c.text(0.04, 0.96, '(c)', transform=ax_c.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (d) H2O VMR vertical profiles for selected Ts values ---
    colors_d = ['C0', 'C2', 'C1', 'C3']
    for (ts_p, (pmid, vmr)), col in zip(sorted(profiles.items()), colors_d):
        # Mask out any unphysically large VMR near p_top (TOA, p < 10 Pa)
        mask = pmid > 5.0
        ax_d.loglog(vmr[mask], pmid[mask] / 100., color=col, lw=1.5,
                    label=f'$T_s$ = {ts_p:.0f} K')
    ax_d.axvline(MOIST_GH_VMR, color='k', lw=0.8, ls='--',
                 label=f'Moist GH threshold')
    ax_d.invert_yaxis()
    ax_d.set_xlabel('H₂O VMR')
    ax_d.set_ylabel('Pressure (hPa)')
    ax_d.legend(fontsize=7, framealpha=0.9)
    ax_d.set_facecolor('white')
    ax_d.text(0.04, 0.04, '(d)', transform=ax_d.transAxes,
              va='bottom', fontsize=9, fontweight='bold')

    fig.tight_layout()
    out_dir = os.path.dirname(os.path.abspath(__file__))
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(out_dir, f'hz_inner.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
