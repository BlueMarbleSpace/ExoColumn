#!/usr/bin/env python3
"""
hz_inner.py  —  Inner habitable zone OLR/Seff sweep (Kopparapu+2013 approach).

For each prescribed surface temperature Ts, ExoColumn builds a moist-adiabat
column (cold start: saturated troposphere, isothermal stratosphere at t_strato),
calls ExoRT radiation once (flux_only mode), and reads the fluxes.

Diagnostics
-----------
OLR  = LWUP[0]          (TOA upwelling LW, W/m²)
ASR  = SWDN[0]-SWUP[0]  (absorbed stellar radiation, W/m²)
Seff = OLR / ASR        (effective stellar flux relative to reference, dimensionless)

Seff > 1 means the column emits more than it absorbs at this Ts under the
reference stellar flux (msdist=1); the star must be brighter / closer.
The runaway greenhouse limit is where Seff(Ts) reaches its maximum.

Kopparapu+2013 Table 1 (G2V, T_eff=5780 K, 1 M_Earth):
  Moist greenhouse : Seff = 1.0129   (H2O stratospheric VMR > 3e-3)
  Runaway GH       : Seff = 1.0140   (OLR saturates, Simpson-Nakajima ~282 W/m²)

Usage
-----
  cd /hugespace/models/ExoColumn
  source /opt/intel/oneapi/setvars.sh
  cd build && make && cd ..
  python tools/hz_inner.py
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

# ---------------------------------------------------------------------------
# Sweep settings
# ---------------------------------------------------------------------------
TS_VALUES = np.arange(200, 425, 5, dtype=float)   # K; 200 to 420 inclusive

# Surface albedo: matches ExoColumn Earth reference case
ALBEDO = 0.2736

# Stratospheric temperature cap (Kopparapu-style isothermal stratosphere)
T_STRATO = 200.0   # K

# Kopparapu+2013 Table 1 inner-edge Seff (G2V star, 1 M_Earth) — approximate
KOPP_RUNAWAY = 1.0140
KOPP_MOIST   = 1.0129
SN_LIMIT     = 282.0   # W/m²  Simpson-Nakajima OLR saturation limit

# ---------------------------------------------------------------------------
# Namelist template
# ---------------------------------------------------------------------------
NML_TEMPLATE = """\
&exocol_nml
  flux_only  = .true.
  o3_profile = 'none'
  msdist     = 1.0
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
  n2_vmr  = 0.9996
  o2_vmr  = 0.0
  ch4_vmr = 0.0
  ar_vmr  = 0.0
/
"""

# ---------------------------------------------------------------------------

def run_one(ts):
    """
    Run ExoColumn in flux_only mode at surface temperature ts.
    Returns (olr, asr, alpha_p) in W/m², or None on failure.
    """
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO)

    # Save original nml (restored in finally block)
    orig = None
    if os.path.exists(NML_PATH):
        with open(NML_PATH) as f:
            orig = f.read()

    try:
        with open(NML_PATH, 'w') as f:
            f.write(nml)

        result = subprocess.run(
            [EXE], cwd=ROOT,
            capture_output=True, text=True, timeout=180
        )
        if result.returncode != 0:
            print(f"  FAIL Ts={ts:.0f}K  returncode={result.returncode}")
            print(result.stderr[-400:] if result.stderr else "(no stderr)")
            return None

        with nc.Dataset(OUT_NC) as ds:
            lwup = ds['LWUP'][:]   # pverp; index 0 = TOA (Fortran index 1)
            swup = ds['SWUP'][:]
            swdn = ds['SWDN'][:]

        olr   = float(lwup[0])
        asr   = float(swdn[0] - swup[0])
        alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
        return olr, asr, alpha

    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"ExoColumn executable not found: {EXE}\n"
                                "Run 'make' in the build/ directory first.")

    print("ExoColumn inner-HZ OLR sweep")
    print(f"  Ts range  : {TS_VALUES[0]:.0f}–{TS_VALUES[-1]:.0f} K  ({len(TS_VALUES)} points)")
    print(f"  t_strato  : {T_STRATO:.0f} K  (isothermal stratosphere)")
    print(f"  rh_init   : 1.0  (saturated troposphere)")
    print(f"  albedo    : {ALBEDO}")
    print(f"  mode      : flux_only (single radiation call per point)")
    print(f"  composition: N2 + CO2 (4e-4) + H2O only  [O2=O3=CH4=Ar=0, Kopparapu-style]")
    print()

    ts_list   = []
    olr_list  = []
    asr_list  = []
    seff_list = []
    alpha_list = []

    for ts in TS_VALUES:
        r = run_one(ts)
        if r is None:
            print(f"  Ts={ts:5.1f} K  →  skipped (run failed)")
            continue
        olr, asr, alpha = r
        seff = olr / asr if asr > 1e-6 else np.nan
        print(f"  Ts={ts:5.1f} K  OLR={olr:7.2f}  ASR={asr:7.2f}  α={alpha:.3f}  Seff={seff:.4f}")
        ts_list.append(ts)
        olr_list.append(olr)
        asr_list.append(asr)
        seff_list.append(seff)
        alpha_list.append(alpha)

    if not ts_list:
        print("No successful runs — nothing to plot.")
        return

    ts_arr    = np.array(ts_list)
    olr_arr   = np.array(olr_list)
    asr_arr   = np.array(asr_list)
    seff_arr  = np.array(seff_list)
    alpha_arr = np.array(alpha_list)

    _plot(ts_arr, olr_arr, asr_arr, seff_arr, alpha_arr)


def _plot(ts, olr, asr, seff, alpha):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 4.0), dpi=300)
    fig.patch.set_facecolor('white')

    # --- Panel 1: OLR and ASR vs Ts ---
    ax = axes[0]
    ax.plot(ts, olr, color='C3', lw=1.5, label='OLR')
    ax.plot(ts, asr, color='C0', lw=1.5, label='ASR (msdist=1)')
    ax.axhline(SN_LIMIT, color='k', lw=0.8, ls='--', label=f'S-N limit ({SN_LIMIT:.0f} W m⁻²)')
    ax.set_xlabel('$T_s$ (K)')
    ax.set_ylabel('Flux (W m⁻²)')
    ax.set_title('OLR and absorbed SW')
    ax.set_xlim(ts[0], ts[-1])
    ax.legend(fontsize=8, framealpha=0.9)
    ax.set_facecolor('white')

    # --- Panel 2: Seff vs Ts ---
    ax = axes[1]
    ax.plot(ts, seff, color='C1', lw=1.5, label='ExoColumn')
    ax.axhline(KOPP_RUNAWAY, color='C3', lw=0.9, ls='--',
               label=f'Kopparapu runaway ({KOPP_RUNAWAY:.4f})')
    ax.axhline(KOPP_MOIST, color='C0', lw=0.9, ls=':',
               label=f'Kopparapu moist GH ({KOPP_MOIST:.4f})')
    ax.axhline(1.0, color='gray', lw=0.6, ls='-', alpha=0.4)
    ax.set_xlabel('$T_s$ (K)')
    ax.set_ylabel('$S_{\\rm eff}$ (relative to modern solar)')
    ax.set_title('Effective stellar flux')
    ax.set_xlim(ts[0], ts[-1])
    ax.legend(fontsize=8, framealpha=0.9)
    ax.set_facecolor('white')

    # --- Panel 3: Planetary albedo vs Ts ---
    ax = axes[2]
    ax.plot(ts, alpha, color='C2', lw=1.5)
    ax.set_xlabel('$T_s$ (K)')
    ax.set_ylabel('Planetary albedo $\\alpha_p$')
    ax.set_title('Planetary albedo')
    ax.set_xlim(ts[0], ts[-1])
    ax.set_facecolor('white')

    fig.tight_layout()
    out_dir = os.path.dirname(os.path.abspath(__file__))
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(out_dir, f'hz_inner.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
