#!/usr/bin/env python3
"""
Calibrate ExoColumn surface albedo for the Kopparapu Earth composition.

Composition (from Kopparapu, communicated 2026-06-04):
  Ar   = 1.0e-2
  CH4  = ~0 (1e-60)
  CO2  = 3.3e-4
  N2   = 0.78
  O2   = 0.210
  O3   = none

Iterates albedo until converged Ts = 288.0 K (± 0.02 K).
Saves the calibrated config and copies output to this directory.
"""

import os, subprocess, shutil
import netCDF4 as nc
import numpy as np

ROOT     = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')
REF_DIR  = os.path.dirname(os.path.abspath(__file__))

TS_TARGET   = 288.0   # K
TS_TOL      = 0.02    # K
ALPHA_SENS  = -120.0  # K per unit albedo (local slope near 0.27-0.32)
ALPHA_INIT  =  0.305  # first guess: higher than 0.2736 since no O3

NML_TEMPLATE = """\
&exocol_nml
  conv_scheme     = 'sbm'
  tau_conv        = 3600.0
  rh_sbm          = 0.7
  moisture_scheme = 'prognostic'
  o3_profile      = 'none'
  wind_speed      = 5.0
  C_D             = 1.5e-3
  msdist          = 1.0
  dz_slab         = 10.0
  surface_flux    = 'mos'
  z0_rough        = 3.21e-5
  latent_heat_mode = 'phase_aware'
/
&exocol_init
  input_file = ''
  ts         = 288.0
  t_strato   = 200.0
  p_top      = 100.0
  rh_init    = 0.7
  coszrs     = 0.5
  asdir      = {alpha:.6f}
  asdif      = {alpha:.6f}
  aldir      = {alpha:.6f}
  aldif      = {alpha:.6f}
/
&exocol_composition
  ps       = 1.0e5
  co2_vmr  = 3.3e-4
  ch4_vmr  = 0.0
  n2_vmr   = 0.78
  o2_vmr   = 0.210
  ar_vmr   = 0.01
  o3_vmr   = 0.0
/
"""

def run_one(alpha):
    """Write namelist, run ExoColumn, return converged Ts."""
    # save and restore any existing namelist
    orig = None
    if os.path.exists(NML_PATH):
        with open(NML_PATH) as f:
            orig = f.read()
    try:
        with open(NML_PATH, 'w') as f:
            f.write(NML_TEMPLATE.format(alpha=alpha))
        result = subprocess.run([EXE], cwd=ROOT, capture_output=False,
                                text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(f"exocol.exe failed (rc={result.returncode})")
        with nc.Dataset(OUT_NC) as ds:
            ts = float(ds['ts'][:])
        return ts
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def main():
    print("Kopparapu-composition albedo calibration")
    print(f"  Target Ts = {TS_TARGET} K  (tolerance ± {TS_TOL} K)")
    print(f"  Initial α = {ALPHA_INIT:.4f}")
    print()

    alpha = ALPHA_INIT
    for iteration in range(8):
        print(f"  Iteration {iteration+1}: α = {alpha:.6f} ...", flush=True)
        ts = run_one(alpha)
        delta = ts - TS_TARGET
        print(f"    → Ts = {ts:.4f} K  (ΔTs = {delta:+.4f} K)")
        if abs(delta) < TS_TOL:
            print(f"\n  Converged: α = {alpha:.6f}  Ts = {ts:.4f} K")
            break
        alpha -= delta / ALPHA_SENS
        alpha = max(0.10, min(0.50, alpha))   # guard rails
    else:
        print("  WARNING: did not converge within 8 iterations")

    # Write calibrated namelist with documentation comment
    final_nml = NML_TEMPLATE.format(alpha=alpha).replace(
        "  asdir",
        f"  ! Surface albedo calibrated to Ts = {ts:.4f} K (Kopparapu composition,\n"
        f"  ! no O3, MOS surface flux, sbm, prognostic, PVER=70).\n"
        "  asdir"
    )
    cal_nml = os.path.join(REF_DIR, 'exocol_config.nml')
    with open(cal_nml, 'w') as f:
        f.write(final_nml)
    print(f"  Saved: {cal_nml}")

    # Copy the final output NC
    shutil.copy(OUT_NC, os.path.join(REF_DIR, 'exocol_out.nc'))
    print(f"  Saved: {os.path.join(REF_DIR, 'exocol_out.nc')}")

    print("\nDone. Run plot_earth_kopparapu.py to visualise.")

if __name__ == '__main__':
    main()
