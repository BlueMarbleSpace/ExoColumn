#!/usr/bin/env python3
"""
check_zenith_quad.py — validate the sw_zenith_quad / sw_nquad implementation.

Three checks, using flux_only IHZ runs at Ts = 400 K (albedo min) and 1500 K
(plateau), composition/albedo identical to hz_inner.py:

  (1) DEFAULT OFF, coszrs=0.5  must reproduce the pre-change single-angle albedo
      (diag_zenith_albedo "alpha(0.5)": 0.195 @400K, 0.219 @1500K).
  (2) QUAD ON, sw_nquad=6      must reproduce the hemispheric Bond albedo
      (diag_zenith_albedo "A_Bond": ~0.185 @400K, ~0.212 @1500K).
  (3) INSOLATION PRESERVED: SWDN(TOA) for quad-on equals the coszrs=0.5
      single-angle value (sum w_i mu_i = 1/2).

Also prints the 4-point quad to confirm convergence vs 6-point.
The user's exocol_config.nml is backed up and restored.
"""
import os, subprocess
import numpy as np
import netCDF4 as nc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE  = os.path.join(ROOT, 'run', 'exocol.exe')
OUT  = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML  = os.path.join(ROOT, 'exocol_config.nml')

REF_SINGLE = {400.0: 0.195, 1500.0: 0.219}   # diag_zenith_albedo alpha(0.5)
REF_BOND   = {400.0: 0.185, 1500.0: 0.212}   # diag_zenith_albedo A_Bond

TEMPLATE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  ihz_profile    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  h2o_eos        = 'ideal'
  sw_zenith_quad = {quad}
  sw_nquad       = {nq}
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = 200.0
  p_top      = 1.0
  rh_init    = 1.0
  coszrs     = 0.5
  asdir      = 0.3200
  asdif      = 0.3200
  aldir      = 0.3200
  aldif      = 0.3200
/
&exocol_composition
  n2_vmr  = 0.78
  o2_vmr  = 0.210
  ar_vmr  = 0.01
  co2_vmr = 3.3e-4
  ch4_vmr = 0.0
  o3_vmr  = 0.0
/
"""

def run(ts, quad, nq=4):
    with open(NML, 'w') as f:
        f.write(TEMPLATE.format(ts=ts, quad=('.true.' if quad else '.false.'), nq=nq))
    r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        print(f"  FAIL ts={ts} quad={quad} nq={nq} rc={r.returncode}\n{(r.stderr or '')[-300:]}")
        return None
    with nc.Dataset(OUT) as ds:
        swup = float(ds['SWUP'][:][0]); swdn = float(ds['SWDN'][:][0])
    return swup, swdn, swup/swdn

def main():
    orig = open(NML).read() if os.path.exists(NML) else None
    try:
        print("Validation of sw_zenith_quad / sw_nquad\n")
        print(f"  {'Ts':>5}  {'mode':>16}  {'SWDN_toa':>9}  {'albedo':>7}  {'reference':>9}  {'ok?':>4}")
        ok_all = True
        for ts in (400.0, 1500.0):
            swdn_ref = None
            for mode, quad, nq, ref in [
                ('single 0.5', False, 4, REF_SINGLE[ts]),
                ('quad nq=4',  True,  4, REF_BOND[ts]),
                ('quad nq=6',  True,  6, REF_BOND[ts]),
            ]:
                res = run(ts, quad, nq)
                if res is None:
                    ok_all = False; continue
                swup, swdn, alb = res
                if swdn_ref is None:
                    swdn_ref = swdn
                ok = abs(alb - ref) < 0.004
                ins_ok = abs(swdn - swdn_ref) < 0.05
                ok_all = ok_all and ok and ins_ok
                print(f"  {ts:>5.0f}  {mode:>16}  {swdn:>9.3f}  {alb:>7.4f}  {ref:>9.3f}"
                      f"  {'OK' if ok else 'XX':>4}   insol {'OK' if ins_ok else 'DIFF'}")
            print()
        print("ALL CHECKS PASS" if ok_all else "*** SOME CHECKS FAILED ***")
    finally:
        if orig is not None:
            open(NML, 'w').write(orig)
        elif os.path.exists(NML):
            os.remove(NML)
    print("Original exocol_config.nml restored.")

if __name__ == '__main__':
    main()
