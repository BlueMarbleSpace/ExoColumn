#!/usr/bin/env python3
"""Run a focused Ts subset through the CURRENT exocol build (whatever PVER it was
compiled with) in flux_only mode and write Ts,OLR,ASR,albedo to a CSV.
Usage: python tools/diag_pver_sweep.py <tmin> <tmax> <tstep> <out.csv>
"""
import os, sys, subprocess
import numpy as np
import netCDF4 as nc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML = os.path.join(ROOT, 'exocol_config.nml')
ALBEDO = 0.24229

TEMPLATE = """\
&exocol_nml
  flux_only=.true.
  variable_ps=.true.
  ihz_profile=.true.
  o3_profile='none'
  msdist=1.0
/
&exocol_init
  input_file=''
  ts={ts:.2f}
  t_strato=200.0
  p_top=1.0
  rh_init=1.0
  coszrs=0.5
  asdir={a:.5f}
  asdif={a:.5f}
  aldir={a:.5f}
  aldif={a:.5f}
/
&exocol_composition
  n2_vmr=0.78
  o2_vmr=0.210
  ar_vmr=0.01
  co2_vmr=3.3e-4
  ch4_vmr=0.0
  o3_vmr=0.0
/
"""

def run_one(ts):
    with open(NML, 'w') as f:
        f.write(TEMPLATE.format(ts=ts, a=ALBEDO))
    r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        return None
    with nc.Dataset(OUT_NC) as ds:
        olr = float(ds['LWUP'][0])
        asr = float(ds['SWDN'][0] - ds['SWUP'][0])
        alb = float(ds['SWUP'][0] / ds['SWDN'][0])
    return olr, asr, alb

def main():
    tmin, tmax, tstep, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    rows = []
    for ts in np.arange(tmin, tmax + 0.1, tstep):
        res = run_one(ts)
        if res is None:
            print(f"  FAIL Ts={ts:.0f}"); continue
        rows.append((ts, *res))
    np.savetxt(out, np.array(rows), header='Ts OLR ASR albedo', fmt='%.4f')
    print(f"wrote {out} ({len(rows)} rows)")

if __name__ == '__main__':
    main()
