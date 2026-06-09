#!/usr/bin/env python3
"""Quick diagnostic: model-top altitude (zint[0]) of the IHZ nonideal cold-start
column as a function of (Ts, p_top).  Answers 'how low must p_top be so every
panel-(d) profile reaches >=100 km'.  Coldest profiled Ts (280 K) is the binding
case.  Does NOT touch hz_inner figures."""
import os, subprocess, numpy as np, netCDF4 as nc

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE  = os.path.join(ROOT, 'run', 'exocol.exe')
OUT  = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML  = os.path.join(ROOT, 'exocol_config.nml')

TEMPLATE = """\
&exocol_nml
  flux_only = .true.
  variable_ps = .true.
  ihz_profile = .true.
  o3_profile = 'none'
  msdist = 1.0
  h2o_eos = 'nonideal'
  sw_zenith_quad = .true.
  sw_nquad = 6
/
&exocol_init
  input_file = ''
  ts = {ts:.2f}
  t_strato = 200.0
  p_top = {ptop:.6g}
  rh_init = 1.0
  coszrs = 0.5
  asdir = 0.32
  asdif = 0.32
  aldir = 0.32
  aldif = 0.32
/
&exocol_composition
  n2_vmr = 0.78
  o2_vmr = 0.210
  ar_vmr = 0.01
  co2_vmr = 3.3e-4
  ch4_vmr = 0.0
  o3_vmr = 0.0
/
"""

def run(ts, ptop):
    orig = open(NML).read() if os.path.exists(NML) else None
    try:
        open(NML, 'w').write(TEMPLATE.format(ts=ts, ptop=ptop))
        r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(f"  Ts={ts} ptop={ptop} FAILED rc={r.returncode}\n{r.stderr[-300:]}")
            return None
        with nc.Dataset(OUT) as ds:
            zint = np.array(ds['zint'][:])
        return zint[0] / 1000.0   # km
    finally:
        if orig is not None:
            open(NML, 'w').write(orig)
        elif os.path.exists(NML):
            os.remove(NML)

if __name__ == '__main__':
    print("Ts[K]   p_top[Pa]   top[km]")
    for ptop in (0.01, 0.005, 0.002, 0.001):
        for ts in (280.0, 300.0, 380.0):
            z = run(ts, ptop)
            if z is not None:
                print(f"{ts:6.0f} {ptop:10.4g} {z:9.1f}")
