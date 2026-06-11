#!/usr/bin/env python3
"""Localize the IHZ panel-(d) residual: overlay ExoColumn vs CLIMA T(P), H2O VMR(P),
and altitude z(P) at a few surface temperatures.

Panel (d) of the moist_runaway figure plots H2O VMR vs altitude; a difference there
can come from (i) the moist pseudoadiabat / cold-point placement (shows up in T(P)
and VMR(P)) or (ii) the hydrostatic altitude mapping (shows up in z(P) at matched
pressure).  This script separates the two by comparing all three against Kopparapu's
CLIMA profiles (clima_last.tab), using the SAME config as the reference figure
(pure N2, liquid cold trap, BPS continuum, nonideal EOS, 6-pt zenith, p_top=0.002).

Requires an ExoColumn binary at PVER>=200.  Run from anywhere (paths via __file__).
"""
import os, subprocess, numpy as np, netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXE  = os.path.join(ROOT, 'run', 'exocol.exe')
NML  = os.path.join(ROOT, 'exocol_config.nml')
OUT  = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
CLIMA = os.path.join(ROOT, 'reference', 'moist_runaway', 'clima_last.tab')
PNG  = os.path.join(HERE, 'diag_ptz_clima.png')
MW_H2O = 18.015
TS_LIST = (340.0, 360.0, 380.0)

NML_T = """&exocol_nml
  flux_only=.true.
  variable_ps=.true.
  ihz_profile=.true.
  o3_profile='none'
  msdist=1.0
  h2o_eos='nonideal'
  h2o_continuum='bps'
  cold_trap_phase='liquid'
  sw_zenith_quad=.true.
  sw_nquad=6
/
&exocol_init
  input_file=''
  ts={ts:.2f}
  t_strato=200.0
  p_top=0.002
  rh_init=1.0
  coszrs=0.5
  asdir=0.32
  asdif=0.32
  aldir=0.32
  aldif=0.32
/
&exocol_composition
  n2_vmr=0.99967
  o2_vmr=0.0
  ar_vmr=0.0
  co2_vmr=3.3e-4
  ch4_vmr=0.0
  o3_vmr=0.0
/
"""


def run_exo(ts):
    open(NML, 'w').write(NML_T.format(ts=ts))
    r = subprocess.run([EXE], cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stderr[-400:]
    with nc.Dataset(OUT) as ds:
        tmid = np.array(ds['tmid'][:]); pmid = np.array(ds['pmid'][:])
        h2ommr = np.array(ds['h2ommr'][:]); mw = float(ds['mw'][:])
        zint = np.array(ds['zint'][:])
    zmid = 0.5 * (zint[:-1] + zint[1:]) / 1000.0
    vmr = h2ommr * mw / MW_H2O
    return dict(T=tmid, P=pmid, z=zmid, vmr=vmr)


def clima_blocks():
    lines = open(CLIMA).read().splitlines()
    hdr = [i for i, l in enumerate(lines) if 'ALT' in l] + [len(lines)]
    cl = {}
    for a, b in zip(hdr[:-1], hdr[1:]):
        rows = []
        for l in lines[a + 1:b]:
            p = l.split()
            if len(p) >= 4:
                try: rows.append([float(x) for x in p[:4]])
                except ValueError: pass
        if len(rows) < 3: continue
        d = np.array(rows)            # alt[km] P[bar] T[K] FH2O
        cl[round(float(d[-1, 2]))] = d
    return cl


def main():
    orig = open(NML).read() if os.path.exists(NML) else None
    cl = clima_blocks()
    res = {}
    try:
        for ts in TS_LIST:
            res[ts] = run_exo(ts)
    finally:
        if orig is not None:
            open(NML, 'w').write(orig)

    fig, ax = plt.subplots(len(TS_LIST), 3, figsize=(11.0, 9.0), dpi=200)
    fig.patch.set_facecolor('white')
    for row, ts in enumerate(TS_LIST):
        e = res[ts]; c = cl[round(ts)]
        cz, cP, cT, cf = c[:, 0], c[:, 1] * 1e5, c[:, 2], c[:, 3]
        # T(P)
        a = ax[row, 0]
        a.semilogy(e['T'], e['P'], '-', color='C3', lw=1.3, label='ExoColumn')
        a.semilogy(cT, cP, '--', color='k', lw=1.0, label='CLIMA')
        a.invert_yaxis(); a.set_ylabel(f'Ts={ts:.0f} K\nP [Pa]')
        if row == 0: a.set_title('T(P)')
        if row == len(TS_LIST) - 1: a.set_xlabel('T [K]')
        if row == 0: a.legend(fontsize=7, frameon=False)
        # VMR(P)
        a = ax[row, 1]
        a.loglog(e['vmr'], e['P'], '-', color='C0', lw=1.3)
        a.loglog(cf, cP, '--', color='k', lw=1.0)
        a.invert_yaxis()
        if row == 0: a.set_title('H2O VMR(P)')
        if row == len(TS_LIST) - 1: a.set_xlabel('H2O VMR')
        # z(P)
        a = ax[row, 2]
        a.semilogy(e['z'], e['P'], '-', color='C2', lw=1.3)
        a.semilogy(cz, cP, '--', color='k', lw=1.0)
        a.invert_yaxis()
        if row == 0: a.set_title('altitude(P)')
        if row == len(TS_LIST) - 1: a.set_xlabel('z [km]')
        for a in ax[row]:
            a.set_facecolor('white')

        # quantitative cold-point comparison
        # ExoColumn cold point: shallowest (lowest-P) level still on the adiabat
        # before the isothermal cap == the tropopause; approximate via min VMR's P
        ecp = e['P'][np.argmin(e['vmr'])]
        # CLIMA cold point: deepest (max-P) 200K level
        m200 = np.abs(cT - cT.min()) < 0.5
        ccp = cP[m200].max()
        print(f"Ts={ts:.0f}: strat VMR exo={e['vmr'].min():.3e} CLIMA={cf.min():.3e} "
              f"(ratio {e['vmr'].min()/cf.min():.2f}); cold-point P  exo~{ecp:.2f} Pa  CLIMA~{ccp:.2f} Pa")

    fig.tight_layout()
    fig.savefig(PNG, facecolor='white')
    print(f"wrote {PNG}")


if __name__ == '__main__':
    main()
