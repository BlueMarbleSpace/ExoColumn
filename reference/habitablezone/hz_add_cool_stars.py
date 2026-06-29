#!/usr/bin/env python3
"""
hz_add_cool_stars.py — add the new ULTRA-COOL (n68-core) BT-Settl stars to the
HZ figure caches without re-running the already-cached M2600..G-Sun + n84 F set.

The four BT-Settl SEDs btsettl_T{2000,2200,2400,2600}_g4.5_m0.0_n68.nc extend the
host-star ladder below the previous 2600 K floor.  Per the agreed convention the
existing 2600 K point (bt-settl_2600_logg4.5_FeH0_n68.nc) is LEFT in place, so
this helper only adds the three genuinely-new cool stars (2000/2200/2400 K).

Mirrors hz_add_n84_stars.py: it runs the SAME inner (Ts) + outer (pCO2) grids as
hz_figure6, derives the Fig-7 limits (moist-GH at the star-independent Ts=344 K,
runaway = inner-branch peak, max-GH = outer minimum), and MERGES into both
caches — INSERTING the cool stars at the FRONT so the fig6 cache stays
index-aligned with hz_figure6/7's STARS list (cool ... Sun), with the n84 F 7200 K
star preserved at the end.

Run with the default n68 core built (PVER>=200) and the Intel runtime on PATH:
    cd build && make clean && make PVER=200          # if not already built
    cd .. && source /opt/intel/oneapi/setvars.sh && \
      python3 reference/habitablezone/hz_add_cool_stars.py
then re-plot:
    HZ_REPLOT=1 python3 reference/habitablezone/hz_figure6.py
    HZ_REPLOT=1 python3 reference/habitablezone/hz_figure7.py
"""

import os
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
F6 = os.path.join(HERE, 'hz_figure6.npz')
F7 = os.path.join(HERE, 'hz_figure7.npz')

TS = np.arange(200., 2200. + 20., 20.)        # matches hz_figure6 inner grid
PCO2 = np.geomspace(1.0, 34.7, 30)            # matches hz_figure6 outer grid
MOIST_TS = 344.0                              # star-independent moist-GH Ts
RUN_LO, RUN_HI = 280.0, 700.0                 # runaway = peak of Seff(Ts) here

# New ultra-cool stars (n68 core + n68 SED): (label, solar_file, teff, color).
# 2600 K is intentionally NOT here (the existing bt-settl_2600 point is kept).
COOL = [
    ('M 2000 K', 'btsettl_T2000_g4.5_m0.0_n68.nc', 2000, '#800026'),
    ('M 2200 K', 'btsettl_T2200_g4.5_m0.0_n68.nc', 2200, '#a50f15'),
    ('M 2400 K', 'btsettl_T2400_g4.5_m0.0_n68.nc', 2400, '#cb181d'),
]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run_cool():
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')
    out = []
    for label, sf, teff, color in COOL:
        print(f"\n=== {label}  [{sf}] (n68 core) ===", flush=True)
        ai = np.full(len(TS), np.nan); si = np.full(len(TS), np.nan)
        for j, ts in enumerate(TS):
            r = inner.run_one(ts, solar_file=sf)
            if r is not None:
                ai[j], si[j] = r['alpha'], r['seff']
            if j % 25 == 0:
                print(f"    Ts={ts:6.0f}  alb={ai[j]:.3f}  Seff={si[j]:.3f}", flush=True)
        ao = np.full(len(PCO2), np.nan); so = np.full(len(PCO2), np.nan)
        for j, p in enumerate(PCO2):
            r = outer.run_one(p, solar_file=sf)
            if r is not None:
                ao[j], so[j] = r['alpha'], r['seff']
        br = (TS >= RUN_LO) & (TS <= RUN_HI)
        moist = float(np.interp(MOIST_TS, TS, si))
        runaway = float(np.nanmax(si[br]))
        mx = float(np.nanmin(so)) if np.any(np.isfinite(so)) else np.nan
        print(f"  -> moist={moist:.3f}  runaway={runaway:.3f}  max-GH={mx:.3f}", flush=True)
        out.append(dict(label=label, teff=teff, color=color, alb_i=ai, seff_i=si,
                        alb_o=ao, seff_o=so, moist=moist, runaway=runaway, mx=mx))
    return out


def merge_fig6(cool):
    d = dict(np.load(F6, allow_pickle=True))
    n = int(d['nstars'])
    existing = [str(d[f's{i}_label']) for i in range(n)]
    # Skip any cool star already present (idempotent re-runs).
    add = [c for c in cool if c['label'] not in existing]
    new = {'ts': d['ts'], 'pco2': d['pco2']}
    k = 0
    for c in add:                                  # cool stars FIRST
        for key in ('alb_i', 'seff_i', 'alb_o', 'seff_o'):
            new[f's{k}_{key}'] = c[key]
        new[f's{k}_label'] = c['label']; new[f's{k}_teff'] = c['teff']
        new[f's{k}_color'] = c['color']
        k += 1
    for i in range(n):                             # then the existing stars
        for key in ('alb_i', 'seff_i', 'alb_o', 'seff_o', 'label', 'teff', 'color'):
            new[f's{k}_{key}'] = d[f's{i}_{key}']
        k += 1
    new['nstars'] = k
    np.savez(F6, **new)
    print(f"\nmerged hz_figure6.npz: {k} stars (added {len(add)} cool)")


def merge_fig7(cool):
    d = dict(np.load(F7, allow_pickle=True))
    lab = [str(x) for x in d['label']]
    add = [c for c in cool if c['label'] not in lab]
    cols = dict(teff=list(np.asarray(d['teff'])),
                label=list(np.asarray(d['label'])),
                color=list(np.asarray(d['color'])),
                seff_inner=list(np.asarray(d['seff_inner'])),
                ts_inner=list(np.asarray(d['ts_inner'])),
                seff_outer=list(np.asarray(d['seff_outer'])),
                seff_runaway=list(np.asarray(d['seff_runaway'])))
    for c in add:
        cols['teff'].append(c['teff']); cols['label'].append(c['label'])
        cols['color'].append(c['color']); cols['seff_inner'].append(c['moist'])
        cols['ts_inner'].append(MOIST_TS); cols['seff_outer'].append(c['mx'])
        cols['seff_runaway'].append(c['runaway'])
    np.savez(F7, teff=np.array(cols['teff'], float), label=np.array(cols['label']),
             color=np.array(cols['color']), seff_inner=np.array(cols['seff_inner'], float),
             ts_inner=np.array(cols['ts_inner'], float),
             seff_outer=np.array(cols['seff_outer'], float),
             seff_runaway=np.array(cols['seff_runaway'], float))
    print(f"merged hz_figure7.npz: {len(cols['teff'])} stars (added {len(add)} cool)")


if __name__ == '__main__':
    cool = run_cool()
    merge_fig6(cool)
    merge_fig7(cool)
