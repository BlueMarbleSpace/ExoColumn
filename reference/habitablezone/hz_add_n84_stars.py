#!/usr/bin/env python3
"""
hz_add_n84_stars.py — add the HOT (n84-core) stars to the HZ figure caches
without re-running the cool n68 stars.

The cool stars (M2600..G-Sun) are n68-core results already cached in
hz_figure6.npz / hz_figure7.npz.  The hot stars need the n84 core: the 10000 K
(A0V) and 20000 K (B2V) BT-Settl SEDs exist only at n84, and the 7200 K F star
is more correct at n84 (resolves its strong shortwave — see
tools/compare_7200_n68_n84.py: n68 underestimates albedo/Seff by ~2%).  So we
replace the n68 F-star with the n84 F-star and append A0V + B2V.

Run with the n84 core built and the Intel runtime on PATH:
    cd build && make clean && make SPEC_DIR=/models/ExoRT/source/src.n84equiv PVER=200
    cd .. && source /opt/intel/oneapi/setvars.sh && \
      python3 reference/habitablezone/hz_add_n84_stars.py
then restore the default n68 core: cd build && make clean && make PVER=200.

It runs the inner (Ts) + outer (pCO2) sweeps on the SAME grids as hz_figure6,
derives the Fig-7 limits (moist-GH at the star-independent Ts=344 K, runaway =
inner-branch peak, max-GH = outer minimum), and MERGES into both caches.
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

# Hot stars (n84 core + n84 SED): (label, solar_file, teff, color)
# Only the F 7200 K star is retained.  The A 10000 K and B 20000 K BT-Settl SEDs
# were tested but DROPPED: their UV/blue-peaked spectra drive ExoRT's shortwave
# two-stream solver into the near-conservative-Rayleigh regime and it returns
# planetary albedo > 1 (unphysical), worse the hotter the star (A: above ~620 K
# surface; B: above ~400 K), corrupting the HZ boundaries (Seff spikes/NaNs).
# The 1-D RT is not valid for hot blue hosts; re-add them only if that is fixed.
HOT = [
    ('F 7200 K',  'bt-settl_7200_logg4.5_FeH0_n84.nc',  7200,  '#7a0177'),
]
REPLACE = 'F 7200 K'   # drop the existing (n68) F-star before appending the n84 one


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run_hot():
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')
    out = []
    for label, sf, teff, color in HOT:
        print(f"\n=== {label}  [{sf}] (n84 core) ===", flush=True)
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


def merge_fig6(hot):
    d = dict(np.load(F6, allow_pickle=True))
    n = int(d['nstars'])
    keep = [i for i in range(n) if str(d[f's{i}_label']) != REPLACE]
    new = {'ts': d['ts'], 'pco2': d['pco2']}
    k = 0
    for i in keep:
        for key in ('alb_i', 'seff_i', 'alb_o', 'seff_o', 'label', 'teff', 'color'):
            new[f's{k}_{key}'] = d[f's{i}_{key}']
        k += 1
    for h in hot:
        for key in ('alb_i', 'seff_i', 'alb_o', 'seff_o', 'label', 'teff', 'color'):
            new[f's{k}_{key}'] = h[key]
        k += 1
    new['nstars'] = k
    np.savez(F6, **new)
    print(f"\nmerged hz_figure6.npz: {k} stars")


def merge_fig7(hot):
    d = dict(np.load(F7, allow_pickle=True))
    lab = [str(x) for x in d['label']]
    keep = [i for i in range(len(lab)) if lab[i] != REPLACE]
    cols = dict(teff=list(np.asarray(d['teff'])[keep]),
                label=list(np.asarray(d['label'])[keep]),
                color=list(np.asarray(d['color'])[keep]),
                seff_inner=list(np.asarray(d['seff_inner'])[keep]),
                ts_inner=list(np.asarray(d['ts_inner'])[keep]),
                seff_outer=list(np.asarray(d['seff_outer'])[keep]),
                seff_runaway=list(np.asarray(d['seff_runaway'])[keep]))
    for h in hot:
        cols['teff'].append(h['teff']); cols['label'].append(h['label'])
        cols['color'].append(h['color']); cols['seff_inner'].append(h['moist'])
        cols['ts_inner'].append(MOIST_TS); cols['seff_outer'].append(h['mx'])
        cols['seff_runaway'].append(h['runaway'])
    np.savez(F7, teff=np.array(cols['teff'], float), label=np.array(cols['label']),
             color=np.array(cols['color']), seff_inner=np.array(cols['seff_inner'], float),
             ts_inner=np.array(cols['ts_inner'], float),
             seff_outer=np.array(cols['seff_outer'], float),
             seff_runaway=np.array(cols['seff_runaway'], float))
    print(f"merged hz_figure7.npz: {len(cols['teff'])} stars")


if __name__ == '__main__':
    hot = run_hot()
    merge_fig6(hot)
    merge_fig7(hot)
