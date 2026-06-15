#!/usr/bin/env python3
"""
compare_7200_n68_n84.py — does stellar-spectrum spectral resolution (n68 vs n84)
change the HZ-relevant radiation for the hot 7200 K F star?

A hot F star peaks in the blue/near-UV and carries strong shortwave features
that the 68-band binning may smear relative to the 84-band binning.  Because the
solar file's band count must match the radiation core, this is NOT a solar-file
swap — it needs the matching ExoColumn build:

    n68 result : n68equiv binary + bt-settl_7200_logg4.5_FeH0_n68.nc
    n84 result : n84equiv binary + bt-settl_7200_logg4.5_FeH0_n84.nc

Build the n84 core (ExoColumn src is band-count-agnostic; just point at the tree):
    cd build && make clean && \
      make SPEC_DIR=/models/ExoRT/source/src.n84equiv PVER=200
and rebuild the default n68 core afterwards with `make clean && make PVER=200`.

Planetary albedo and S_eff are band-integrated scalars, so they are directly
comparable across the two cores.  We run the SAME 7200 K column (flux_only,
nonideal EOS — the hz_inner config) at a set of surface temperatures through
both cores via hz_inner.run_one, so only the radiation treatment differs.

Usage (run once per core, with that core's binary on run/exocol.exe):
    CMP_SOLAR=bt-settl_7200_logg4.5_FeH0_n68.nc CMP_LABEL=n68 python3 tools/compare_7200_n68_n84.py
    CMP_SOLAR=bt-settl_7200_logg4.5_FeH0_n84.nc CMP_LABEL=n84 python3 tools/compare_7200_n68_n84.py
then (no env vars) print the comparison table:
    python3 tools/compare_7200_n68_n84.py
"""

import os
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUTDIR = os.path.join(ROOT, 'reference', 'habitablezone')
TS = [250., 288., 320., 350., 400., 500.]   # K, across the inner branch


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, rel))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(solar, label):
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    rows = []
    print(f"=== 7200 K, {label} core+SED ({solar}) ===", flush=True)
    for ts in TS:
        r = inner.run_one(ts, solar_file=solar)
        if r is None:
            print(f"  Ts={ts:.0f} FAILED", flush=True)
            rows.append([ts, np.nan, np.nan, np.nan, np.nan])
            continue
        rows.append([ts, r['alpha'], r['seff'], r['olr'], r['asr']])
        print(f"  Ts={ts:6.0f}  albedo={r['alpha']:.4f}  Seff={r['seff']:.4f}"
              f"  OLR={r['olr']:7.2f}  ASR={r['asr']:7.2f}", flush=True)
    np.savez(os.path.join(OUTDIR, f'cmp7200_{label}.npz'), rows=np.array(rows))
    print(f"  wrote cmp7200_{label}.npz")


def compare():
    a = np.load(os.path.join(OUTDIR, 'cmp7200_n68.npz'))['rows']
    b = np.load(os.path.join(OUTDIR, 'cmp7200_n84.npz'))['rows']
    print("\n  7200 K F star — n68 vs n84 (planetary albedo and S_eff):")
    print(f"  {'Ts[K]':>5} | {'alb_n68':>8} {'alb_n84':>8} {'Δalb':>8} | "
          f"{'Seff_n68':>8} {'Seff_n84':>8} {'ΔSeff':>8}")
    for i in range(len(a)):
        ts = a[i, 0]
        an, ax = a[i, 1], b[i, 1]
        sn, sx = a[i, 2], b[i, 2]
        print(f"  {ts:5.0f} | {an:8.4f} {ax:8.4f} {ax-an:+8.4f} | "
              f"{sn:8.4f} {sx:8.4f} {sx-sn:+8.4f}")


if __name__ == '__main__':
    _sol, _lab = os.environ.get('CMP_SOLAR'), os.environ.get('CMP_LABEL')
    if _sol and _lab:
        run(_sol, _lab)
    else:
        compare()
