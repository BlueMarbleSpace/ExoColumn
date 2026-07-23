#!/usr/bin/env python3
"""
hz_figure6.py  —  ExoColumn analogue of Kopparapu et al. (2013) Figure 6.

Kopparapu+2013 Fig. 6 is a 4-panel habitable-zone summary, one curve per host
star (their F 7200 K, G 5800 K/Sun, K 4800 K, M 3800 K, M 2600 K):

  (a) planetary albedo            vs surface temperature   [inner edge]
  (b) effective stellar flux Seff vs surface temperature   [inner edge]
  (c) planetary albedo            vs CO2 partial pressure   [outer edge]
  (d) Seff                        vs CO2 partial pressure   [outer edge]

The stellar-temperature dependence enters through the planetary albedo: cooler
stars emit more in the near-IR, where H2O (inner edge) and CO2 (outer edge)
absorb strongly, so the atmosphere absorbs more sunlight, the albedo drops, and
the HZ limits shift to lower Seff.  This is exactly what the new runtime
&exocol_nml::solar_file option lets ExoColumn reproduce: one binary, a different
n68 BT-Settl SED per star (ExoRT renormalises each SED to the fixed S0 and
re-optimises the SW bands, so only the spectral SHAPE enters).

Each inner-edge point reuses reference/moist_runaway/hz_inner.run_one (the
validated nonideal-EOS / ice cold-trap / CLIMA-convention IHZ config); each
outer-edge point reuses reference/max_greenhouse/hz_outer.run_one (CO2
condensation + T-dependent cp).  The host star is the only thing that changes
between curves, so this is apples-to-apples across spectral type.

COOL-HALF SET (this first deliverable).  Kopparapu's F 7200 K and exact K
4800 K are deferred until those n68 SEDs are generated; the n68 grid currently
tops out at the Sun with BT-Settl spectra at <= 4500 K.  We map onto the three
coolest Kopparapu types with the spectra we have and add the Sun:

    M 2600 K   bt-settl_2600   (= Kopparapu's 2600 K exactly)
    M 3700 K   bt-settl_3700   (~ Kopparapu's 3800 K)
    K 4500 K   bt-settl_4500   (~ Kopparapu's 4800 K)
    G 5780 K   G2V_SUN (Sun)   (~ Kopparapu's 5800 K; '' = compile-time default)

Requires the binary built at PVER>=200 with the solar_file option
(make clean && make PVER=200) and the Intel OneAPI runtime on the path
(source /opt/intel/oneapi/setvars.sh) so exocol.exe finds libimf.so.

Results are cached to hz_figure6.npz; set HZ_REPLOT=1 to re-plot without
re-running the sweeps.  Coarsen the grids for a quick test with
HZ_TS_STEP / HZ_PCO2_N.
"""

import os
import sys
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
CACHE = os.path.join(HERE, 'hz_figure6.npz')
FIG_PNG = os.path.join(HERE, 'hz_figure6.png')
FIG_PDF = os.path.join(HERE, 'hz_figure6.pdf')


def _load(modname, relpath):
    """Import a sibling reference script by file path (clean, no sys.path edits)."""
    spec = importlib.util.spec_from_file_location(
        modname, os.path.join(ROOT, relpath))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')

# --------------------------------------------------------------------------
# Host-star set (label, n68 solar_file, Teff[K], colour).  '' => Sun (default).
# Colours follow Kopparapu's convention: M red -> K orange -> G green.
# --------------------------------------------------------------------------
STARS = [
    ('M 2000 K', 'btsettl_T2000_g4.5_m0.0_n68.nc',    2000, '#800026'),
    ('M 2200 K', 'btsettl_T2200_g4.5_m0.0_n68.nc',    2200, '#a50f15'),
    ('M 2400 K', 'btsettl_T2400_g4.5_m0.0_n68.nc',    2400, '#cb181d'),
    ('M 2600 K', 'bt-settl_2600_logg4.5_FeH0_n68.nc', 2600, '#d62728'),
    ('M 3000 K', 'bt-settl_3000_logg4.5_FeH0_n68.nc', 3000, '#ff7f0e'),
    ('M 3300 K', 'bt-settl_3300_logg4.5_FeH0_n68.nc', 3300, '#e6b800'),
    ('M 3700 K', 'bt-settl_3700_logg4.5_FeH0_n68.nc', 3700, '#2ca02c'),
    ('K 4000 K', 'bt-settl_4000_logg4.5_FeH0_n68.nc', 4000, '#17becf'),
    ('K 4500 K', 'bt-settl_4500_logg4.5_FeH0_n68.nc', 4500, '#1f77b4'),
    ('K 4800 K', 'bt-settl_4800_logg4.5_FeH0_n68.nc', 4800, '#3b4cc0'),
    ('G 5780 K (Sun)', '',                            5780, '#9467bd'),
]
# These are the n68-core stars (run by this script's sweep).  The F 7200 K star
# is added SEPARATELY at n84 by hz_add_n84_stars.py (the n84 SED resolves its
# shortwave better — n68 underestimates albedo/Seff ~2%).  Hotter A/B stars were
# tried and dropped (unphysical albedo>1; see hz_add_n84_stars.py).
# Kopparapu's 5 types are covered (F7200/G-Sun/K4800/M3700~3800/M2600) plus
# infill. 3000 K ~ Proxima Cen; 3300 K ~ mid-M (AD Leo/GJ 1214); 4000 K ~ K7/M0.
# (Curves are recoloured by Teff in plot(); the hex above is a fallback.)

# Inner-edge surface-temperature grid (Kopparapu Fig 6 a/b span 200-2200 K).
TS_MIN  = float(os.environ.get('HZ_TS_MIN',  '200'))
TS_MAX  = float(os.environ.get('HZ_TS_MAX',  '2200'))
TS_STEP = float(os.environ.get('HZ_TS_STEP', '20'))
TS_VALUES = np.arange(TS_MIN, TS_MAX + TS_STEP, TS_STEP, dtype=float)

# Outer-edge CO2 grid [bar].  Upper end = psat_CO2(273 K) = 34.7 bar (the same
# physical endpoint as hz_outer; beyond it CO2 condenses at the surface itself).
PCO2_MIN = float(os.environ.get('HZ_PCO2_MIN', '1.0'))
PCO2_MAX = float(os.environ.get('HZ_PCO2_MAX', '34.7'))
PCO2_N   = int(os.environ.get('HZ_PCO2_N',     '30'))
PCO2_VALUES = np.geomspace(PCO2_MIN, PCO2_MAX, PCO2_N)


def _med(y, w=5):
    """Light running-median smoother for the vertical-resolution sawtooth
    (same closure hz_inner uses).  NaN-safe; w must be odd."""
    y = np.asarray(y, float)
    n = len(y)
    if n < w:
        return y
    half = w // 2
    out = y.copy()
    for i in range(half, n - half):
        seg = y[i - half:i + half + 1]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            out[i] = np.median(seg)
    return out


def sweep():
    """Run inner (Ts) and outer (pCO2) sweeps for every star; return a dict."""
    data = {'ts': TS_VALUES, 'pco2': PCO2_VALUES, 'stars': []}
    for label, sf, teff, color in STARS:
        star_id = sf if sf else 'G2V_SUN_n68.nc (default)'
        print(f"\n=== {label}  [{star_id}] ===", flush=True)

        print(f"  inner edge: {len(TS_VALUES)} Ts points 200-2200 K", flush=True)
        alb_i = np.full(len(TS_VALUES), np.nan)
        seff_i = np.full(len(TS_VALUES), np.nan)
        for j, ts in enumerate(TS_VALUES):
            r = inner.run_one(ts, solar_file=sf)
            if r is not None:
                alb_i[j], seff_i[j] = r['alpha'], r['seff']
            if j % 20 == 0:
                print(f"    Ts={ts:6.0f} K  alb={alb_i[j]:.3f}  Seff={seff_i[j]:.3f}",
                      flush=True)

        print(f"  outer edge: {len(PCO2_VALUES)} pCO2 points 1-34.7 bar", flush=True)
        alb_o = np.full(len(PCO2_VALUES), np.nan)
        seff_o = np.full(len(PCO2_VALUES), np.nan)
        for j, pco2 in enumerate(PCO2_VALUES):
            r = outer.run_one(pco2, solar_file=sf)
            if r is not None:
                alb_o[j], seff_o[j] = r['alpha'], r['seff']

        mg = np.nanmin(seff_o) if np.any(np.isfinite(seff_o)) else np.nan
        print(f"  max-greenhouse Seff (outer minimum) = {mg:.3f}", flush=True)
        data['stars'].append(dict(label=label, teff=teff, color=color,
                                  alb_i=alb_i, seff_i=seff_i,
                                  alb_o=alb_o, seff_o=seff_o))
    return data


def save_cache(data):
    flat = {'ts': data['ts'], 'pco2': data['pco2']}
    for i, s in enumerate(data['stars']):
        for k in ('alb_i', 'seff_i', 'alb_o', 'seff_o'):
            flat[f's{i}_{k}'] = s[k]
        flat[f's{i}_label'] = s['label']
        flat[f's{i}_teff'] = s['teff']
        flat[f's{i}_color'] = s['color']
    flat['nstars'] = len(data['stars'])
    np.savez(CACHE, **flat)
    print(f"\nCached: {CACHE}")


def load_cache():
    d = np.load(CACHE, allow_pickle=True)
    data = {'ts': d['ts'], 'pco2': d['pco2'], 'stars': []}
    for i in range(int(d['nstars'])):
        data['stars'].append(dict(
            label=str(d[f's{i}_label']), teff=int(d[f's{i}_teff']),
            color=str(d[f's{i}_color']),
            alb_i=d[f's{i}_alb_i'], seff_i=d[f's{i}_seff_i'],
            alb_o=d[f's{i}_alb_o'], seff_o=d[f's{i}_seff_o']))
    return data


def plot(data):
    ts, pco2 = data['ts'], data['pco2']
    fig, ax = plt.subplots(2, 2, figsize=(9.0, 7.2))
    (a, b), (c, dax) = ax
    fig.patch.set_facecolor('white')

    # Colour the curves continuously by Teff (cool = red ... warm = blue,
    # Kopparapu style) and key them with a single colorbar (added below) rather
    # than per-curve labels.  A linear Normalize on Teff makes the colorbar a
    # faithful Teff axis; the turbo segment is truncated to avoid the near-black
    # extremes.
    from matplotlib.colors import Normalize, LinearSegmentedColormap
    import matplotlib.cm as mcm
    teffs = np.array([s['teff'] for s in data['stars']], float)
    norm = Normalize(vmin=teffs.min(), vmax=teffs.max())
    cmap = LinearSegmentedColormap.from_list(
        'hz_teff', plt.get_cmap('turbo')(np.linspace(0.92, 0.10, 256)))
    for s in data['stars']:
        s['color'] = cmap(norm(s['teff']))

    for s in data['stars']:
        col = s['color']
        a.plot(ts, _med(s['alb_i']),  color=col, lw=1.6)
        b.plot(ts, _med(s['seff_i']), color=col, lw=1.6)
        c.plot(pco2, s['alb_o'],  color=col, lw=1.6)
        dax.plot(pco2, s['seff_o'], color=col, lw=1.6)

    # Inner-edge panels (a, b): vs surface temperature
    for axx in (a, b):
        axx.set_xlim(200, 2200)
        axx.set_xlabel('Surface temperature [K]')
    a.set_ylabel('Planetary albedo')
    a.set_ylim(0, 0.40)
    a.set_title('Inner edge', loc='left', fontsize=10)
    b.set_ylabel(r'Effective stellar flux  $S_{\rm eff}$')
    b.set_ylim(0.5, 1.8)   # Kopparapu Fig 6b framing; the nonideal-steam
    # supercritical branch (Ts >~ 1600 K) climbs off-scale to Seff ~ 3.5.
    b.set_title('Inner edge', loc='left', fontsize=10)

    # Outer-edge panels (c, d): vs CO2 partial pressure
    for axx in (c, dax):
        axx.set_xscale('log')
        axx.set_xlim(1, 35)
        axx.set_xlabel('CO$_2$ partial pressure [bar]')
    c.set_ylabel('Planetary albedo')
    c.set_title('Outer edge', loc='left', fontsize=10)
    dax.set_ylabel(r'Effective stellar flux  $S_{\rm eff}$')
    dax.set_title('Outer edge', loc='left', fontsize=10)

    # Single Teff colorbar on the right (replaces the per-curve labels).
    sm = mcm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig.subplots_adjust(left=0.07, right=0.87, top=0.96, bottom=0.07,
                        wspace=0.26, hspace=0.30)
    cax = fig.add_axes((0.895, 0.10, 0.018, 0.80))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label(r'Stellar effective temperature  $T_{\rm eff}$  [K]')
    cb.set_ticks([2000, 3000, 4000, 5000, 6000, 7000])
    fig.savefig(FIG_PNG, dpi=300)
    fig.savefig(FIG_PDF)
    print(f"Wrote: {FIG_PNG}\n       {FIG_PDF}")


def main():
    if os.environ.get('HZ_REPLOT') == '1':
        if not os.path.exists(CACHE):
            raise FileNotFoundError(f"No cache to re-plot: {CACHE}")
        print(f"HZ_REPLOT: re-plotting from {CACHE}")
        data = load_cache()
    else:
        data = sweep()
        save_cache(data)
    plot(data)


if __name__ == '__main__':
    main()
