#!/usr/bin/env python3
"""
hz_figure7.py  —  ExoColumn analogue of Kopparapu et al. (2013) Figure 7(a).

Kopparapu+2013 Fig. 7 plots the habitable-zone boundaries as a function of
stellar effective temperature: panel (a) the effective stellar flux S_eff, and
panel (b) the corresponding orbital distance.  Per the user we omit the Selsis
et al. (2007) comparison curves.

This script produces panel (a) — S_eff(T_eff) for the two boundaries used in
Fig. 7:

    inner edge  = MOIST greenhouse   (stratospheric H2O VMR reaches 3e-3)
    outer edge  = MAXIMUM greenhouse (S_eff minimum of the dense-CO2 sweep)

both fully from ExoColumn.  Panel (b) (distance vs stellar mass) needs an
external T_eff <-> mass <-> luminosity main-sequence relation to turn S_eff
into AU, so it is handled separately once that relation is chosen.

METHOD NOTE (why this is cheap and self-consistent).  In flux_only cold-start
mode the column T/H2O profile is built from Ts + composition ALONE — the host
star never enters it.  So OLR and the stratospheric-H2O moist-GH trigger are
star-INDEPENDENT; the moist-greenhouse limit falls at the SAME Ts for every
star, and only S_eff = OLR/ASR differs (via the SED-dependent albedo).  We
therefore:
  * INNER: run a focused Ts sweep (reusing the validated hz_inner.run_one,
    which returns the cold-trap H2O VMR), find the Ts where it crosses 3e-3,
    and read S_eff there for each star.
  * OUTER: read the max-greenhouse S_eff (= min of seff_o) straight out of the
    Figure-6 cache (reference/habitablezone/hz_figure6.npz) — no re-run.

COOL-HALF set (M2600/M3700/K4500/G-Sun); F 7200 K and K 4800 K append once
those n68 SEDs exist.  Requires the PVER>=200 binary with &exocol_nml::solar_file
and the Intel runtime on PATH (source /opt/intel/oneapi/setvars.sh).

Cache: hz_figure7.npz (HZ_REPLOT=1 re-plots without re-running the inner sweep).
"""

import os
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
FIG6_CACHE = os.path.join(HERE, 'hz_figure6.npz')
CACHE = os.path.join(HERE, 'hz_figure7.npz')
FIG_PNG = os.path.join(HERE, 'hz_figure7.png')
FIG_PDF = os.path.join(HERE, 'hz_figure7.pdf')

MOIST_GH_VMR = 3.0e-3   # Kopparapu Sec 3.1 moist-greenhouse stratospheric H2O VMR

# Same star set / order as hz_figure6 (index-aligned with its cache).
STARS = [
    ('M 2600 K', 'bt-settl_2600_logg4.5_FeH0_n68.nc', 2600, '#d62728'),
    ('M 3700 K', 'bt-settl_3700_logg4.5_FeH0_n68.nc', 3700, '#e377c2'),
    ('K 4500 K', 'bt-settl_4500_logg4.5_FeH0_n68.nc', 4500, '#ff7f0e'),
    ('G 5780 K (Sun)', '',                            5780, '#2ca02c'),
]

# --- Baraffe et al. (1998) solar-metallicity 5 Gyr isochrone --------------
# Per Kopparapu (2013, p.12) the HZ distances use the Baraffe+1998
# solar-metallicity evolutionary models at 5 Gyr, fed into Eq. (3):
#     d = sqrt( (L/Lsun) / S_eff )  AU.
# Rows (Teff[K], Mass[Msun], Mbol[mag]) from VizieR J/A+A/337/403 (tab1-3,
# [M/H]=0, Age~5 Gyr).  Baraffe gives two mixing-length branches: Lmix=1.0
# (cool stars, full mass range) and Lmix=1.9 (solar-calibrated, M>=0.7); per
# their prescription we use Lmix=1.9 for Teff>=4400 K and Lmix=1.0 below.
# Luminosity is recovered from Mbol with Baraffe's adopted Mbol_sun = 4.64
# (self-consistent: their 1 Msun/5 Gyr model -> L=1.03 Lsun):
#     log10(L/Lsun) = 0.4 * (Mbol_sun - Mbol).
MBOL_SUN = 4.64
BARAFFE_TEFF_SPLIT = 4400.0   # >= -> Lmix=1.9 branch; below -> Lmix=1.0

# Lmix=1.0 branch: Teff[K], Mass[Msun], Mbol
_BAR_L10 = np.array([
    [2003, 0.075, 14.46], [2313, 0.080, 13.66], [2641, 0.090, 12.81],
    [2812, 0.100, 12.32], [2921, 0.110, 11.97], [3055, 0.130, 11.47],
    [3151, 0.150, 11.09], [3235, 0.175, 10.72], [3292, 0.200, 10.42],
    [3376, 0.250, 9.93],  [3436, 0.300, 9.54],  [3475, 0.350, 9.23],
    [3522, 0.400, 8.93],  [3582, 0.450, 8.62],  [3649, 0.500, 8.30],
    [3804, 0.570, 7.81],  [3893, 0.600, 7.58],  [4239, 0.700, 6.84],
    [4654, 0.800, 6.11],  [5044, 0.900, 5.41],  [5399, 1.000, 4.70],
])
# Lmix=1.9 branch (solar-calibrated, M>=0.7): Teff[K], Mass[Msun], Mbol
_BAR_L19 = np.array([
    [4418, 0.700, 6.74], [4910, 0.800, 6.01],
    [5370, 0.900, 5.31], [5814, 1.000, 4.61],
])


def baraffe_mass_lum(teff):
    """(mass[Msun], L/Lsun) on the Baraffe+1998 5 Gyr [M/H]=0 isochrone."""
    tab = _BAR_L19 if teff >= BARAFFE_TEFF_SPLIT else _BAR_L10
    mass = float(np.interp(teff, tab[:, 0], tab[:, 1]))
    mbol = float(np.interp(teff, tab[:, 0], tab[:, 2]))
    return mass, 10.0 ** (0.4 * (MBOL_SUN - mbol))


# Focused inner-edge Ts grid bracketing the moist-GH crossing (~350 K).
TS_MIN  = float(os.environ.get('HZ7_TS_MIN',  '300'))
TS_MAX  = float(os.environ.get('HZ7_TS_MAX',  '420'))
TS_STEP = float(os.environ.get('HZ7_TS_STEP', '4'))
TS_VALUES = np.arange(TS_MIN, TS_MAX + TS_STEP, TS_STEP, dtype=float)


def _load(modname, relpath):
    spec = importlib.util.spec_from_file_location(
        modname, os.path.join(ROOT, relpath))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _moist_gh_seff(ts, strat_vmr, seff):
    """S_eff at the Ts where the cold-trap H2O VMR first reaches 3e-3
    (log-linear interpolation in Ts)."""
    ts, strat_vmr, seff = map(np.asarray, (ts, strat_vmr, seff))
    ok = np.isfinite(strat_vmr) & np.isfinite(seff) & (strat_vmr > 0)
    ts, strat_vmr, seff = ts[ok], strat_vmr[ok], seff[ok]
    above = np.where(strat_vmr >= MOIST_GH_VMR)[0]
    if above.size == 0 or above[0] == 0:
        return np.nan, np.nan
    i = above[0]
    lo, hi = i - 1, i
    f = (np.log(MOIST_GH_VMR) - np.log(strat_vmr[lo])) / \
        (np.log(strat_vmr[hi]) - np.log(strat_vmr[lo]))
    ts_mg = ts[lo] + f * (ts[hi] - ts[lo])
    seff_mg = seff[lo] + f * (seff[hi] - seff[lo])
    return ts_mg, seff_mg


def sweep():
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    fig6 = np.load(FIG6_CACHE)   # outer max-GH from the Figure-6 sweep

    teffs, colors, labels = [], [], []
    seff_inner, ts_inner, seff_outer = [], [], []
    for i, (label, sf, teff, color) in enumerate(STARS):
        # OUTER: max-greenhouse = minimum of the dense-CO2 S_eff sweep (cached).
        seff_o = fig6[f's{i}_seff_o']
        mg_out = float(np.nanmin(seff_o)) if np.any(np.isfinite(seff_o)) else np.nan

        # INNER: focused Ts sweep for the moist-GH crossing.
        print(f"\n=== {label} (inner moist-GH sweep) ===", flush=True)
        sv = np.full(len(TS_VALUES), np.nan)
        se = np.full(len(TS_VALUES), np.nan)
        for j, ts in enumerate(TS_VALUES):
            r = inner.run_one(ts, solar_file=sf)
            if r is not None:
                sv[j], se[j] = r['strat_vmr'], r['seff']
        ts_mg, mg_in = _moist_gh_seff(TS_VALUES, sv, se)
        print(f"  moist-GH: Ts={ts_mg:.1f} K  S_eff(inner)={mg_in:.3f}   "
              f"max-GH: S_eff(outer)={mg_out:.3f}", flush=True)

        teffs.append(teff); colors.append(color); labels.append(label)
        seff_inner.append(mg_in); ts_inner.append(ts_mg); seff_outer.append(mg_out)

    return dict(teff=np.array(teffs, float), color=np.array(colors),
                label=np.array(labels), seff_inner=np.array(seff_inner),
                ts_inner=np.array(ts_inner), seff_outer=np.array(seff_outer))


def plot(d):
    order = np.argsort(d['teff'])
    teff = d['teff'][order]
    si = d['seff_inner'][order]
    so = d['seff_outer'][order]

    # Baraffe+1998 5 Gyr mass & luminosity -> HZ distances via Kopparapu Eq. (3).
    mass = np.array([baraffe_mass_lum(t)[0] for t in teff])
    lum  = np.array([baraffe_mass_lum(t)[1] for t in teff])
    d_in  = np.sqrt(lum / si)   # inner (moist GH) distance [AU]
    d_out = np.sqrt(lum / so)   # outer (max GH) distance [AU]

    RED, BLUE, GRN = '#d62728', '#1f6feb', '#b8e6b8'
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.0, 5.0))
    fig.patch.set_facecolor('white')

    # ---- Panel (a): HZ fluxes — S_eff vs Teff --------------------------------
    axa.plot(si, teff, 'o-', color=RED, lw=1.8, ms=6)
    axa.plot(so, teff, 's-', color=BLUE, lw=1.8, ms=6)
    axa.fill_betweenx(teff, so, si, color=GRN, alpha=0.45, zorder=0)
    axa.set_xlim(1.25, 0.15)            # reversed: high flux (inner) on the left
    axa.set_ylim(2400, 6100)
    axa.set_xlabel(r'Effective flux incident on the planet  $S/S_0$')
    axa.set_ylabel(r'Stellar effective temperature  $T_{\rm eff}$  [K]')
    axa.set_title('(a)  HZ fluxes', loc='left', fontsize=10)
    axa.grid(alpha=0.25, lw=0.5)
    axa.text(si[-1] + 0.05, teff[-1] + 30, 'Inner edge\n(moist greenhouse)',
             color=RED, fontsize=8.5, ha='left', va='bottom')
    axa.text(so[-1] - 0.02, teff[-1] + 30, 'Outer edge\n(maximum greenhouse)',
             color=BLUE, fontsize=8.5, ha='right', va='bottom')
    axa.text(0.63, 3550, 'Habitable\nzone', color='#2f7d2f', fontsize=9.5,
             ha='center', va='center', weight='bold')

    # ---- Panel (b): HZ distances — distance vs stellar mass ------------------
    axb.plot(d_in, mass, 'o-', color=RED, lw=1.8, ms=6)
    axb.plot(d_out, mass, 's-', color=BLUE, lw=1.8, ms=6)
    axb.fill_betweenx(mass, d_in, d_out, color=GRN, alpha=0.45, zorder=0)
    axb.set_xscale('log'); axb.set_yscale('log')
    axb.set_xlim(0.02, 2.2)
    axb.set_ylim(0.07, 1.25)
    axb.set_yticks([0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
    axb.set_yticklabels(['0.1', '0.2', '0.3', '0.5', '0.7', '1.0'])
    axb.set_xticks([0.03, 0.1, 0.3, 1.0, 2.0])
    axb.set_xticklabels(['0.03', '0.1', '0.3', '1.0', '2.0'])
    axb.set_xlabel('Distance [AU]')
    axb.set_ylabel(r'Stellar mass  [$M_\odot$]')
    axb.set_title('(b)  HZ distances  (Baraffe et al. 1998, 5 Gyr)',
                  loc='left', fontsize=10)
    axb.grid(alpha=0.25, lw=0.5, which='both')
    axb.text(d_in[-1], mass[-1] * 1.04, 'inner', color=RED, fontsize=8.5,
             ha='center', va='bottom')
    axb.text(d_out[-1], mass[-1] * 1.04, 'outer', color=BLUE, fontsize=8.5,
             ha='center', va='bottom')

    fig.suptitle('ExoColumn habitable zone vs host star — '
                 'Kopparapu (2013) Fig. 7 analogue (cool-half set)', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG_PNG, dpi=300)
    fig.savefig(FIG_PDF)
    print(f"\nWrote: {FIG_PNG}\n       {FIG_PDF}")
    print("\nBaraffe 5 Gyr mapping + Eq.(3) distances:")
    for t, m, l, di, do in zip(teff, mass, lum, d_in, d_out):
        print(f"  Teff={t:5.0f} K  M={m:.3f} Msun  L/Lsun={l:.4g}  "
              f"d_inner={di:.3f} AU  d_outer={do:.3f} AU")


def main():
    if os.environ.get('HZ_REPLOT') == '1' and os.path.exists(CACHE):
        print(f"HZ_REPLOT: re-plotting from {CACHE}")
        z = np.load(CACHE)
        d = {k: z[k] for k in z.files}
    else:
        d = sweep()
        np.savez(CACHE, **d)
        print(f"Cached: {CACHE}")
    plot(d)


if __name__ == '__main__':
    main()
