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
]  # n68-core stars (index-aligned with hz_figure6.py); the F 7200 K star is
   # added at n84 by hz_add_n84_stars.py (A/B hotter stars dropped, alpha>1).

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


# The Baraffe+1998 5 Gyr isochrone tops out at 1.0 Msun / ~5814 K.  Hotter stars
# (the F 7200 K ~ 1.6 Msun) have left the main sequence by 5 Gyr, so the 5 Gyr
# isochrone has no point for them; for those we fall back to PRESENT-DAY
# main-sequence values (Pecaut & Mamajek 2013), (Teff[K], mass[Msun], logL).
# This mixes a present-day MS point into an otherwise-5 Gyr panel for the warm
# end — noted in the caption.  F0V = 7220 K.
BARAFFE_TEFF_MAX = 5814.0
# Present-day MS sequence (Pecaut & Mamajek 2013) for the hot stars beyond the
# Baraffe 5 Gyr range: (Teff[K], mass[Msun], logL).  F0V .. B2V.
_PM_HOT = np.array([
    [7220.0,  1.61, 0.86],    # F0V  (7200 K case)
    [9700.0,  2.18, 1.58],    # A0V
    [10400.0, 2.68, 1.80],    # B9.5V  (brackets the 10000 K case)
    [17000.0, 5.40, 2.99],    # B3V
    [18500.0, 6.10, 3.20],    # B2.5V
    [20600.0, 7.30, 3.43],    # B2V   (brackets the 20000 K case)
])


def star_mass_lum(teff):
    """(mass[Msun], L/Lsun) for the host star: Baraffe+1998 5 Gyr isochrone up to
    5814 K, present-day MS (Pecaut & Mamajek 2013) above it."""
    if teff > BARAFFE_TEFF_MAX:
        m = float(np.interp(teff, _PM_HOT[:, 0], _PM_HOT[:, 1]))
        ll = float(np.interp(teff, _PM_HOT[:, 0], _PM_HOT[:, 2]))
        return m, 10.0 ** ll
    tab = _BAR_L19 if teff >= BARAFFE_TEFF_SPLIT else _BAR_L10
    mass = float(np.interp(teff, tab[:, 0], tab[:, 1]))
    mbol = float(np.interp(teff, tab[:, 0], tab[:, 2]))
    return mass, 10.0 ** (0.4 * (MBOL_SUN - mbol))


# --- Kopparapu et al. (2013) Table 3 parametric fit (their Eq. 2) ----------
# S_eff = Seff_sun + a*T* + b*T*^2 + c*T*^3 + d*T*^4,  with T* = Teff - 5780 K.
# Used to overlay their published boundaries as dashed comparison curves.
KOPP_COEFFS = {  # (Seff_sun, a, b, c, d)
    'runaway': (1.0512, 1.3242e-4, 1.5418e-8, -7.9895e-12, -1.8328e-15),
    'moist':   (1.0140, 8.1774e-5, 1.7063e-9, -4.3241e-12, -6.6462e-16),
    'maxgh':   (0.3438, 5.8942e-5, 1.6558e-9, -3.0045e-12, -5.2983e-16),
}


def kopp_seff(teff, boundary):
    """Kopparapu (2013) Eq. (2) S_eff at the given boundary, for Teff array."""
    s0, a, b, c, dd = KOPP_COEFFS[boundary]
    t = np.asarray(teff, float) - 5780.0
    return s0 + a * t + b * t**2 + c * t**3 + dd * t**4


# Runaway-greenhouse S_eff = the MAXIMUM (Simpson-Nakajima) value of the
# inner-branch S_eff(Ts) — the peak flux a saturated atmosphere can radiate
# away before runaway.  Searched over the moist branch (Ts 280-440 K, the peak
# is at ~300-340 K), below the supercritical-steam rise.  The upper bound is
# 440 K (NOT 700 K) to stay consistent with reference/planet_mass/hz_mass.py's
# runaway definition: for the hot F 7200 K star the supercritical branch begins
# by ~600 K, so a 700 K bracket would grab the supercritical climb (Seff~1.35 at
# 680 K) instead of the SN knee (~1.24), making the 1 M_E runaway curve disagree
# between the two figures.  Cooler stars are unaffected (their curve declines
# after the ~320 K peak; supercritical onset is ~1750 K).
RUNAWAY_TS_LO, RUNAWAY_TS_HI = 280.0, 440.0

# Main-sequence spectral-type boundaries, in Teff [K] (panel a) and mass [Msun]
# (panel b), with the type letter for each band (M..A within the plotted range).
# Topmost plotted type is F (the 7200 K / 1.61 Msun host); the A band is omitted.
SPT_TEFF_BOUNDS  = [3900., 5300., 6000.]
SPT_TEFF_CENTERS = [3100., 4600., 5650., 6650.]
SPT_MASS_BOUNDS  = [0.60, 0.88, 1.04]
SPT_MASS_CENTERS = [0.182, 0.727, 0.957, 1.330]
SPT_LABELS = ['M', 'K', 'G', 'F']


def add_spectral_axis(ax, bounds, centers, labels):
    """Right-hand secondary y-axis: main-sequence spectral-type bands.  The type
    letter sits at each band centre; short ticks mark the type boundaries."""
    axr = ax.secondary_yaxis('right')
    axr.set_yticks(centers)
    axr.set_yticklabels(labels)
    axr.set_yticks(bounds, minor=True)
    axr.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())  # no numbers
    axr.tick_params(axis='y', which='major', length=0)   # letters only, no mark
    axr.tick_params(axis='y', which='minor', length=5)    # boundary ticks
    axr.set_ylabel('Spectral type')


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

    fig6_ts = fig6['ts']         # inner-edge Ts grid from the Figure-6 sweep

    teffs, colors, labels = [], [], []
    seff_inner, ts_inner, seff_outer, seff_runaway = [], [], [], []
    for i, (label, sf, teff, color) in enumerate(STARS):
        # OUTER: max-greenhouse = minimum of the dense-CO2 S_eff sweep (cached).
        seff_o = fig6[f's{i}_seff_o']
        mg_out = float(np.nanmin(seff_o)) if np.any(np.isfinite(seff_o)) else np.nan
        # RUNAWAY: peak of the inner-edge S_eff(Ts) curve (Fig-6 inner sweep).
        # Star-independent profile -> same peak Ts, S_eff differs only via albedo.
        seff_i_full = fig6[f's{i}_seff_i']
        br = (fig6_ts >= RUNAWAY_TS_LO) & (fig6_ts <= RUNAWAY_TS_HI)
        run_in = float(np.nanmax(seff_i_full[br]))

        # INNER (moist GH): focused Ts sweep for the 3e-3 crossing.
        print(f"\n=== {label} (inner moist-GH sweep) ===", flush=True)
        sv = np.full(len(TS_VALUES), np.nan)
        se = np.full(len(TS_VALUES), np.nan)
        for j, ts in enumerate(TS_VALUES):
            r = inner.run_one(ts, solar_file=sf)
            if r is not None:
                sv[j], se[j] = r['strat_vmr'], r['seff']
        ts_mg, mg_in = _moist_gh_seff(TS_VALUES, sv, se)
        print(f"  moist-GH: Ts={ts_mg:.1f} K  S_eff={mg_in:.3f}   "
              f"runaway-GH (peak): S_eff={run_in:.3f}   "
              f"max-GH: S_eff={mg_out:.3f}", flush=True)

        teffs.append(teff); colors.append(color); labels.append(label)
        seff_inner.append(mg_in); ts_inner.append(ts_mg)
        seff_outer.append(mg_out); seff_runaway.append(run_in)

    return dict(teff=np.array(teffs, float), color=np.array(colors),
                label=np.array(labels), seff_inner=np.array(seff_inner),
                ts_inner=np.array(ts_inner), seff_outer=np.array(seff_outer),
                seff_runaway=np.array(seff_runaway))


def plot(d):
    from matplotlib.lines import Line2D
    order = np.argsort(d['teff'])
    teff = d['teff'][order]
    s_moist = d['seff_inner'][order]      # moist greenhouse (inner edge)
    s_max   = d['seff_outer'][order]      # maximum greenhouse (outer edge)
    if 'seff_runaway' in d:
        s_run = np.asarray(d['seff_runaway'], float)[order]
    else:                                  # fallback for an older cache
        f6 = np.load(FIG6_CACHE); t6 = f6['ts']
        br = (t6 >= RUNAWAY_TS_LO) & (t6 <= RUNAWAY_TS_HI)
        s_run = np.array([np.nanmax(f6[f's{i}_seff_i'][br])
                          for i in range(len(STARS))])[order]

    # Baraffe+1998 5 Gyr mass & luminosity -> HZ distances via Kopparapu Eq. (3).
    mass = np.array([star_mass_lum(t)[0] for t in teff])
    lum  = np.array([star_mass_lum(t)[1] for t in teff])
    d_run, d_moist, d_max = (np.sqrt(lum / s) for s in (s_run, s_moist, s_max))

    # Smooth Teff grid for the dashed Kopparapu (2013) Eq.(2)+Table-3 overlay in
    # panel (a); panel (b) uses per-star points (avoids the G->F stellar gap).
    # The Kopparapu parametric fit is stated valid only for 2600-7200 K, so the
    # dashed overlay is capped at 2600 K even though our curves run down to 2000 K.
    KOPP_TEFF_FLOOR = 2600.0
    tk = np.linspace(max(KOPP_TEFF_FLOOR, teff.min()), teff.max(), 120)

    RUN, MOIST, MAX, GRN = '#ff7f0e', '#d62728', '#1f6feb', '#b8e6b8'
    BNDS = [('runaway', RUN, '^'), ('moist', MOIST, 'o'), ('maxgh', MAX, 's')]
    fig, (axa, axb) = plt.subplots(2, 1, figsize=(7.5, 8.5))
    fig.patch.set_facecolor('white')

    # ---- Panel (a): HZ fluxes — S_eff vs Teff --------------------------------
    axa.fill_betweenx(teff, s_max, s_moist, color=GRN, alpha=0.45, zorder=0)
    for s, (_, col, _m) in zip((s_run, s_moist, s_max), BNDS):
        axa.plot(s, teff, ls='-', color=col, lw=1.8)
    for bnd, col, _ in BNDS:                      # Kopparapu 2013 (dashed)
        axa.plot(kopp_seff(tk, bnd), tk, '--', color=col, lw=1.3)
    axa.set_xlim(1.40, 0.15)            # reversed: high flux (inner) on the left
    axa.set_ylim(1900, 7500)            # extended to the 2000 K BT-Settl floor
    axa.set_xlabel(r'Effective flux incident on the planet  $S/S_0$')
    axa.set_ylabel(r'Stellar effective temperature  $T_{\rm eff}$  [K]')
    axa.text(0.60, 4700, 'Habitable\nzone', color='#2f7d2f', fontsize=9.5,
             ha='center', va='center', weight='bold')
    # Boundary identity labelled directly on the plot (colour-matched leaders);
    # runaway/maximum above the F-star (7200 K) curve tops, moist below the curves.
    axa.text(1.36, 7300, 'Runaway greenhouse',
             color=RUN, ha='left', va='center', fontsize=8.5)
    axa.annotate('Moist greenhouse', xy=(0.90, 2740), xytext=(1.03, 2470),
                 color=MOIST, ha='center', va='center', fontsize=8.5,
                 arrowprops=dict(arrowstyle='->', color=MOIST, lw=0.7))
    axa.text(0.46, 7300, 'Maximum greenhouse',
             color=MAX, ha='center', va='center', fontsize=8.5)
    add_spectral_axis(axa, SPT_TEFF_BOUNDS, SPT_TEFF_CENTERS, SPT_LABELS)

    # ---- Panel (b): HZ distances — distance vs stellar mass ------------------
    _fm = np.isfinite(mass)   # drop stars with no Baraffe 5 Gyr point (F 7200 K)
    axb.fill_betweenx(mass[_fm], d_moist[_fm], d_max[_fm],
                      color=GRN, alpha=0.45, zorder=0)
    for dd, (_, col, _m) in zip((d_run, d_moist, d_max), BNDS):
        axb.plot(dd[_fm], mass[_fm], ls='-', color=col, lw=1.8)
    _km = _fm & (teff >= 2600.0)   # Kopparapu fit valid 2600-7200 K -> capped
    for bnd, col, _ in BNDS:                      # Kopparapu 2013 (dashed, per-star)
        axb.plot(np.sqrt(lum[_km] / kopp_seff(teff[_km], bnd)), mass[_km],
                 '--', color=col, lw=1.3)
    # Reference planetary systems: filled circles sized ~proportional to planet
    # radius (diameter = MS_E * R/R_Earth), at the host-star mass and each
    # planet's orbital distance.
    MS_E = 9.0          # marker diameter [pt] for a 1 R_Earth planet
    PLANET_C = '0.3'    # dark grey for all planet markers + labels
    # Solar system (Sun, 1 Msun): (a[AU], R[R_Earth]) — markers only, one group
    # label (individual planet names omitted).
    for _a, _r in [(0.387, 0.383), (0.723, 0.949), (1.000, 1.000), (1.524, 0.532)]:
        axb.plot(_a, 1.0, 'o', color=PLANET_C, ms=MS_E * _r, zorder=8)
    axb.text(1.00, 1.30, 'Solar System', ha='center', va='bottom',
             fontsize=9.5, color=PLANET_C, zorder=8)
    # TRAPPIST-1 (M = 0.0898 Msun, Teff 2566 K; Agol et al. 2021): planets b..h.
    # Its mass matches our M2600 case, so e/f/g fall in the computed HZ strip.
    T1M = 0.0898
    for _a, _r in [(0.01154, 1.116), (0.01580, 1.097), (0.02227, 0.788),
                   (0.02925, 0.920), (0.03853, 1.045), (0.04688, 1.129),
                   (0.06193, 0.755)]:
        axb.plot(_a, T1M, 'o', color=PLANET_C, ms=MS_E * _r, zorder=9)
    axb.text(0.0145, 0.069, 'TRAPPIST-1', ha='center', va='top',
             fontsize=9.5, color=PLANET_C, zorder=9)
    # Tidal-lock radius (Kasting, Whitmire & Reynolds 1993, Eq. 10, all-CGS):
    #   r_T = 0.027 (P0 t / Q)^(1/6) M^(1/3)  [cm],
    # with their adopted P0 = 13.5 hr, t = 4.5 Gyr, Q = 100 -> r_T(1 Msun)=0.459 AU
    # (Mercury, 0.39 AU, sits just inside it, as the paper notes).  Planets to the
    # LEFT of this dotted line are expected to be tidally locked.
    _P0 = 13.5 * 3600.0          # original rotation period [s]
    _TAGE = 4.5e9 * 3.15576e7    # 4.5 Gyr [s]
    _Q = 100.0
    _MSUN_G = 1.989e33           # solar mass [g]
    _AU_CM = 1.495979e13
    _mg = np.geomspace(0.055, 1.9, 100)
    _rT = 0.027 * (_P0 * _TAGE / _Q) ** (1.0 / 6.0) \
        * (_mg * _MSUN_G) ** (1.0 / 3.0) / _AU_CM
    axb.plot(_rT, _mg, ls=':', color='0.35', lw=1.5, zorder=3)
    axb.set_xscale('log'); axb.set_yscale('log')
    axb.set_xlim(0.009, 4.5)            # F-star outer edge reaches ~3.9 AU
    axb.set_ylim(0.055, 1.9)            # F-star mass = 1.61 Msun
    # On-curve label for the tidal-lock line, rotated to match it in display space.
    _il = int(np.argmin(np.abs(_mg - 0.22)))
    _p1 = axb.transData.transform((_rT[_il], _mg[_il]))
    _p2 = axb.transData.transform((_rT[_il + 6], _mg[_il + 6]))
    _ang = np.degrees(np.arctan2(_p2[1] - _p1[1], _p2[0] - _p1[0]))
    axb.text(_rT[_il], _mg[_il], 'Tidal lock radius  ', rotation=_ang,
             rotation_mode='anchor', ha='right', va='bottom',
             fontsize=9.5, color='0.35', zorder=4)
    axb.set_yticks([0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5])
    axb.set_yticklabels(['0.1', '0.2', '0.3', '0.5', '0.7', '1.0', '1.5'])
    axb.set_xticks([0.01, 0.03, 0.1, 0.3, 1.0, 3.0])
    axb.set_xticklabels(['0.01', '0.03', '0.1', '0.3', '1.0', '3.0'])
    add_spectral_axis(axb, SPT_MASS_BOUNDS, SPT_MASS_CENTERS, SPT_LABELS)
    axb.set_xlabel('Distance [AU]')
    axb.set_ylabel(r'Stellar mass  [$M_\odot$]')
    # HZ distances use the Baraffe et al. (1998) 5 Gyr isochrone (see caption/README).

    # Legend distinguishes source only (boundary identity is labelled on-curve);
    # placed on panel (b) upper-left, which is empty (high mass + small distance).
    handles = [Line2D([0], [0], color='0.35', lw=2, ls='-',  label='ExoColumn'),
               Line2D([0], [0], color='0.35', lw=2, ls='--', label='Clima (Kopparapu et al. 2013)')]
    axb.legend(handles=handles, fontsize=8, loc='upper left', frameon=False)

    fig.tight_layout()
    fig.savefig(FIG_PNG, dpi=300)
    fig.savefig(FIG_PDF)
    print(f"\nWrote: {FIG_PNG}\n       {FIG_PDF}")
    print("\n  Teff   M[Msun] L/Lsun  | S_eff run/moist/max  | d[AU] run/moist/max"
          "   (Kopp run/moist/max)")
    for k, t in enumerate(teff):
        kr, km, kx = (kopp_seff(t, b) for b in ('runaway', 'moist', 'maxgh'))
        print(f"  {t:5.0f}  {mass[k]:.3f}  {lum[k]:.4g} | "
              f"{s_run[k]:.3f}/{s_moist[k]:.3f}/{s_max[k]:.3f} | "
              f"{d_run[k]:.3f}/{d_moist[k]:.3f}/{d_max[k]:.3f}"
              f"   ({kr:.3f}/{km:.3f}/{kx:.3f})")


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
