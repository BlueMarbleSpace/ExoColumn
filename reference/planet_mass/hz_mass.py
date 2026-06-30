#!/usr/bin/env python3
"""
hz_mass.py  —  ExoColumn analogue of Kopparapu et al. (2014) Figure 3.

Kopparapu et al. (2014, ApJL 787, L29) extended the Kopparapu et al. (2013)
habitable-zone (HZ) limits to a range of planetary masses (0.1, 1, 5 M_Earth).
A more massive planet has higher surface gravity, which compresses the H2O
(inner edge) and CO2 (outer edge) column depths.  Their Figure 3 plots the HZ
limits versus stellar effective temperature for the three masses:

    inner edge = RUNAWAY greenhouse   (the Simpson-Nakajima OLR/Seff peak;
                 Kopparapu+2014 adopt this rather than the moist-greenhouse
                 limit, since the two differ by <2% and the runaway limit is
                 less sensitive to the assumed tropopause temperature)
    outer edge = MAXIMUM greenhouse   (the Seff minimum of a dense-CO2 sweep;
                 nearly mass-independent)

This script reproduces that figure from ExoColumn.  The mass enters through two
quantities, both following Kopparapu+2014's prescription:

  (1) SURFACE GRAVITY  g.  ExoRT bakes gravity in at COMPILE time (exo_g ->
      SHR_CONST_G; the radiation core's column amount is pdel/g).  So each mass
      needs its own binary, built via the Makefile's EXO_G override
      (make PVER=200 EXO_G=<g>).  The build/sweep loop is driven by
      sweep_all_masses.sh (or `python hz_mass.py all`).

  (2) BACKGROUND N2 PRESSURE  p_N2.  Scaled with mass per their Eq. (3),
      p_N2 ∝ R^2.40 ∝ M^0.75 (case 3, "N2 scaled with planet radius").  This is
      a runtime namelist value, passed to hz_inner/hz_outer.run_one(n2_bar=...).

MASS-RADIUS RELATION (Kopparapu+2014, from exoplanets.org):
    M/M_E = 0.968 (R/R_E)^3.2            (M < 5 M_E)
  =>  g(M)   = g_E * M^(1 - 2/3.2) = g_E * M^0.375
      p_N2(M) = 1 bar * M^(2.40/3.2)  = 1 bar * M^0.75
We ANCHOR the 1 M_E case to Earth exactly (g_E = 9.80616 m/s^2, p_N2 = 1 bar),
so the 1 M_E curve reproduces ExoColumn's validated inner/outer HZ reference
(reference/moist_runaway, reference/max_greenhouse) and the multi-stellar
Figure-6/7 Sun point.  Kopparapu's empirical fit coefficients (0.968, 0.937)
put g(1 M_E) ~2% off Earth; anchoring removes that while keeping their mass
TREND (the exact powers M^0.375, M^0.75 are independent of the normalisation).

  M [M_E]   g [m/s^2]   p_N2 [bar]
    0.1       4.135       0.178
    1.0       9.806       1.000
    5.0      17.933       3.344

The host-star set is the n68 cool/solar ladder (M 2600 K -> G Sun) shared with
hz_figure6/7, PLUS the hot F 7200 K endpoint on the n84 core (its SED exists
only at n84, which also resolves the strong F-star shortwave better; a separate
n84 build per mass, as in hz_add_n84_stars.py).  So each mass curve is the n68
M->G points + the n84 F point, spanning Kopparapu's full 2600-7200 K range.

USAGE (the build steps need the Intel OneAPI runtime; source setvars.sh first):
    source /opt/intel/oneapi/setvars.sh
    bash reference/planet_mass/sweep_all_masses.sh      # n68 sweep + n84 F + plot
  or, equivalently, from the project root:
    python reference/planet_mass/hz_mass.py all          # n68 M->G sweep (×3 mass)
    python reference/planet_mass/hz_mass.py addf         # + n84 F 7200 K (×3 mass)
  individual steps:
    python reference/planet_mass/hz_mass.py gravity 5   # print g for 5 M_E
    python reference/planet_mass/hz_mass.py sweep 5      # n68 sweep current binary
    python reference/planet_mass/hz_mass.py addf1 5      # merge F into 5 M_E cache
    python reference/planet_mass/hz_mass.py plot         # plot from caches

Per-mass results are cached to hz_mass_m{mass}.npz; `plot` (or HZ_REPLOT=1)
re-plots without re-running.
"""

import os
import sys
import subprocess
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
BUILD_DIR = os.path.join(ROOT, 'build')
FIG_PNG = os.path.join(HERE, 'hz_mass.png')
FIG_PDF = os.path.join(HERE, 'hz_mass.pdf')

PVER = int(os.environ.get('HZ_PVER', '200'))   # hz_inner/outer want PVER>=200

# --------------------------------------------------------------------------
# Kopparapu et al. (2014) planetary-mass parameterization.
# --------------------------------------------------------------------------
G_EARTH = 9.80616          # m/s^2
MR_COEFF, MR_EXP = 0.968, 3.2    # M/M_E = 0.968 (R/R_E)^3.2
MASSES = [0.1, 1.0, 5.0]   # M_Earth (Kopparapu+2014 range)


def radius_of_mass(m):
    """R/R_E from the Kopparapu+2014 mass-radius fit."""
    return (m / MR_COEFF) ** (1.0 / MR_EXP)


def gravity_of_mass(m):
    """Surface gravity [m/s^2], anchored so 1 M_E = Earth exactly.
    g(M)/g(1) = M / R(M)^2 / [1 / R(1)^2] = M^(1 - 2/3.2) = M^0.375."""
    return G_EARTH * m ** (1.0 - 2.0 / MR_EXP)


def n2_bar_of_mass(m):
    """N2 background partial pressure [bar], Kopparapu+2014 Eq. (3):
    p_N2 ∝ R^2.40 ∝ M^(2.40/3.2) = M^0.75, anchored to 1 bar at 1 M_E."""
    return m ** (2.40 / MR_EXP)


def mass_tag(m):
    """Filename/label tag for a mass (e.g. 0.1 -> '0.1', 5.0 -> '5.0')."""
    return f'{m:g}'


def cache_path(m):
    return os.path.join(HERE, f'hz_mass_m{mass_tag(m)}.npz')


# --------------------------------------------------------------------------
# Host-star set (label, n68 solar_file, Teff[K]).  '' => Sun (compile default).
# Index-/order-aligned with hz_figure6/7's n68 ladder.
# --------------------------------------------------------------------------
STARS = [
    ('M 2000 K',       'btsettl_T2000_g4.5_m0.0_n68.nc',    2000),
    ('M 2200 K',       'btsettl_T2200_g4.5_m0.0_n68.nc',    2200),
    ('M 2400 K',       'btsettl_T2400_g4.5_m0.0_n68.nc',    2400),
    ('M 2600 K',       'bt-settl_2600_logg4.5_FeH0_n68.nc', 2600),
    ('M 3000 K',       'bt-settl_3000_logg4.5_FeH0_n68.nc', 3000),
    ('M 3300 K',       'bt-settl_3300_logg4.5_FeH0_n68.nc', 3300),
    ('M 3700 K',       'bt-settl_3700_logg4.5_FeH0_n68.nc', 3700),
    ('K 4000 K',       'bt-settl_4000_logg4.5_FeH0_n68.nc', 4000),
    ('K 4500 K',       'bt-settl_4500_logg4.5_FeH0_n68.nc', 4500),
    ('K 4800 K',       'bt-settl_4800_logg4.5_FeH0_n68.nc', 4800),
    ('G 5780 K (Sun)', '',                                  5780),
]

# Hot F 7200 K endpoint — n84 core only.  The BT-Settl 7200 K SED exists at n84,
# which also resolves the strong F-star shortwave better than n68 (~2% in
# albedo/Seff; see reference/habitablezone/hz_add_n84_stars.py).  Added per mass
# on a separate n84 build via add_f_star() / the `addf` driver mode, then merged
# into the per-mass cache (so each curve = n68 M->G points + the n84 F point).
# A/B hotter stars give unphysical albedo>1 and are excluded.
F_STAR = ('F 7200 K', 'bt-settl_7200_logg4.5_FeH0_n84.nc', 7200)
N84_SPEC_DIR = '/models/ExoRT/source/src.n84equiv'

# Ultra-cool BT-Settl stars (n68 core/SED) extending the ladder below 2600 K.
# These are the leading entries of STARS; the `addcool` driver mode runs ONLY
# these and merges into existing per-mass caches (no full re-sweep), mirroring
# the F-star `addf` workflow.  The 2600 K point is unchanged (kept in STARS).
COOL_STARS = [
    ('M 2000 K', 'btsettl_T2000_g4.5_m0.0_n68.nc', 2000),
    ('M 2200 K', 'btsettl_T2200_g4.5_m0.0_n68.nc', 2200),
    ('M 2400 K', 'btsettl_T2400_g4.5_m0.0_n68.nc', 2400),
]

# Mass curve colours (Kopparapu+2014 convention: 0.1 blue, 1 green, 5 red).
MASS_COLOR = {0.1: '#1f77b4', 1.0: '#2ca02c', 5.0: '#d62728'}

# Main-sequence spectral-type bands on the Teff axis (same mapping as
# hz_figure7.py's top panel; topmost plotted type is F, A omitted).
SPT_TEFF_BOUNDS  = [3900., 5300., 6000.]
SPT_TEFF_CENTERS = [3100., 4600., 5650., 6650.]
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

# Runaway (Simpson-Nakajima) peak bracket: the inner-edge Seff(Ts) maximum sits
# at Ts ~ 300-330 K; sample 280-440 K finely enough to resolve it.
TS_RUNAWAY = np.arange(280.0, 440.0 + 1e-6,
                       float(os.environ.get('HZ_TS_STEP', '5')))
# Outer-edge CO2 grid [bar] (same physical endpoint as hz_outer: psat_CO2(273 K)).
PCO2_VALUES = np.geomspace(1.0, 34.7, int(os.environ.get('HZ_PCO2_N', '25')))

# Kopparapu et al. (2014) Table 1 -> Eq. (4): Seff = S0 + a T* + b T*^2 +
# c T*^3 + d T*^4, T* = Teff - 5780 K.  Runaway coeffs are mass-specific;
# maximum greenhouse is common to all masses.
KOPP2014 = {
    'runaway': {
        0.1: (0.99,  1.209e-4, 1.404e-8, -7.418e-12, -1.713e-15),
        1.0: (1.107, 1.332e-4, 1.58e-8,  -8.308e-12, -1.931e-15),
        5.0: (1.188, 1.433e-4, 1.707e-8, -8.968e-12, -2.084e-15),
    },
    'maxgh': (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
}


def kopp2014_seff(teff, coeffs):
    s0, a, b, c, dd = coeffs
    t = np.asarray(teff, float) - 5780.0
    return s0 + a * t + b * t**2 + c * t**3 + dd * t**4


def _load(modname, relpath):
    """Import a sibling reference script by file path (no sys.path edits)."""
    spec = importlib.util.spec_from_file_location(
        modname, os.path.join(ROOT, relpath))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------
# Build + sweep
# --------------------------------------------------------------------------
def build_binary(m):
    """Rebuild ExoColumn (default n68 core) at this mass's surface gravity
    (needs the Intel OneAPI runtime already on PATH / in the environment)."""
    g = gravity_of_mass(m)
    print(f"\n### BUILD  M = {m:g} M_E  ->  EXO_G = {g:.5f} m/s^2  "
          f"(PVER={PVER}) ###", flush=True)
    subprocess.run(['make', 'clean'], cwd=BUILD_DIR, check=True)
    subprocess.run(['make', f'PVER={PVER}', f'EXO_G={g:.5f}'],
                   cwd=BUILD_DIR, check=True)


def build_binary_n84(m):
    """Rebuild the n84 core at this mass's gravity, for the hot F 7200 K star
    (its SED and better-resolved shortwave need n84)."""
    g = gravity_of_mass(m)
    print(f"\n### BUILD (n84)  M = {m:g} M_E  ->  EXO_G = {g:.5f} m/s^2  "
          f"(PVER={PVER}) ###", flush=True)
    subprocess.run(['make', 'clean'], cwd=BUILD_DIR, check=True)
    subprocess.run(['make', f'SPEC_DIR={N84_SPEC_DIR}', f'PVER={PVER}',
                    f'EXO_G={g:.5f}'], cwd=BUILD_DIR, check=True)


def restore_earth():
    """Rebuild the default n68 core at Earth gravity so run/exocol.exe is the
    validated reference binary again."""
    print("\n### RESTORE default n68 Earth-gravity binary (EXO_G=9.80616) ###",
          flush=True)
    subprocess.run(['make', 'clean'], cwd=BUILD_DIR, check=True)
    subprocess.run(['make', f'PVER={PVER}', f'EXO_G={G_EARTH:.5f}'],
                   cwd=BUILD_DIR, check=True)


def add_f_star(m):
    """Run the F 7200 K star (n84 core/SED) for mass m on the CURRENTLY BUILT
    n84 binary and merge it into the per-mass cache (replacing any existing F).
    Same runaway/max-GH definitions as sweep_one_mass."""
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')
    n2 = n2_bar_of_mass(m)
    label, sf, te = F_STAR
    print(f"\n=== ADD F  M = {m:g} M_E   p_N2 = {n2:.4f} bar  [{sf}, n84] ===",
          flush=True)
    si = np.full(len(TS_RUNAWAY), np.nan)
    for j, ts in enumerate(TS_RUNAWAY):
        r = inner.run_one(ts, solar_file=sf, n2_bar=n2)
        if r is not None:
            si[j] = r['seff']
    run_pk = float(np.nanmax(si)) if np.any(np.isfinite(si)) else np.nan
    so = np.full(len(PCO2_VALUES), np.nan)
    for j, pc in enumerate(PCO2_VALUES):
        r = outer.run_one(pc, solar_file=sf, n2_bar=n2)
        if r is not None:
            so[j] = r['seff']
    max_gh = float(np.nanmin(so)) if np.any(np.isfinite(so)) else np.nan
    print(f"  {label:16s} Teff={te:5d} K   runaway Seff={run_pk:.3f}   "
          f"max-GH Seff={max_gh:.3f}", flush=True)

    if not os.path.exists(cache_path(m)):
        raise FileNotFoundError(
            f"No base cache {cache_path(m)} — run the n68 sweep first.")
    z = dict(np.load(cache_path(m), allow_pickle=True))
    lab = np.asarray(z['label'])
    keep = lab != label    # drop any previous F point before re-appending
    np.savez(cache_path(m),
             mass=z['mass'], g=z['g'], n2_bar=z['n2_bar'],
             teff=np.append(np.asarray(z['teff'])[keep], te),
             label=np.append(lab[keep], label),
             seff_runaway=np.append(np.asarray(z['seff_runaway'])[keep], run_pk),
             seff_maxgh=np.append(np.asarray(z['seff_maxgh'])[keep], max_gh))
    print(f"  merged F 7200 K into {cache_path(m)}", flush=True)


def add_cool_stars(m):
    """Run the ultra-cool stars (COOL_STARS, n68 core/SED) for mass m on the
    CURRENTLY BUILT n68 binary (which must be compiled at this mass's gravity)
    and merge them into the per-mass cache.  Same runaway/max-GH definitions as
    sweep_one_mass; idempotent (drops any existing entry with the same label)."""
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')
    n2 = n2_bar_of_mass(m)
    if not os.path.exists(cache_path(m)):
        raise FileNotFoundError(
            f"No base cache {cache_path(m)} — run the n68 sweep first.")
    z = dict(np.load(cache_path(m), allow_pickle=True))
    te_arr, lab_arr = np.asarray(z['teff']), np.asarray(z['label'])
    run_arr, max_arr = np.asarray(z['seff_runaway']), np.asarray(z['seff_maxgh'])
    for label, sf, te in COOL_STARS:
        print(f"\n=== ADD COOL  M = {m:g} M_E   p_N2 = {n2:.4f} bar  "
              f"[{sf}] ===", flush=True)
        si = np.full(len(TS_RUNAWAY), np.nan)
        for j, ts in enumerate(TS_RUNAWAY):
            r = inner.run_one(ts, solar_file=sf, n2_bar=n2)
            if r is not None:
                si[j] = r['seff']
        run_pk = float(np.nanmax(si)) if np.any(np.isfinite(si)) else np.nan
        so = np.full(len(PCO2_VALUES), np.nan)
        for j, pc in enumerate(PCO2_VALUES):
            r = outer.run_one(pc, solar_file=sf, n2_bar=n2)
            if r is not None:
                so[j] = r['seff']
        max_gh = float(np.nanmin(so)) if np.any(np.isfinite(so)) else np.nan
        print(f"  {label:16s} Teff={te:5d} K   runaway Seff={run_pk:.3f}   "
              f"max-GH Seff={max_gh:.3f}", flush=True)
        keep = lab_arr != label    # idempotent: drop a prior entry for this star
        te_arr = np.append(te_arr[keep], te)
        lab_arr = np.append(lab_arr[keep], label)
        run_arr = np.append(run_arr[keep], run_pk)
        max_arr = np.append(max_arr[keep], max_gh)
    np.savez(cache_path(m), mass=z['mass'], g=z['g'], n2_bar=z['n2_bar'],
             teff=te_arr, label=lab_arr,
             seff_runaway=run_arr, seff_maxgh=max_arr)
    print(f"  merged {len(COOL_STARS)} cool stars into {cache_path(m)}", flush=True)


def sweep_one_mass(m):
    """Run the inner (runaway) + outer (max-GH) multi-stellar sweep for mass m
    on the CURRENTLY BUILT binary (which must have been compiled at this mass's
    gravity).  Cache and return the per-star HZ limits."""
    inner = _load('hz_inner', 'reference/moist_runaway/hz_inner.py')
    outer = _load('hz_outer', 'reference/max_greenhouse/hz_outer.py')
    g, n2 = gravity_of_mass(m), n2_bar_of_mass(m)
    print(f"\n=== SWEEP  M = {m:g} M_E   g = {g:.4f} m/s^2   "
          f"p_N2 = {n2:.4f} bar ===", flush=True)

    teff, run_seff, max_seff, labels = [], [], [], []
    for label, sf, te in STARS:
        # INNER runaway = peak of Seff(Ts) over the Simpson-Nakajima bracket.
        si = np.full(len(TS_RUNAWAY), np.nan)
        for j, ts in enumerate(TS_RUNAWAY):
            r = inner.run_one(ts, solar_file=sf, n2_bar=n2)
            if r is not None:
                si[j] = r['seff']
        run_pk = float(np.nanmax(si)) if np.any(np.isfinite(si)) else np.nan
        # OUTER max greenhouse = minimum of Seff(pCO2).
        so = np.full(len(PCO2_VALUES), np.nan)
        for j, pc in enumerate(PCO2_VALUES):
            r = outer.run_one(pc, solar_file=sf, n2_bar=n2)
            if r is not None:
                so[j] = r['seff']
        max_gh = float(np.nanmin(so)) if np.any(np.isfinite(so)) else np.nan

        print(f"  {label:16s} Teff={te:5d} K   runaway Seff={run_pk:.3f}   "
              f"max-GH Seff={max_gh:.3f}", flush=True)
        teff.append(te); run_seff.append(run_pk)
        max_seff.append(max_gh); labels.append(label)

    out = dict(mass=m, g=g, n2_bar=n2,
               teff=np.array(teff, float), label=np.array(labels),
               seff_runaway=np.array(run_seff), seff_maxgh=np.array(max_seff))
    np.savez(cache_path(m), **out)
    print(f"  cached: {cache_path(m)}", flush=True)
    return out


# --------------------------------------------------------------------------
# Plot — Kopparapu et al. (2014) Figure 3 analogue
# --------------------------------------------------------------------------
def plot():
    from matplotlib.lines import Line2D
    caches = {m: cache_path(m) for m in MASSES if os.path.exists(cache_path(m))}
    if not caches:
        raise FileNotFoundError(
            "No hz_mass_m*.npz caches found — run the sweep first "
            "(python hz_mass.py all, or sweep_all_masses.sh).")

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    fig.patch.set_facecolor('white')

    # Smooth Teff grid for the dashed Kopparapu (2014) Table-1 overlay
    # (drawn across the full 2600-7200 K range the parametric fit covers).
    tk = np.linspace(2600.0, 7200.0, 200)

    data = {}
    for m in MASSES:
        if m not in caches:
            continue
        z = np.load(caches[m], allow_pickle=True)
        data[m] = z
        col = MASS_COLOR.get(m, '0.3')
        order = np.argsort(z['teff'])
        te = z['teff'][order]
        s_run = z['seff_runaway'][order]
        # ExoColumn runaway (inner edge): solid, mass-coloured.
        ax.plot(s_run, te, '-', color=col, lw=2.0, zorder=5,
                label=f'{m:g} ' + r'$M_\oplus$')
        # Kopparapu+2014 runaway (Table 1): dashed, same colour.
        ax.plot(kopp2014_seff(tk, KOPP2014['runaway'][m]), tk, '--',
                color=col, lw=1.2, zorder=4)

    # Outer edge (maximum greenhouse): nearly mass-independent.  Plot the
    # ExoColumn 1 M_E outer edge as the representative curve + Kopparapu dashed.
    m_ref = 1.0 if 1.0 in data else sorted(data)[0]
    zr = data[m_ref]
    order = np.argsort(zr['teff'])
    MAXGH_C = '#6a3d9a'    # purple, distinct from the 0.1 M_E inner-edge blue
    ax.plot(zr['seff_maxgh'][order], zr['teff'][order], '-',
            color=MAXGH_C, lw=2.0, zorder=5, label='Maximum greenhouse')
    ax.plot(kopp2014_seff(tk, KOPP2014['maxgh']), tk, '--',
            color=MAXGH_C, lw=1.2, zorder=4)

    # Solar-system reference points at the Sun's Teff (Seff = S/S0 = 1/a^2).
    for nm, a_au in [('Venus', 0.723), ('Earth', 1.000), ('Mars', 1.524)]:
        s = 1.0 / a_au**2
        ax.plot(s, 5780, 'o', color='0.25', ms=5, zorder=8)
        ax.annotate(nm, xy=(s, 5780), xytext=(0, 7), textcoords='offset points',
                    ha='center', va='bottom', fontsize=8, color='0.25', zorder=8)

    ax.set_xlim(1.40, 0.15)          # reversed; matches Fig. 7 top-panel S/S0 range
    ax.set_ylim(2000, 7200)          # extended to the 2000 K BT-Settl floor
    # NB: the dashed Kopparapu+2014 overlay (tk above) stays capped at 2600 K,
    # the stated lower validity bound of their parametric fit.
    ax.set_xlabel(r'Effective flux incident on the planet  $S/S_0$')
    ax.set_ylabel(r'Stellar effective temperature  $T_{\rm eff}$  [K]')
    add_spectral_axis(ax, SPT_TEFF_BOUNDS, SPT_TEFF_CENTERS, SPT_LABELS)

    # Legends: mass identity (colours) + source key (solid/dashed).
    leg1 = ax.legend(loc='lower left', fontsize=8.5, title='Inner edge (runaway)',
                     framealpha=0.92)
    leg1.get_title().set_fontsize(8.5)
    ax.add_artist(leg1)
    src = [Line2D([0], [0], color='0.35', lw=2, ls='-', label='ExoColumn'),
           Line2D([0], [0], color='0.35', lw=1.2, ls='--',
                  label='Kopparapu et al. 2014')]
    ax.legend(handles=src, loc='upper right', fontsize=8.5, framealpha=0.92)

    fig.tight_layout()
    fig.savefig(FIG_PNG, dpi=300)
    fig.savefig(FIG_PDF)
    print(f"\nWrote: {FIG_PNG}\n       {FIG_PDF}")

    # Console summary at the Sun (the canonical Kopparapu+2014 numbers).
    print("\n  Inner edge (runaway greenhouse) Seff at the Sun (Teff=5780 K):")
    print("    M [M_E]   ExoColumn   Kopparapu+2014")
    for m in MASSES:
        if m not in data:
            continue
        z = data[m]
        k = int(np.argmin(np.abs(z['teff'] - 5780)))
        kp = KOPP2014['runaway'][m][0]
        print(f"      {m:4g}      {z['seff_runaway'][k]:.3f}        {kp:.3f}")


# --------------------------------------------------------------------------
def main(argv):
    cmd = argv[1] if len(argv) > 1 else ('plot' if os.environ.get('HZ_REPLOT')
                                         == '1' else 'all')
    if cmd == 'gravity':
        m = float(argv[2])
        print(f'{gravity_of_mass(m):.5f}')
    elif cmd == 'n2':
        m = float(argv[2])
        print(f'{n2_bar_of_mass(m):.5f}')
    elif cmd == 'build':
        build_binary(float(argv[2]))
    elif cmd == 'sweep':
        sweep_one_mass(float(argv[2]))
    elif cmd == 'addf1':
        add_f_star(float(argv[2]))
    elif cmd == 'addcool1':
        add_cool_stars(float(argv[2]))
    elif cmd == 'plot':
        plot()
    elif cmd == 'all':
        for m in MASSES:
            build_binary(m)
            sweep_one_mass(m)
        restore_earth()
        plot()
    elif cmd == 'addf':
        # Append the hot F 7200 K endpoint (n84 core) to each per-mass cache,
        # then restore the default n68 Earth binary and re-plot.  Run AFTER the
        # n68 sweep (`all`) has produced the base caches.
        for m in MASSES:
            build_binary_n84(m)
            add_f_star(m)
        restore_earth()
        plot()
    elif cmd == 'addcool':
        # Append the ultra-cool stars (COOL_STARS, n68 core) to each per-mass
        # cache, then restore the default n68 Earth binary and re-plot.  Run
        # AFTER the n68 sweep (`all`) has produced the base caches.  Each mass is
        # rebuilt at its own gravity (compile-time), like `all`.
        for m in MASSES:
            build_binary(m)
            add_cool_stars(m)
        restore_earth()
        plot()
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == '__main__':
    main(sys.argv)
