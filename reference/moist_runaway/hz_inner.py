#!/usr/bin/env python3
"""
hz_inner.py  —  ExoColumn version of Kopparapu+2013 Figure 3 (inner HZ, G2V star).

REFERENCE CASE: moist-/runaway-greenhouse inner HZ, NON-IDEAL water EOS.
=========================================================================
This is the self-contained generator for the reference/moist_runaway case.
It uses the Kasting (1988) Appendix-A non-ideal moist pseudoadiabat with the
native IAPWS-95 EOS (h2o_eos='nonideal', which forces the Wagner-Pruss steam
saturation pressure).  The ideal-gas variant has been retired — this case is
the result we validate against.

To re-sweep (full ~14 min) or re-plot (instant), see README.md in this directory.
The full per-run configuration is embedded in NML_TEMPLATE below, so the sweep is
reproducible given an ExoColumn binary built at PVER>=200.

For each prescribed surface temperature Ts, ExoColumn builds a moist-adiabat
column (cold start: fully saturated at all levels including stratosphere,
isothermal stratosphere at t_strato = 200 K), calls ExoRT radiation once
(flux_only mode), and reads the output profiles and fluxes.

Four-panel figure matching Kopparapu+2013 Fig. 3:
  (a) OLR and absorbed SW vs Ts
  (b) Planetary albedo vs Ts
  (c) Effective stellar flux Seff = OLR / ASR vs Ts
  (d) H2O VMR vertical profiles for selected Ts values

Composition: N2=0.78 + O2=0.21 + Ar=0.01 + CO2=3.3e-4 + H2O (Kopparapu 2013 Earth tuning; no O3/CH4).

Shortwave solar zenith angle: 6-point Gauss-Legendre HEMISPHERIC quadrature
(sw_zenith_quad=.true., sw_nquad=6), matching Kopparapu+2013's 6-angle averaging,
rather than a single 60deg.  A single 60deg overweights grazing incidence and
biases the planetary albedo (panel b) HIGH by ~0.006-0.010; the hemispheric Bond
average removes that bias (see tools/diag_zenith_albedo.py).  The residual offset
vs Kopparapu (~0.015-0.02) is near-IR H2O shortwave absorption (ExoRT n68 vs their
HITEMP-2010 + BPS continuum) — a spectral-data difference, not a methodology one.

Model top: p_top = 0.002 Pa (lowered below the former 0.01 Pa) so that even the
COLDEST profiled column (280 K) reaches >=100 km — at 0.01 Pa it topped out at only
~96 km (the warm columns already reached ~110-150 km).  Panel (d) then clips at 100 km
(PANEL_ZTOP_KM) so every H2O profile spans the full altitude axis, as in Kopparapu
Fig 3d.  The shortwave is converged wrt the top (mass above 0.002 Pa is ~1e-8 of the
column; see tools/diag_ptop_sensitivity.py).  Build at PVER>=200 so the larger
log-pressure span does not coarsen the troposphere (preserves the panel a-c sawtooth
behaviour); the extra ~0.7 decade of (isothermal, radiatively quiet) stratosphere
above the cold point only marginally raises the sawtooth, which the median smoother
removes.

NOTE on the residual sawtooth in panels (a)-(c): this is a vertical-resolution
discretization artifact, not a physics error.  In the runaway regime OLR/ASR are set by
the emission/absorption level near the tropopause kink (moist adiabat meeting the
isothermal t_strato cap).  As Ts rises, ps grows (variable_ps) and the tropopause
migrates across the fixed log-spaced layers; each layer-boundary crossing re-samples that
kink and steps OLR/ASR by one tooth.  Confirmed by PVER doubling: 70->140 doubles the
number of teeth and roughly halves their amplitude (~1/N convergence).  The teeth track
the tropopause-layer index (k_top_conv) at 92%, and are NOT the cold-trap water (~5e-6,
radiatively negligible) nor the 10-bar k-table crossing.  See tools/diag_staircase*.py.
Independent caveat (see hz_roadmap): the ExoRT k-tables end at 500 K, so the high-Ts
(>~500 K) branch is extrapolated and should not be over-interpreted.
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# This script lives in reference/moist_runaway/, two levels below the project
# root (which holds run/, iofiles/, exocol_config.nml).  Outputs (figure + cache)
# are written next to this script.
HERE     = os.path.dirname(os.path.abspath(__file__))
ROOT     = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

# Kopparapu et al. (2013) reference data (provided directly by R. Kopparapu) for the
# direct model-vs-model comparison overlaid on the figure:
#   • waterloss_IHZ_present.dat — the inner-HZ sweep vs surface temperature TGO,
#     columns: TGO  SEFF  PALB  FH2O  FTIR(1)[=OLR]  FTSO(1)[=absorbed SW].
#     Feeds panels (a) (OLR/ASR), (b) (planetary albedo), (c) (Seff).
#   • clima_last.tab — one CLIMA water-loss vertical profile; col 1 = altitude [km],
#     col 4 = H2O volume mixing ratio.  Feeds panel (d).
KOPP_SWEEP = os.path.join(HERE, 'waterloss_IHZ_present.dat')
KOPP_CLIMA = os.path.join(HERE, 'clima_last.tab')

# ---------------------------------------------------------------------------
# Ts grid (env-overridable for quick tests: HZ_TS_MIN/HZ_TS_MAX/HZ_TS_STEP).
# The sweep spans the full inner-HZ runaway (Kasting 1988 Fig 1; Kopparapu 2013
# Fig 3): cold → moist greenhouse → runaway PLATEAU → critical point → supercritical
# runaway → POST-RUNAWAY RISE.
#   • Above the H2O critical point (Tc=647.1 K) no liquid ocean can exist, so the
#     whole water inventory is in the atmosphere and the cold start builds a DRY
#     (unsaturated) adiabat at the hot base, switching to the saturated pseudoadiabat
#     aloft where the rising parcel first cools to saturation (below Tc).  The water
#     column is then fixed at the inventory (~271 bar) for all Ts>Tc, so the cool
#     radiatively active upper atmosphere — and hence OLR — is Ts-INDEPENDENT: the
#     runaway plateau (~295 W/m², flat from ~400 to ~1600 K).
#   • At high Ts the Wien peak shifts into the near-IR/visible H2O WINDOWS, which
#     transmit the hot lower atmosphere to space → OLR rises steeply again (the
#     post-runaway branch; Kasting 1988, Kopparapu 2013 both show this — it is NOT
#     a flat plateau all the way up).  Our curve reproduces it (onset ~1800 K).
#   CAVEAT on the rise: ExoRT's n68 k-tables have a hardcoded T-grid ceiling of 500 K
#   (radgrid.F90, read-only); the deep >500 K layers use clamped opacity, so the
#   QUANTITATIVE steepness of the rise above ~1600 K (shaded) carries extrapolation
#   uncertainty (likely too steep — missing HITEMP hot bands).  Onset/shape are
#   physical; absolute high-Ts values would need HITEMP-class k-data.
TS_MIN  = float(os.environ.get('HZ_TS_MIN',  '200'))
TS_MAX  = float(os.environ.get('HZ_TS_MAX',  '2500'))
TS_STEP = float(os.environ.get('HZ_TS_STEP', '5'))
TS_VALUES = np.arange(TS_MIN, TS_MAX + TS_STEP, TS_STEP, dtype=float)   # K
H2O_TCRIT = 647.1   # K, water critical point (saturated → supercritical dry base)
# CAVEAT: the ExoRT n68 k-tables have a hardcoded T-grid ceiling of 500 K
# (read-only radgrid.F90); above ~1600 K the deep layers use clamped opacity, so
# the post-runaway OLR rise is qualitatively right but quantitatively extrapolated
# (would need HITEMP-class k-data).

# Water-vapour equation of state for the moist pseudoadiabat (FIXED for this
# reference case): Kasting (1988) Appendix-A general non-ideal pseudoadiabat with
# the native IAPWS-95 EOS (beta, real-gas entropies/cp), which forces the
# Wagner-Pruss steam saturation pressure.  The retired ideal-gas variant produced
# hz_inner.{pdf,png}; this case is hz_inner_nonideal.* (TAG kept for that lineage).
H2O_EOS = 'nonideal'
TAG     = '_nonideal'

ALBEDO   = 0.32      # surface albedo = Kopparapu et al. (2013) value, for a DIRECT
                     # comparison (their cloud-free 1-D model uses Ad=0.32 to mimic
                     # cloud reflection).  Surface albedo affects only absorbed SW
                     # (hence ASR & Seff), not the LW-driven cold-start T profile/OLR.
                     # [Our ExoRT-tuned Ts=288 K value was 0.24229; kept in git history.]
T_STRATO = 200.0     # isothermal stratosphere cap [K]

# Sub-freezing saturation phase.  'ice' — the model default — is ALSO CLIMA's
# actual convention: below 273.16 K Kopparapu's CLIMA saturates over ice
# everywhere (SATRAT: CC with the sublimation latent heat) and its sub-freezing
# moist adiabat (convec.f label 13) is the ice-saturated two-component
# pseudoadiabat that ExoColumn reproduces via exocol_convadj::twocomp_dlnTdlnP
# (sub-273 K T(P) agrees with clima_last.tab to ±0.07 K at Ts=380 K).  The
# earlier 'liquid' setting here was based on a misreading of their convention;
# its apparent hot-case agreement came from two compensating errors (wetter
# supercooled-liquid esat vs the too-steep textbook malr fallback).
# Env-overridable: HZ_COLD_TRAP=liquid for the sensitivity variant.
COLD_TRAP_PHASE = os.environ.get('HZ_COLD_TRAP', 'ice')

MW_H2O   = 18.015    # g/mol
P_DRY    = 1.0e5     # Pa; the sweep's dry-N2 inventory (cold-start DEFAULT_PS).
                     # With variable_ps the model surface pressure is
                     # ps = P_DRY + e_sfc, so e_sfc = ps − P_DRY exactly.

# The OLR/albedo/Seff curves carry a small high-frequency SAWTOOTH: as Ts sweeps,
# the tropopause kink (moist adiabat meeting the isothermal t_strato) migrates
# across the fixed log-spaced layers and the RT re-samples it at each crossing.
# It is a vertical-resolution sampling artifact, NOT physics (PVER 70→140 halves
# its amplitude — see project_hz_staircase).  A light centered running-median
# (window HZ_SMOOTH points, default 5 ≈ 25 K) suppresses it while preserving the
# physical shape (plateau, albedo dip, post-runaway rise); set HZ_SMOOTH=1 for the
# raw curve.  Higher PVER halves the sawtooth at the source (1/N convergence:
# 0.008→0.004 p-p in albedo); the window-5 median then removes the residual.
# The published figure is generated at PVER>=200 (build with `make clean &&
# make PVER=200`): the extended p_top=0.01 Pa top spreads the log-pressure grid
# over ~2 more decades, so PVER must rise from the former 140 to ~200 to keep the
# tropospheric resolution (hence the sawtooth amplitude) unchanged.
SMOOTH_WIN = int(os.environ.get('HZ_SMOOTH', '5'))

MOIST_GH_VMR = 3.0e-3   # moist-greenhouse stratospheric H2O VMR threshold (Kopparapu Sec 3.1)

# Ts values at which to save the full H2O vertical profile for panel (d)
PROFILE_TS = [280.0, 300.0, 320.0, 340.0, 360.0, 380.0]

# Of the profiled curves, only these get an on-curve label in panel (d); the
# others follow the same monotonic ordering and are left unlabelled to reduce
# clutter.  LABEL_TS_PREFIX gets the explicit "Ts = " prefix so the reader knows
# the numbers are surface temperatures (the rest are bare "NNN K").
LABEL_TS = [320.0, 340.0, 360.0]
LABEL_TS_PREFIX = 340.0

# Altitude clip for panel (d) [km].  p_top=0.002 Pa guarantees every profiled column
# reaches this; the panel is capped here so the warm (taller) columns do not stretch
# the axis and every curve fills it (cf. Kopparapu Fig 3d).
PANEL_ZTOP_KM = 100.0

# Results cache so the figure can be re-plotted (label tweaks, clip, styling) without
# re-running the ~14-min sweep.  Set HZ_REPLOT=1 to load this and skip the runs.
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'hz_inner{TAG}.npz')

# BPS H2O-continuum overlay.  The Wolf BPS continuum is the one Kopparapu's CLIMA
# used; ExoColumn defaults to MT_CKD 3.3 (8-gpt).  Overlaying both isolates the
# continuum's contribution to the residual albedo/Seff offset vs Kopparapu (see
# memory project_bps_continuum: BPS closes ~half the albedo gap; OLR is continuum-
# insensitive, so the OLR/Seff residual is line data).  ON by default; HZ_BPS=0
# reproduces the MT_CKD-only figure.  In flux_only mode the continuum changes ONLY
# the radiation diagnostics (OLR/ASR/albedo/Seff) — the cold-start T and H2O
# profiles are identical for both continua — so the overlay appears on panels
# (a)-(c) only; panel (d) is continuum-independent.
BPS_OVERLAY = os.environ.get('HZ_BPS', '1') != '0'
CACHE_BPS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         f'hz_inner{TAG}_bps.npz')

NML_TEMPLATE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  ihz_profile    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  h2o_eos        = '{h2o_eos}'
  h2o_continuum  = '{continuum}'
  cold_trap_phase = '{cold_trap}'
  sw_zenith_quad = .true.
  sw_nquad       = 6
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = {t_strato:.1f}
  p_top      = 0.002
  rh_init    = 1.0
  coszrs     = 0.5
  asdir      = {albedo:.4f}
  asdif      = {albedo:.4f}
  aldir      = {albedo:.4f}
  aldif      = {albedo:.4f}
/
&exocol_composition
  n2_vmr   = 0.99967
  o2_vmr   = 0.0
  ar_vmr   = 0.0
  co2_vmr  = 3.3e-4
  ch4_vmr  = 0.0
  o3_vmr   = 0.0
/
"""
# Background composition: pure N2 (1 bar) + 330 ppm CO2 + saturated H2O — matching
# Kopparapu et al. (2013)'s IHZ "Earth" model (their Fig 3: FCO2=3.3e-4, N2 background;
# no O2/O3/CH4) and consistent with the OHZ sweep.  O2+Ar (the earlier Earth-air
# mix) shifted Seff by only ~+0.003 and albedo by ~+0.005, but pure N2 is the
# apples-to-apples choice and isolates the H2O line-data residual.

# ---------------------------------------------------------------------------

def run_one(ts, continuum='mtckd'):
    """
    Run ExoColumn flux_only at surface temperature ts with the chosen H2O
    continuum ('mtckd' | 'bps').
    Returns dict with scalar diagnostics and profiles, or None on failure.
    """
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO,
                              h2o_eos=H2O_EOS, continuum=continuum,
                              cold_trap=COLD_TRAP_PHASE)
    orig = None
    if os.path.exists(NML_PATH):
        with open(NML_PATH) as f:
            orig = f.read()
    try:
        with open(NML_PATH, 'w') as f:
            f.write(nml)
        result = subprocess.run([EXE], cwd=ROOT, capture_output=True,
                                text=True, timeout=180)
        if result.returncode != 0:
            print(f"  FAIL Ts={ts:.0f} K  rc={result.returncode}")
            print(result.stderr[-300:] if result.stderr else '')
            return None

        with nc.Dataset(OUT_NC) as ds:
            lwup   = ds['LWUP'][:]
            swup   = ds['SWUP'][:]
            swdn   = ds['SWDN'][:]
            tmid   = ds['tmid'][:]
            pmid   = ds['pmid'][:]
            h2ommr = np.array(ds['h2ommr'][:])
            mwdry  = float(ds['mw'][:])
            zint   = np.array(ds['zint'][:])   # interface heights [m]
            ps     = float(ds['ps'][:])

        olr   = float(lwup[0])
        asr   = float(swdn[0] - swup[0])
        alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
        seff  = olr / asr if asr > 1e-6 else np.nan

        # H2O volume mixing ratio profile — exact two-component conversion from
        # the moist specific humidity q.  The dilute form q·mwdry/18 overshoots
        # by ~25% at the 380 K surface (q ≈ 0.45); CLIMA's FH2O is a true mole
        # fraction, so the exact form is required for the panel-(d) overlay.
        h2o_vmr = (h2ommr / MW_H2O) / (h2ommr / MW_H2O + (1.0 - h2ommr) / mwdry)
        # Stratospheric H2O VMR = the cold-trap minimum (driest point in the
        # column); used to locate the moist-greenhouse / water-loss limit
        # (the Ts where this reaches MOIST_GH_VMR = 3e-3; Kopparapu Sec 3.1).
        strat_vmr = float(np.nanmin(h2o_vmr)) if np.any(np.isfinite(h2o_vmr)) else np.nan
        # Layer-midpoint altitude [km] from interface heights (TOA=0 … srf=pverp)
        zmid_km = 0.5 * (zint[:-1] + zint[1:]) / 1000.
        # Prepend the true surface point (z = 0): with variable_ps the column is
        # p_dry (1 bar) + e_sfc, so the surface H2O mole fraction is just
        # (ps − p_dry)/ps.  CLIMA's profiles include z = 0; without this point
        # our plotted curves start at the lowest midpoint (~200 m), which at
        # Ts = 280–300 K reads ~5–10% dry of the surface value.
        sfc_vmr = (ps - P_DRY) / ps
        if 0.0 < sfc_vmr < 1.0:
            zmid_km = np.concatenate([zmid_km, [0.0]])
            h2o_vmr = np.concatenate([h2o_vmr, [sfc_vmr]])
            pmid    = np.concatenate([np.array(pmid), [ps]])

        return dict(olr=olr, asr=asr, alpha=alpha, seff=seff, strat_vmr=strat_vmr,
                    h2o_vmr=h2o_vmr, pmid=np.array(pmid),
                    tmid=np.array(tmid), zmid_km=zmid_km)
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def _save_cache(arr, profiles, path=CACHE):
    """Persist sweep results so the figure can be re-plotted without re-running.
    `profiles` may be empty (the BPS overlay reuses the MT_CKD profiles for panel d)."""
    ks = sorted(profiles)
    np.savez(path, rows=arr,
             prof_ts=np.array(ks, dtype=float),
             prof_z=np.array([profiles[k][0] for k in ks]),
             prof_v=np.array([profiles[k][1] for k in ks]),
             prof_p=np.array([profiles[k][2] for k in ks]))
    print(f"Cached: {path}")


def _load_cache(path=CACHE):
    d = np.load(path)
    arr = d['rows']
    profiles = {float(t): (z, v, p) for t, z, v, p in
                zip(d['prof_ts'], d['prof_z'], d['prof_v'], d['prof_p'])}
    return arr, profiles


def _sweep(continuum):
    """Run the full Ts sweep for one continuum.  Returns (rows_array, profiles).
    Profiles (for panel d) are only meaningful for the first/MT_CKD pass."""
    rows = []
    profiles = {}
    for ts in TS_VALUES:
        r = run_one(ts, continuum)
        if r is None:
            print(f"  [{continuum}] Ts={ts:5.1f} K  SKIPPED")
            continue
        print(f"  [{continuum}] Ts={ts:5.1f} K  OLR={r['olr']:7.2f}  ASR={r['asr']:7.2f}"
              f"  α={r['alpha']:.3f}  Seff={r['seff']:.4f}")
        rows.append((ts, r['olr'], r['asr'], r['alpha'], r['seff'], r['strat_vmr']))
        if any(abs(ts - ts_p) < 0.5 for ts_p in PROFILE_TS):
            profiles[ts] = (r['zmid_km'], r['h2o_vmr'], r['pmid'])
    return np.array(rows), profiles


def main():
    # Re-plot from the cached sweep (instant; for label/clip/styling tweaks).
    if os.environ.get('HZ_REPLOT'):
        if not os.path.exists(CACHE):
            raise FileNotFoundError(f"No cache to re-plot: {CACHE} (run the sweep first)")
        print(f"HZ_REPLOT: re-plotting from {CACHE}")
        arr, profiles = _load_cache()
        ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a = arr.T
        bps = None
        if BPS_OVERLAY and os.path.exists(CACHE_BPS):
            bps, _ = _load_cache(CACHE_BPS)
            print(f"HZ_REPLOT: BPS overlay from {CACHE_BPS}")
        elif BPS_OVERLAY:
            print("  (no BPS cache yet — plotting MT_CKD only; run the sweep to build it)")
        _plot(ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a, profiles, bps=bps)
        return

    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"Executable not found: {EXE}")

    print("ExoColumn inner-HZ sweep  —  Figure 3 (Kopparapu+2013 style)")
    print(f"  H2O EOS     : {H2O_EOS}" +
          ("  (Kasting 1988 Appendix-A, IAPWS-95)" if H2O_EOS == 'nonideal' else ""))
    print(f"  H2O contin. : MT_CKD" + ("  +  BPS overlay" if BPS_OVERLAY else ""))
    print(f"  cold trap   : saturate over {COLD_TRAP_PHASE}" +
          ("  (= CLIMA convention; see COLD_TRAP_PHASE note)" if COLD_TRAP_PHASE == 'ice'
           else "  (sensitivity variant; CLIMA + model default = ice)"))
    print(f"  Ts range    : {TS_VALUES[0]:.0f}–{TS_VALUES[-1]:.0f} K")
    print(f"  t_strato    : {T_STRATO:.0f} K  (isothermal, fully saturated)")
    print(f"  variable_ps : ps = p_N2(1 bar) + esat(Ts)")
    print(f"  composition : N2 + CO2=3.3e-4 + H2O  [no O2/Ar/O3/CH4]")
    print()

    arr, profiles = _sweep('mtckd')
    if arr.size == 0:
        print("No successful runs.")
        return
    _save_cache(arr, profiles, CACHE)
    ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a = arr.T

    bps = None
    if BPS_OVERLAY:
        print("\n  --- BPS continuum pass (radiation diagnostics only) ---")
        barr, _ = _sweep('bps')
        if barr.size:
            _save_cache(barr, {}, CACHE_BPS)   # profiles are continuum-independent
            bps = barr

    _plot(ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a, profiles, bps=bps)


def _smooth(y, w=SMOOTH_WIN):
    """Centered running median (window w points) to suppress the vertical-
    resolution sawtooth (a numerical sampling artifact, not physics).  Preserves
    monotonic trends and broad features; w<=1 returns the raw curve."""
    y = np.asarray(y, float)
    if w <= 1 or y.size < 3:
        return y
    h = w // 2
    out = y.copy()
    for i in range(y.size):
        out[i] = np.median(y[max(0, i - h):min(y.size, i + h + 1)])
    return out


def _load_kopp_sweep():
    """Kopparapu+2013 inner-HZ sweep (waterloss_IHZ_present.dat).
    Returns dict of columns or None if the file is absent."""
    if not os.path.isfile(KOPP_SWEEP):
        return None
    # 3-line header (title, column names, '====' rule) then numeric rows.
    d = np.genfromtxt(KOPP_SWEEP, skip_header=3)
    return dict(tgo=d[:, 0], seff=d[:, 1], palb=d[:, 2], fh2o=d[:, 3],
                olr=d[:, 4], asr=d[:, 5])


def _load_kopp_clima():
    """Kopparapu CLIMA water-loss vertical profiles (clima_last.tab).
    The file concatenates one ALT/P/T/FH2O/... block PER surface temperature
    (220, 240, ... K); each block runs top-of-atmosphere -> surface.  Columns are
    ALT[km]  P[bar]  T[K]  FH2O[VMR]  ...  Returns {surface_T: (altitude_km,
    h2o_vmr)} keyed by the bottom-level (alt=0) temperature, or None if absent."""
    if not os.path.isfile(KOPP_CLIMA):
        return None
    with open(KOPP_CLIMA) as f:
        lines = f.readlines()
    hdr = [i for i, ln in enumerate(lines) if 'ALT' in ln] + [len(lines)]
    blocks = {}
    for i0, i1 in zip(hdr[:-1], hdr[1:]):
        d = np.genfromtxt(lines[i0 + 1:i1])
        if d.ndim != 2 or d.shape[0] < 2:
            continue
        surf_t = round(float(d[-1, 2]))   # T at the alt=0 (surface) row
        blocks[surf_t] = (d[:, 0], d[:, 3])
    return blocks


def _plot(ts, olr, asr, alpha, seff, strat_vmr, profiles, bps=None):
    # Suppress the vertical-resolution sawtooth (see SMOOTH_WIN note above).
    olr, asr, alpha, seff = (_smooth(olr), _smooth(asr),
                             _smooth(alpha), _smooth(seff))

    # BPS-continuum overlay (radiation diagnostics only; same Ts grid).
    bts = bolr = basr = balp = bseff = None
    if bps is not None and bps.size:
        bts, bolr, basr, balp, bseff = (bps[:, 0], _smooth(bps[:, 1]),
                                        _smooth(bps[:, 2]), _smooth(bps[:, 3]),
                                        _smooth(bps[:, 4]))

    # ExoColumn's own moist-greenhouse (water-loss) and runaway-greenhouse Seff:
    #   moist GH = Seff at the first Ts where the stratospheric H2O VMR reaches
    #             MOIST_GH_VMR (3e-3, the water-loss threshold; Kopparapu Sec 3.1).
    #   runaway  = Seff at the H2O critical temperature Tc = 647.1 K — above Tc no
    #             liquid ocean can survive, so this is the canonical runaway point
    #             (the plateau is flat through here, ~the runaway-greenhouse limit).
    svmr = _smooth(strat_vmr)
    exo_moist = exo_runaway = moist_ts = runaway_ts = np.nan
    _wet = np.where(svmr >= MOIST_GH_VMR)[0]
    if _wet.size:
        exo_moist = float(seff[_wet[0]]); moist_ts = float(ts[_wet[0]])
    if ts[0] <= H2O_TCRIT <= ts[-1]:
        runaway_ts  = H2O_TCRIT
        exo_runaway = float(np.interp(H2O_TCRIT, ts, seff))
    print(f"  ExoColumn moist-GH Seff = {exo_moist:.3f} (Ts={moist_ts:.0f} K)"
          f" ;  runaway Seff = {exo_runaway:.3f} (Tc={runaway_ts:.0f} K)")

    # BPS-continuum greenhouse-limit Seff.  The cold-start T/H2O profiles (hence
    # strat_vmr → moist_ts, and Tc) are continuum-independent, so the BPS limits
    # are simply the BPS Seff curve read at the SAME Ts as the MT_CKD limits.
    bps_moist = bps_runaway = np.nan
    if bseff is not None:
        if np.isfinite(moist_ts):
            bps_moist = float(np.interp(moist_ts, bts, bseff))
        if np.isfinite(runaway_ts):
            bps_runaway = float(np.interp(runaway_ts, bts, bseff))
        print(f"  BPS       moist-GH Seff = {bps_moist:.3f} (Ts={moist_ts:.0f} K)"
              f" ;  runaway Seff = {bps_runaway:.3f} (Tc={runaway_ts:.0f} K)")

    # Kopparapu+2013 reference data for the model-vs-model overlay (dashed grey).
    kopp = _load_kopp_sweep()
    kclima = _load_kopp_clima()
    GREY    = '0.3'         # single neutral grey for all non-data annotation
                            # elements (markers, reference line, legend keys, labels)
    KOPP_KW = dict(lw=0.9, ls='--', zorder=2)   # colour set per panel (matches the
                                                # ExoColumn quantity); thin dashed

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.0), dpi=300)
    fig.patch.set_facecolor('white')
    (ax_a, ax_b), (ax_c, ax_d) = axes

    kw = dict(xlim=(200., 2200.))

    # --- (a) OLR and ASR vs Ts ---
    BPS_KW = dict(lw=1.1, ls=':', zorder=3)   # dotted = ExoColumn BPS continuum
    ax_a.plot(ts, olr, color='C3', lw=1.5)
    ax_a.plot(ts, asr, color='C0', lw=1.5)
    if bolr is not None:
        ax_a.plot(bts, bolr, color='C3', **BPS_KW)
        ax_a.plot(bts, basr, color='C0', **BPS_KW)
    if kopp is not None:
        ax_a.plot(kopp['tgo'], kopp['olr'], color='C3', **KOPP_KW)
        ax_a.plot(kopp['tgo'], kopp['asr'], color='C0', **KOPP_KW)
    ax_a.set_ylabel('Flux (W m⁻²)')
    ax_a.set_xlabel('Surface temperature (K)')
    ax_a.set_ylim(100., 450.)      # focus on the plateau + onset of the rise
    # Label the curves directly (no legend): F_IR = outgoing longwave (OLR),
    # F_SOL = absorbed shortwave.  Anchored on the runaway plateau.
    ax_a.text(800., 315., '$F_\\mathrm{IR}$', color='C3', fontsize=10,
              fontweight='bold', ha='center', va='bottom')
    ax_a.text(800., 248., '$F_\\mathrm{SOL}$', color='C0', fontsize=10,
              fontweight='bold', ha='center', va='top')
    # Single style legend (applies to every panel): solid = this work, dashed
    # = Kopparapu+2013.  Colour continues to encode the quantity in each panel.
    if kopp is not None:
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], color=GREY, lw=1.5,
                          label='ExoColumn (MT_CKD)' if bolr is not None else 'ExoColumn')]
        if bolr is not None:
            handles.append(Line2D([0], [0], color=GREY, lw=1.1, ls=':',
                                  label='ExoColumn (BPS)'))
        handles.append(Line2D([0], [0], color=GREY, lw=0.9, ls='--',
                              label='Clima (Kopparapu et al. 2013)'))
        ax_a.legend(handles=handles, loc='upper left', fontsize=8, frameon=False)
    ax_a.set(**kw)
    ax_a.set_facecolor('white')

    # --- (b) Planetary albedo vs Ts ---
    ax_b.plot(ts, alpha, color='C2', lw=1.5)
    if balp is not None:
        ax_b.plot(bts, balp, color='C2', **BPS_KW)
    if kopp is not None:
        ax_b.plot(kopp['tgo'], kopp['palb'], color='C2', **KOPP_KW)
    ax_b.axhline(ALBEDO, color=GREY, lw=0.9, ls=':')
    ax_b.set_ylabel('Planetary albedo $\\alpha_p$')
    ax_b.set_xlabel('Surface temperature (K)')
    ax_b.set_ylim(0.15, 0.35)
    # Label the curve and the reference line directly (no legend).
    ax_b.text(1400., ALBEDO - 0.006, f'Surface albedo = {ALBEDO:.2f}',
              color=GREY, fontsize=8, ha='center', va='top')
    ax_b.set(**kw)
    ax_b.set_facecolor('white')

    # --- (c) Seff vs Ts ---
    ax_c.plot(ts, seff, color='C1', lw=1.5)
    if bseff is not None:
        ax_c.plot(bts, bseff, color='C1', **BPS_KW)
    if kopp is not None:
        ax_c.plot(kopp['tgo'], kopp['seff'], color='C1', **KOPP_KW)
    # Markers on the curve at the two greenhouse limits: moist-GH (water-loss; the
    # Ts where stratospheric H2O VMR first reaches 3e-3) and runaway-GH (at the
    # H2O critical temperature Tc).  Circle = MT_CKD, diamond = BPS, square = Clima.
    if np.isfinite(exo_moist):
        ax_c.plot(moist_ts, exo_moist, 'o', color=GREY, ms=3,
                  markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)
    if np.isfinite(exo_runaway):
        ax_c.plot(runaway_ts, exo_runaway, 'o', color=GREY, ms=3,
                  markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)
    if np.isfinite(bps_moist):
        ax_c.plot(moist_ts, bps_moist, 'D', color=GREY, ms=3,
                  markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)
    if np.isfinite(bps_runaway):
        ax_c.plot(runaway_ts, bps_runaway, 'D', color=GREY, ms=3,
                  markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)

    # Kopparapu+2013 IHZ limits on their own (dashed) Seff curve, derived with the
    # same physical definitions used for ExoColumn: moist greenhouse where the
    # stratospheric H2O VMR (FH2O column) first reaches MOIST_GH_VMR (interpolated
    # on Kopparapu's coarser 20 K grid), runaway greenhouse at the H2O critical
    # temperature Tc.  Drawn as grey SQUARES (ExoColumn uses circles); each also
    # sits on the dashed-vs-solid curve of its own model.
    kopp_moist_seff = kopp_run_seff = np.nan
    if kopp is not None:
        kt, ksf, kfh = kopp['tgo'], kopp['seff'], kopp['fh2o']
        wet = np.where(kfh >= MOIST_GH_VMR)[0]
        if wet.size and wet[0] > 0:
            j = wet[0]
            kopp_moist_ts = float(np.interp(np.log10(MOIST_GH_VMR),
                                            np.log10(kfh[j - 1:j + 1]), kt[j - 1:j + 1]))
            kopp_moist_seff = float(np.interp(kopp_moist_ts, kt, ksf))
            ax_c.plot(kopp_moist_ts, kopp_moist_seff, 's', color=GREY, ms=3,
                      markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)
        if kt[0] <= H2O_TCRIT <= kt[-1]:
            kopp_run_seff = float(np.interp(H2O_TCRIT, kt, ksf))
            ax_c.plot(H2O_TCRIT, kopp_run_seff, 's', color=GREY, ms=3,
                      markeredgecolor=GREY, markeredgewidth=0.5, zorder=6)

    ax_c.set_ylabel('$S_{\\rm eff}$')
    ax_c.set_xlabel('Surface temperature (K)')
    ax_c.set_ylim(0.4, 1.8)
    # Label each limit (with the implied IHZ orbital distance d = 1/sqrt(Seff) [AU]),
    # giving both models' Seff: ExoColumn (circle, solid curve) and Kopparapu+2013
    # (square, dashed curve).  The two limits are close in Ts (~350 vs ~647 K), so
    # moist sits above its marker (arrow) and runaway just to the right of its own.
    # Lead each model's line with its plot-marker glyph (● circle = ExoColumn,
    # ■ square = Clima) so the label doubles as the marker key.  Same grey as the
    # plotted markers and text.
    CIRCLE, DIAMOND, SQUARE = '•', '◆', '▪'   # glyphs ~matching the ms=3 markers
    def _limit_text(name, exo_s, bps_s, kopp_s):
        def _au(s):
            return f'$S_\\mathrm{{eff}}$ = {s:.3f} ({1.0 / np.sqrt(s):.3f} AU)'
        lines = [name, f'{CIRCLE} ExoColumn (MT_CKD):  {_au(exo_s)}']
        if np.isfinite(bps_s):
            lines.append(f'{DIAMOND} ExoColumn (BPS):       {_au(bps_s)}')
        if np.isfinite(kopp_s):
            lines.append(f'{SQUARE} Clima:                        {_au(kopp_s)}')
        return '\n'.join(lines)

    if np.isfinite(exo_moist):
        ax_c.annotate(_limit_text('MOIST GREENHOUSE', exo_moist, bps_moist, kopp_moist_seff),
                      xy=(moist_ts, exo_moist), xytext=(2, 66),
                      textcoords='offset points', color=GREY, fontsize=8,
                      ha='left', va='top',
                      arrowprops=dict(arrowstyle='->', color=GREY, lw=0.6,
                                      shrinkB=3))
    if np.isfinite(exo_runaway):
        ax_c.annotate(_limit_text('RUNAWAY GREENHOUSE', exo_runaway, bps_runaway, kopp_run_seff),
                      xy=(runaway_ts, exo_runaway), xytext=(40, -48),
                      textcoords='offset points', color=GREY, fontsize=8,
                      ha='left', va='top',
                      arrowprops=dict(arrowstyle='->', color=GREY, lw=0.6,
                                      shrinkB=3))
    ax_c.set(**kw)
    ax_c.set_facecolor('white')

    # --- (d) H2O VMR vertical profiles vs altitude for selected Ts values ---
    items_d = sorted(profiles.items())
    colors_d = plt.cm.plasma(np.linspace(0.05, 0.85, max(len(items_d), 1)))
    lab = []   # (ts_p, x_anchor, color) for de-collided labelling below
    for (ts_p, (zmid, vmr, pmid)), col in zip(items_d, colors_d):
        # Drop only genuine steam-cap layers (vmr -> ~1 where es = 0.99*p); the
        # cold-trap stratosphere is a well-behaved CONSTANT vmr that we keep all the
        # way to the model top.  (The old pmid>0.03 Pa cut chopped ~15 km of good
        # constant-vmr strat off the cold columns, stopping them ~90 km short of the
        # clip.)  With p_top = 0.002 Pa every column tops out >100 km, so the 100 km
        # axis clip frames all curves and hides any cap layer above it (cf. Kopp 3d).
        mask = vmr < 0.99
        z = zmid[mask]; v = vmr[mask]
        ax_d.semilogx(v, z, color=col, lw=1.5)
        # Overlay the Kopparapu CLIMA profile at the SAME surface temperature,
        # dashed and colour-matched (solid = ExoColumn, dashed = Kopparapu).
        if kclima is not None:
            kc = kclima.get(round(ts_p))
            if kc is not None:
                kz, kv = kc
                km = kz <= PANEL_ZTOP_KM
                ax_d.semilogx(kv[km], kz[km], color=col, lw=0.9, ls='--')
        # Label only the selected curves (the rest follow the same pattern).
        if not any(abs(ts_p - t) < 0.5 for t in LABEL_TS):
            continue
        # Anchor the label on the curve's (vertical, constant-vmr) stratospheric
        # segment, at the highest point still inside the clip.
        vis = np.where(z <= PANEL_ZTOP_KM)[0]
        ilbl = vis[np.argmax(z[vis])] if vis.size else int(np.argmax(z))
        lab.append((ts_p, v[ilbl], col))

    # Place labels just below the top edge; nudge any pair too close in VMR down
    # in altitude so they read (the selected 320/340/360 K are well separated).
    order = sorted(range(len(lab)), key=lambda j: lab[j][1])
    ylab = {}
    last_x = last_y = None
    for j in order:
        x = lab[j][1]
        y = PANEL_ZTOP_KM - 2.0
        if last_x is not None and abs(np.log10(x / last_x)) < 0.45 and abs(y - last_y) < 9.0:
            y = last_y - 11.0
        ylab[j] = y
        last_x, last_y = x, y
    for j, (ts_p, x, col) in enumerate(lab):
        txt = (f'$T_s$ = {ts_p:.0f} K' if abs(ts_p - LABEL_TS_PREFIX) < 0.5
               else f'{ts_p:.0f} K')
        ax_d.annotate(txt, xy=(x, ylab[j]), xytext=(3, -10),
                      textcoords='offset points', color=col, fontsize=10,
                      fontweight='bold', ha='left', va='top')

    ax_d.set_xscale('log')
    ax_d.set_xlabel('H₂O volume mixing ratio')
    ax_d.set_ylabel('Altitude (km)')
    ax_d.set_ylim(0., PANEL_ZTOP_KM)
    ax_d.set_facecolor('white')
    # Known residual in the overlay: Kopparapu's tabulated profiles freeze-dry at
    # their cold trap's GRID level, which sits one coarse level (ΔlnP ≈ 0.10-0.19)
    # below the 200 K cap (T* = 202-206 K for Ts ≤ 340 K), inflating their
    # stratospheric H2O by the es_ice(T*)/es_ice(200 K) factor of ~1.3-2.1.
    # ExoColumn freeze-dries at the interpolated 200 K crossing (deliberately —
    # grid-snapped cold traps caused the panel (a)-(c) staircase).  Verified:
    # our ratios x their grid factor = 0.97-0.98 at every profiled Ts.
    if COLD_TRAP_PHASE == 'ice':
        ax_d.text(0.03, 0.03,
                  'Kopparapu profiles freeze-dry at their last grid level\n'
                  'below the 200 K cold trap (202–206 K for $T_s\\leq$ 340 K):\n'
                  'their strat. H₂O reads ×1.3–2 high there',
                  transform=ax_d.transAxes, color=GREY, fontsize=7,
                  ha='left', va='bottom')
    else:
        ax_d.text(0.03, 0.03,
                  'cold trap: sat. over supercooled liquid\n'
                  '(sensitivity variant; CLIMA + model default = ice)',
                  transform=ax_d.transAxes, color=GREY, fontsize=7,
                  ha='left', va='bottom')

    fig.tight_layout(w_pad=2.0, h_pad=2.0)
    out_dir = os.path.dirname(os.path.abspath(__file__))
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(out_dir, f'hz_inner{TAG}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
