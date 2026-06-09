#!/usr/bin/env python3
"""
hz_inner.py  —  ExoColumn version of Kopparapu+2013 Figure 3 (inner HZ, G2V star).

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

Model top: p_top = 0.01 Pa (extended 2 decades below the former 1 Pa) so the H2O
profiles in panel (d) reach ~100-130 km as in Kopparapu Fig 3d.  The shortwave is
converged wrt the top (mass above 0.01 Pa is ~1e-7 of the column; see
tools/diag_ptop_sensitivity.py).  Build at PVER>=200 so the larger log-pressure
span does not coarsen the troposphere (preserves the panel a-c sawtooth behaviour).

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
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

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
# Above this Ts the ExoRT n68 k-tables (T-grid ceiling 500 K, hardcoded in the
# read-only radgrid.F90) under-absorb the near-IR emission of the hot dry base, so
# OLR leaks spuriously off the runaway plateau.  Sweeps that go higher (set
# HZ_TS_MAX) shade this region as opacity-limited (would need HITEMP k-data).
T_KTABLE_RELIABLE = 1600.0   # K

# Water-vapour equation of state for the moist pseudoadiabat:
#   'ideal'    — dilute ideal-gas moist adiabat (malr); the historical default.
#   'nonideal' — Kasting (1988) Appendix-A general non-ideal pseudoadiabat with
#                the native IAPWS-95 EOS (beta, real-gas entropies/cp).  Forces
#                the Wagner-Pruss steam saturation pressure.
# Override with HZ_H2O_EOS=nonideal.  Outputs are tagged so the two modes do not
# overwrite each other (hz_inner.pdf vs hz_inner_nonideal.pdf).
H2O_EOS = os.environ.get('HZ_H2O_EOS', 'ideal')
TAG     = '' if H2O_EOS == 'ideal' else '_' + H2O_EOS

ALBEDO   = 0.32      # surface albedo = Kopparapu et al. (2013) value, for a DIRECT
                     # comparison (their cloud-free 1-D model uses Ad=0.32 to mimic
                     # cloud reflection).  Surface albedo affects only absorbed SW
                     # (hence ASR & Seff), not the LW-driven cold-start T profile/OLR.
                     # [Our ExoRT-tuned Ts=288 K value was 0.24229; kept in git history.]
T_STRATO = 200.0     # isothermal stratosphere cap [K]

MW_H2O   = 18.015    # g/mol

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

# Kopparapu+2013 inner-edge Seff (G2V, 1 M_Earth)
KOPP_RUNAWAY = 1.06     # runaway greenhouse
KOPP_MOIST   = 1.015    # moist greenhouse
SN_LIMIT     = 282.0    # W/m²  Simpson-Nakajima OLR limit
MOIST_GH_VMR = 3.0e-3   # moist greenhouse stratospheric H2O VMR threshold

# Ts values at which to save the full H2O vertical profile for panel (d)
PROFILE_TS = [280.0, 300.0, 320.0, 340.0, 360.0, 380.0]

NML_TEMPLATE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  ihz_profile    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  h2o_eos        = '{h2o_eos}'
  sw_zenith_quad = .true.
  sw_nquad       = 6
/
&exocol_init
  input_file = ''
  ts         = {ts:.2f}
  t_strato   = {t_strato:.1f}
  p_top      = 0.01
  rh_init    = 1.0
  coszrs     = 0.5
  asdir      = {albedo:.4f}
  asdif      = {albedo:.4f}
  aldir      = {albedo:.4f}
  aldif      = {albedo:.4f}
/
&exocol_composition
  n2_vmr   = 0.78
  o2_vmr   = 0.210
  ar_vmr   = 0.01
  co2_vmr  = 3.3e-4
  ch4_vmr  = 0.0
  o3_vmr   = 0.0
/
"""

# ---------------------------------------------------------------------------

def run_one(ts):
    """
    Run ExoColumn flux_only at surface temperature ts.
    Returns dict with scalar diagnostics and profiles, or None on failure.
    """
    nml = NML_TEMPLATE.format(ts=ts, t_strato=T_STRATO, albedo=ALBEDO,
                              h2o_eos=H2O_EOS)
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
            h2ommr = ds['h2ommr'][:]
            mwdry  = float(ds['mw'][:])
            zint   = np.array(ds['zint'][:])   # interface heights [m]

        olr   = float(lwup[0])
        asr   = float(swdn[0] - swup[0])
        alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
        seff  = olr / asr if asr > 1e-6 else np.nan

        # H2O volume mixing ratio profile: VMR ≈ h2ommr × mwdry / MW_H2O
        h2o_vmr = h2ommr * mwdry / MW_H2O
        # Stratospheric H2O VMR = the cold-trap minimum (driest point in the
        # column); used to locate the moist-greenhouse / water-loss limit
        # (the Ts where this reaches MOIST_GH_VMR = 3e-3; Kopparapu Sec 3.1).
        strat_vmr = float(np.nanmin(h2o_vmr)) if np.any(np.isfinite(h2o_vmr)) else np.nan
        # Layer-midpoint altitude [km] from interface heights (TOA=0 … srf=pverp)
        zmid_km = 0.5 * (zint[:-1] + zint[1:]) / 1000.

        return dict(olr=olr, asr=asr, alpha=alpha, seff=seff, strat_vmr=strat_vmr,
                    h2o_vmr=h2o_vmr, pmid=np.array(pmid),
                    tmid=np.array(tmid), zmid_km=zmid_km)
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def main():
    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"Executable not found: {EXE}")

    print("ExoColumn inner-HZ sweep  —  Figure 3 (Kopparapu+2013 style)")
    print(f"  H2O EOS     : {H2O_EOS}" +
          ("  (Kasting 1988 Appendix-A, IAPWS-95)" if H2O_EOS == 'nonideal' else ""))
    print(f"  Ts range    : {TS_VALUES[0]:.0f}–{TS_VALUES[-1]:.0f} K")
    print(f"  t_strato    : {T_STRATO:.0f} K  (isothermal, fully saturated)")
    print(f"  variable_ps : ps = p_N2(1 bar) + esat(Ts)")
    print(f"  composition : N2=0.78 + O2=0.21 + Ar=0.01 + CO2=3.3e-4 + H2O  [no O3/CH4]")
    print()

    rows = []
    profiles = {}   # ts → (pmid, h2o_vmr) for selected Ts
    for ts in TS_VALUES:
        r = run_one(ts)
        if r is None:
            print(f"  Ts={ts:5.1f} K  SKIPPED")
            continue
        print(f"  Ts={ts:5.1f} K  OLR={r['olr']:7.2f}  ASR={r['asr']:7.2f}"
              f"  α={r['alpha']:.3f}  Seff={r['seff']:.4f}")
        rows.append((ts, r['olr'], r['asr'], r['alpha'], r['seff'], r['strat_vmr']))
        if any(abs(ts - ts_p) < 0.5 for ts_p in PROFILE_TS):
            profiles[ts] = (r['zmid_km'], r['h2o_vmr'], r['pmid'])

    if not rows:
        print("No successful runs.")
        return

    arr  = np.array(rows)
    ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a = arr.T

    _plot(ts_a, olr_a, asr_a, alp_a, seff_a, svmr_a, profiles)


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


def _plot(ts, olr, asr, alpha, seff, strat_vmr, profiles):
    # Suppress the vertical-resolution sawtooth (see SMOOTH_WIN note above).
    olr, asr, alpha, seff = (_smooth(olr), _smooth(asr),
                             _smooth(alpha), _smooth(seff))

    # ExoColumn's own moist-greenhouse (water-loss) and runaway-greenhouse Seff:
    #   moist GH = Seff at the first Ts where the stratospheric H2O VMR reaches
    #             MOIST_GH_VMR (3e-3, the water-loss threshold; Kopparapu Sec 3.1).
    #   runaway  = the Seff plateau (limiting value over 600 K..T_KTABLE_RELIABLE,
    #             below the post-runaway rise).
    svmr = _smooth(strat_vmr)
    exo_moist = exo_runaway = moist_ts = np.nan
    _wet = np.where(svmr >= MOIST_GH_VMR)[0]
    if _wet.size:
        exo_moist = float(seff[_wet[0]]); moist_ts = float(ts[_wet[0]])
    _plat = (ts >= 600.) & (ts <= T_KTABLE_RELIABLE)
    if np.any(_plat):
        exo_runaway = float(np.nanmax(seff[_plat]))
    print(f"  ExoColumn moist-GH Seff = {exo_moist:.3f} (Ts={moist_ts:.0f} K)"
          f" ;  runaway Seff = {exo_runaway:.3f}")

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 5.5), dpi=300)
    fig.patch.set_facecolor('white')
    (ax_a, ax_b), (ax_c, ax_d) = axes

    kw = dict(xlim=(200., 2200.))

    # --- (a) OLR and ASR vs Ts ---
    ax_a.plot(ts, olr, color='C3', lw=1.5, label='OLR')
    ax_a.plot(ts, asr, color='C0', lw=1.5, label='Absorbed SW')
    ax_a.axhline(SN_LIMIT, color='k', lw=0.8, ls='--',
                 label=f'S-N limit ({SN_LIMIT:.0f} W m⁻²)')
    if ts[-1] > H2O_TCRIT:
        ax_a.axvline(H2O_TCRIT, color='0.55', lw=0.7, ls='-.',
                     label=f'$T_c$ = {H2O_TCRIT:.0f} K (→ supercrit.)')
    ax_a.set_ylabel('Flux (W m⁻²)')
    ax_a.set_xlabel('$T_s$ (K)')
    ax_a.set_ylim(100., 450.)      # focus on the plateau + onset of the rise
    ax_a.legend(fontsize=7, framealpha=0.9)
    ax_a.set(**kw)
    ax_a.set_facecolor('white')
    ax_a.text(0.04, 0.96, '(a)', transform=ax_a.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (b) Planetary albedo vs Ts ---
    ax_b.plot(ts, alpha, color='C2', lw=1.5)
    ax_b.set_ylabel('Planetary albedo $\\alpha_p$')
    ax_b.set_xlabel('$T_s$ (K)')
    ax_b.set_ylim(0.15, 0.35)
    ax_b.set(**kw)
    ax_b.set_facecolor('white')
    ax_b.text(0.04, 0.96, '(b)', transform=ax_b.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (c) Seff vs Ts ---
    ax_c.plot(ts, seff, color='C1', lw=1.5, label='ExoColumn (G2V)')
    ax_c.axhline(KOPP_RUNAWAY, color='C3', lw=0.9, ls='--',
                 label=f'Kopparapu runaway ({KOPP_RUNAWAY:.4f})')
    ax_c.axhline(KOPP_MOIST, color='C0', lw=0.9, ls=':',
                 label=f'Kopparapu moist GH ({KOPP_MOIST:.4f})')
    if np.isfinite(exo_runaway):
        ax_c.axhline(exo_runaway, color='C1', lw=1.1, ls='--',
                     label=f'ExoColumn runaway ({exo_runaway:.3f})')
    if np.isfinite(exo_moist):
        ax_c.axhline(exo_moist, color='C4', lw=1.1, ls=':',
                     label=f'ExoColumn moist GH ({exo_moist:.3f})')
        ax_c.plot(moist_ts, exo_moist, 'o', color='C4', ms=4, zorder=5)
    ax_c.axhline(1.0, color='gray', lw=0.5, alpha=0.5)
    if ts[-1] > H2O_TCRIT:
        ax_c.axvline(H2O_TCRIT, color='0.55', lw=0.7, ls='-.')
    ax_c.set_ylabel('$S_{\\rm eff}$')
    ax_c.set_xlabel('$T_s$ (K)')
    ax_c.set_ylim(0.4, 1.8)
    ax_c.legend(fontsize=7, framealpha=0.9, loc='lower right')
    ax_c.set(**kw)
    ax_c.set_facecolor('white')
    ax_c.text(0.04, 0.96, '(c)', transform=ax_c.transAxes,
              va='top', fontsize=9, fontweight='bold')

    # --- (d) H2O VMR vertical profiles vs altitude for selected Ts values ---
    items_d = sorted(profiles.items())
    colors_d = plt.cm.plasma(np.linspace(0.05, 0.85, max(len(items_d), 1)))
    for (ts_p, (zmid, vmr, pmid)), col in zip(items_d, colors_d):
        # Trim only the topmost layer(s) where the saturated cold-start hits the
        # es = 0.99*p steam cap (q -> ~0.98); the rest of the stratosphere is the
        # well-behaved uniform cold-trap value.  With p_top = 0.01 Pa this shows
        # the H2O profile up to ~100-130 km (cf. Kopparapu Fig 3d).
        mask = pmid > 0.03
        z = zmid[mask]; v = vmr[mask]
        ax_d.semilogx(v, z, color=col, lw=1.5)
        # label each curve directly at its stratospheric top (highest altitude),
        # where the curves fan out by altitude (warmer ⇒ deeper ⇒ taller)
        itop = int(np.argmax(z))
        ax_d.annotate(f'{ts_p:.0f} K', xy=(v[itop], z[itop]),
                      xytext=(4, 0), textcoords='offset points',
                      color=col, fontsize=6.5, fontweight='bold',
                      ha='left', va='center')
    ax_d.set_xscale('log')
    ax_d.set_xlabel('H₂O VMR')
    ax_d.set_ylabel('Altitude (km)')
    ax_d.set_ylim(bottom=0.)
    ax_d.set_facecolor('white')
    ax_d.text(0.04, 0.96, '(d)', transform=ax_d.transAxes,
              va='top', fontsize=9, fontweight='bold')

    fig.tight_layout()
    out_dir = os.path.dirname(os.path.abspath(__file__))
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(out_dir, f'hz_inner{TAG}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
