#!/usr/bin/env python3
"""
hz_outer.py  —  ExoColumn version of Kopparapu+2013 Figure 5 (outer HZ, G2V star).

REFERENCE CASE: maximum-greenhouse outer HZ (CO2 condensation + T-dependent cp).
==================================================================================
Follows Kopparapu et al. (2013) Section 3.3 exactly:

  - Earth-like planet, 1 bar N2 background, surface temperature FIXED at 273 K.
  - CO2 partial pressure pCO2 swept upward from 1 bar to 34.7 bar = psat_CO2(273 K),
    the same physical endpoint as Kopparapu's "1 to 35 bar".  Layers deeper than
    10 bar use clamped k-table pressure broadening — see the 10-BAR CAP note below.
  - Fully (H2O-)saturated troposphere (rh_init = 1; trivially little water at
    273 K), moist adiabat from the surface.
  - CO2 CONDENSATION: wherever the ascending adiabat would supersaturate in CO2,
    the profile is pinned to the CO2 saturation curve (co2_condense = .true.;
    Kasting 1991) — Span & Wagner (1996) saturation curve, see src/exocol_co2.F90.
  - T-DEPENDENT cp: the CO2 share of cpdry is evaluated at the local temperature
    (cp_co2_tdep = .true.), as in Kopparapu's Shomate-equation update (their
    Sec. 2.1 item 4: cp_CO2 drops ~30% at low T, steepening the dry adiabat).
  - Isothermal stratosphere at t_strato = 154 K (their CO2-condensation cold-trap
    temperature argument).
  - Surface albedo 0.32 (their cloud-free Earth tuning, same as our inner-HZ
    case), 6-point Gauss-Legendre hemispheric zenith average (their 6-angle
    averaging), present solar constant.

For each pCO2, ExoColumn builds the prescribed column (cold start), calls ExoRT
radiation once (flux_only), and we record F_IR (OLR), F_SOL (absorbed SW),
planetary albedo and Seff = F_IR/F_SOL.  The maximum-greenhouse limit is the
Seff minimum.  Three-panel figure matching Kopparapu+2013 Fig. 5 (house style),
with their curves (pixel-digitized from the paper by
tools/digitize_kopparapu_fig5.py -> kopparapu2013_fig5.npz) overlaid for direct
comparison.  Kopparapu's published maximum-greenhouse limit is d = 1.70 AU,
Seff = 0.343 (Fig 5 caption + Table 1; Kasting+1993: 1.67 AU).  The paper is
internally inconsistent at the last digit (Sec. 3.3 text prints "0.325", the
drawn Fig-5(c) curve bottoms at ~0.337 = what we digitize, and the 2014
parametric fit gives 0.356); we overlay the faithful digitized curves and label
the limit with the citable 0.343/1.70 AU headline (see KOPP_* constants below).

>10-BAR LAYERS — CLAMPED PRESSURE BROADENING: the ExoRT n68equiv k-coefficient
pressure grid ends at 10 bar (radgrid.F90, read-only); above it ExoRT linearly
EXTRAPOLATES k in log10(p), which is unphysical (and can produce negative k).
The build now generates src/calc_opd_mod.F90 from ExoRT's source with the
k-table pressure lookup clamped at the 10-bar table edge (same PVER-style
local-copy pattern; see build/Makefile) — i.e. pressure broadening is frozen at
its 10-bar value for deeper layers, while CIA/continuum amagats keep the true
density.  This is defensible because the >10-bar layers are LW-opaque and
barely sunlit, and the albedo rise is Rayleigh-driven (path mass, not k-table).
Columns staying below 10 bar are bit-identical.  Verified: the 8.99-bar sweep
point reproduces the pre-clamp values exactly.

KNOWN RADIATION OFFSET (dense CO2): against Kopparapu's Figure 1 benchmark
(early Mars: 2 bar, 95% CO2 / 5% N2, Ts = 250 K, 167 K isothermal stratosphere)
ExoRT n68equiv gives OLR = 97-100 W/m2 where CLIMA = 86, SMART (line-by-line) =
88.4, Wordsworth+2010 = 88.2 W/m2.  The leak is spectral: ~+9 W/m2 vs SMART in
the 400-850 cm-1 15-um far wings and the 850-1100 cm-1 window (n68equiv has
CO2-CO2 CIA [GBB+dimer] and sub-Lorentzian chi-factors baked into the HITRAN16
k-table, but no CLIMA-style empirical CO2 window continuum; ExoRT carries a
KCO2CONT file only for its n28 band set).  On the SW side n68equiv absorbs less
near-IR sunlight in dense CO2 (planetary albedo at pCO2 = 1 bar: 0.33 vs their
0.28) — same family as the documented inner-HZ near-IR H2O offset.  Both biases
push our Seff HIGH relative to Fig. 5; they are ExoRT spectral-data limitations
(ExoRT is read-only), quantified here by the overlay and in README.md.

Re-sweep takes ~5 min; set HZ_REPLOT=1 to re-plot instantly from the .npz cache.
The full per-run configuration is embedded in NML_TEMPLATE below (binary built
at PVER>=200).
"""

import os
import subprocess
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
HERE     = os.path.dirname(os.path.abspath(__file__))
ROOT     = os.path.abspath(os.path.join(HERE, '..', '..'))
EXE      = os.path.join(ROOT, 'run', 'exocol.exe')
OUT_NC   = os.path.join(ROOT, 'iofiles', 'exocol_out.nc')
NML_PATH = os.path.join(ROOT, 'exocol_config.nml')

# pCO2 grid [bar].  Upper end = psat_CO2(273 K) = 34.7 bar, the same physical
# endpoint as Kopparapu's "1 to 35 bar (the saturation vapor pressure for CO2
# at that temperature)" — beyond it CO2 condenses at the surface itself.
PCO2_MIN = float(os.environ.get('OHZ_PCO2_MIN', '1.0'))
PCO2_MAX = float(os.environ.get('OHZ_PCO2_MAX', '34.7'))
PCO2_N   = int(os.environ.get('OHZ_PCO2_N',   '40'))
PCO2_VALUES = np.geomspace(PCO2_MIN, PCO2_MAX, PCO2_N)

TS       = 273.0     # K, fixed surface temperature (Kopparapu Sec 3.3)
T_STRATO = 154.0     # K, isothermal stratosphere (their CO2 cold-trap argument)
ALBEDO   = 0.32      # surface albedo = Kopparapu et al. (2013) Earth tuning
P_N2     = 1.0e5     # Pa, fixed N2 background

# pCO2 values [bar] for which the T(p) profile is cached (diagnostics; shows the
# CO2-saturation-pinned upper troposphere).
PROFILE_PCO2 = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]

CACHE = os.path.join(HERE, 'hz_outer.npz')
KOPP  = os.path.join(HERE, 'kopparapu2013_fig5.npz')

# Kopparapu et al. (2013) maximum-greenhouse limit for the Sun, AS PUBLISHED.
# We quote the universally-cited headline value from the Figure 5 caption and
# Table 1: d = 1.70 AU, Seff = 0.343 (self-consistent: 1/sqrt(0.343) = 1.707 AU).
# The paper is internally inconsistent at the decimal level: the Sec. 3.3 TEXT
# prints "Seff = 0.325" (= 1.75 AU, and below the drawn Fig-5(c) curve), the
# drawn Fig-5(c) curve bottoms at ~0.337 (what our pixel digitization recovers,
# verified against the rendered figure), and the 2014 parametric fit (ApJL 787,
# L29, Eq. 4) lists Seff_sun = 0.356 (= 1.68 AU).  The digitized curves overlaid
# below are faithful to Fig 5; this constant is only used to LABEL the limit with
# the citable caption/Table-1 number rather than our ~0.337 pixel read.
KOPP_SEFF_MAX = 0.343
KOPP_PCO2_MAX = 8.0     # bar (their "~8 bar")
KOPP_D_MAX    = 1.70    # AU, as published (1/sqrt(0.343) = 1.707, rounded to 1.70)

NML_TEMPLATE = """\
&exocol_nml
  flux_only      = .true.
  variable_ps    = .true.
  co2_condense   = .true.
  cp_co2_tdep    = .true.
  o3_profile     = 'none'
  msdist         = 1.0
  solar_file     = '{solar_file}'
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
  ps       = {p_dry:.6e}
  co2_vmr  = {co2_vmr:.8f}
  n2_vmr   = {n2_vmr:.8f}
  o2_vmr   = 0.0
  ch4_vmr  = 0.0
  o3_vmr   = 0.0
  ar_vmr   = 0.0
/
"""

# ---------------------------------------------------------------------------

def run_one(pco2_bar, solar_file=''):
    """Run ExoColumn flux_only at CO2 partial pressure pco2_bar [bar].
    solar_file selects the host-star n68 SED ('' => compile-time G2V_SUN_n68.nc,
    the Sun); used by the multi-stellar Figure-6 sweep (hz_figure6.py).
    Returns dict of scalar fluxes + the T(p) profile, or None on failure."""
    p_dry   = P_N2 + pco2_bar * 1.0e5           # dry surface pressure [Pa]
    co2_vmr = pco2_bar * 1.0e5 / p_dry          # dry-air mole fractions
    n2_vmr  = P_N2 / p_dry
    nml = NML_TEMPLATE.format(ts=TS, t_strato=T_STRATO, albedo=ALBEDO,
                              p_dry=p_dry, co2_vmr=co2_vmr, n2_vmr=n2_vmr,
                              solar_file=solar_file)
    orig = None
    if os.path.exists(NML_PATH):
        with open(NML_PATH) as f:
            orig = f.read()
    try:
        with open(NML_PATH, 'w') as f:
            f.write(nml)
        result = subprocess.run([EXE], cwd=ROOT, capture_output=True,
                                text=True, timeout=300)
        if result.returncode != 0:
            print(f"  FAIL pCO2={pco2_bar:.2f} bar  rc={result.returncode}")
            print(result.stderr[-300:] if result.stderr else '')
            return None
        with nc.Dataset(OUT_NC) as ds:
            lwup = ds['LWUP'][:]
            swup = ds['SWUP'][:]
            swdn = ds['SWDN'][:]
            tmid = np.array(ds['tmid'][:])
            pmid = np.array(ds['pmid'][:])
        olr   = float(lwup[0])
        asr   = float(swdn[0] - swup[0])
        alpha = float(swup[0] / swdn[0]) if swdn[0] > 0 else 0.0
        seff  = olr / asr if asr > 1e-6 else np.nan
        return dict(olr=olr, asr=asr, alpha=alpha, seff=seff,
                    tmid=tmid, pmid=pmid)
    finally:
        if orig is not None:
            with open(NML_PATH, 'w') as f:
                f.write(orig)
        elif os.path.exists(NML_PATH):
            os.remove(NML_PATH)


def _save_cache(arr, profiles):
    ks = sorted(profiles)
    np.savez(CACHE, rows=arr,
             prof_pco2=np.array(ks, dtype=float),
             prof_t=np.array([profiles[k][0] for k in ks]),
             prof_p=np.array([profiles[k][1] for k in ks]))
    print(f"Cached: {CACHE}")


def _load_cache():
    d = np.load(CACHE)
    profiles = {float(c): (t, p) for c, t, p in
                zip(d['prof_pco2'], d['prof_t'], d['prof_p'])}
    return d['rows'], profiles


def main():
    if os.environ.get('HZ_REPLOT'):
        if not os.path.exists(CACHE):
            raise FileNotFoundError(f"No cache to re-plot: {CACHE} (run the sweep first)")
        print(f"HZ_REPLOT: re-plotting from {CACHE}")
        arr, _profiles = _load_cache()
        _plot(*arr.T)
        return

    if not os.path.isfile(EXE):
        raise FileNotFoundError(f"Executable not found: {EXE}")

    print("ExoColumn outer-HZ sweep  —  Figure 5 (Kopparapu+2013 style)")
    print(f"  pCO2 range  : {PCO2_VALUES[0]:.2f}–{PCO2_VALUES[-1]:.2f} bar "
          f"({PCO2_N} points; >10 bar layers use clamped k-table broadening)")
    print(f"  Ts          : {TS:.0f} K fixed;  t_strato = {T_STRATO:.0f} K")
    print(f"  composition : 1 bar N2 + pCO2;  rh_init = 1;  albedo = {ALBEDO}")
    print(f"  physics     : co2_condense + cp_co2_tdep;  6-pt zenith quadrature")
    print()

    rows = []
    profiles = {}
    for pco2 in PCO2_VALUES:
        r = run_one(pco2)
        if r is None:
            print(f"  pCO2={pco2:6.2f} bar  SKIPPED")
            continue
        print(f"  pCO2={pco2:6.2f} bar  OLR={r['olr']:7.2f}  ASR={r['asr']:7.2f}"
              f"  α={r['alpha']:.3f}  Seff={r['seff']:.4f}")
        rows.append((pco2, r['olr'], r['asr'], r['alpha'], r['seff']))
        if any(abs(pco2 - c) / c < 0.04 for c in PROFILE_PCO2):
            profiles[pco2] = (r['tmid'], r['pmid'])

    if not rows:
        print("No successful runs.")
        return
    arr = np.array(rows)
    _save_cache(arr, profiles)
    _plot(*arr.T)


def _plot(pco2, olr, asr, alpha, seff):
    # Maximum-greenhouse limit = the Seff minimum (quadratic fit through the
    # three points around the discrete minimum for a smooth location).
    imin = int(np.argmin(seff))
    if 0 < imin < len(seff) - 1:
        x = np.log10(pco2[imin-1:imin+2]); y = seff[imin-1:imin+2]
        c = np.polyfit(x, y, 2)
        x_min = -c[1] / (2 * c[0])
        p_min = 10 ** x_min
        s_min = float(np.polyval(c, x_min))
        at_edge = False
    else:
        p_min, s_min = float(pco2[imin]), float(seff[imin])
        at_edge = True
    d_min = 1.0 / np.sqrt(s_min)
    print(f"\n  ExoColumn maximum greenhouse: Seff = {s_min:.4f} at "
          f"pCO2 = {p_min:.2f} bar  ->  d = {d_min:.3f} AU"
          + ("  [AT SWEEP EDGE — minimum not resolved]" if at_edge else ""))
    print(f"  Kopparapu+2013              : Seff = {KOPP_SEFF_MAX:.3f} at "
          f"pCO2 ~ {KOPP_PCO2_MAX:.0f} bar  ->  d = {KOPP_D_MAX:.2f} AU  "
          f"(Fig 5 caption / Table 1)")

    kp = np.load(KOPP) if os.path.exists(KOPP) else None

    fig, (ax_a, ax_b, ax_c) = plt.subplots(3, 1, figsize=(5.0, 9.0), dpi=300,
                                           sharex=True)
    fig.patch.set_facecolor('white')
    kw_kopp = dict(color='0.55', lw=1.1, ls='--', zorder=2)

    # --- (a) F_IR and F_SOL vs pCO2 ---
    if kp is not None:
        ax_a.plot(kp['fir_p'],  kp['fir'],  **kw_kopp)
        ax_a.plot(kp['fsol_p'], kp['fsol'], **kw_kopp)
    ax_a.plot(pco2, olr, color='C3', lw=1.6, zorder=4)
    ax_a.plot(pco2, asr, color='C0', lw=1.6, zorder=4)
    ax_a.set_ylabel('Flux (W m⁻²)')
    ax_a.set_ylim(50., 260.)
    ax_a.text(2.6, 238., '$F_\\mathrm{SOL}$', color='C0', fontsize=11,
              fontweight='bold', ha='center', va='bottom')
    ax_a.text(2.6, 95., '$F_\\mathrm{IR}$', color='C3', fontsize=11,
              fontweight='bold', ha='center', va='top')
    ax_a.text(16., 152., 'Kopparapu+2013', color='0.45', fontsize=8,
              ha='center', va='bottom', rotation=-18)
    ax_a.set_facecolor('white')

    # --- (b) planetary albedo vs pCO2 ---
    if kp is not None:
        ax_b.plot(kp['alb_p'], kp['alb'], **kw_kopp)
    ax_b.plot(pco2, alpha, color='C2', lw=1.6, zorder=4)
    ax_b.axhline(ALBEDO, color='0.3', lw=0.9, ls=':')
    ax_b.set_ylabel('Planetary albedo $\\alpha_p$')
    ax_b.set_ylim(0.2, 0.7)
    ax_b.text(13., ALBEDO + 0.008, f'Surface albedo = {ALBEDO:.2f}',
              color='0.3', fontsize=8, ha='left', va='bottom')
    ax_b.text(16., 0.475, 'Kopparapu+2013', color='0.45', fontsize=8,
              ha='center', va='top', rotation=14)
    ax_b.set_facecolor('white')

    # --- (c) Seff vs pCO2 ---
    if kp is not None:
        ax_c.plot(kp['seff_p'], kp['seff'], **kw_kopp)
        ks = kp['seff']; kpp = kp['seff_p']
        ki = int(np.argmin(ks))
        # Mark the drawn-curve minimum (the dot sits on the digitized curve),
        # but LABEL it with Kopparapu's published headline value (Fig 5 caption /
        # Table 1: Seff = 0.343, d = 1.70 AU) rather than our ~0.337 pixel read.
        ax_c.plot(kpp[ki], ks[ki], 'o', color='0.55', ms=3, zorder=3)
        ax_c.annotate(f'Kopparapu+2013\n$S_\\mathrm{{eff}}$ = '
                      f'{KOPP_SEFF_MAX:.3f} ({KOPP_D_MAX:.2f} AU)',
                      xy=(kpp[ki], ks[ki]), xytext=(40, -4),
                      textcoords='offset points', color='0.45', fontsize=8,
                      ha='left', va='top')
    ax_c.plot(pco2, seff, color='C1', lw=1.6, zorder=4)
    if at_edge:
        # Seff is still falling at the sweep edge: quote the edge value as a
        # bound on the maximum-greenhouse limit, not a minimum.
        ax_c.text(10.3, 0.435, 'Maximum greenhouse\n'
                  f'$S_\\mathrm{{eff}} \\leq$ {s_min:.3f} '
                  f'($d \\geq$ {d_min:.2f} AU)\nat the sweep edge',
                  color='0.25', fontsize=8, ha='left', va='bottom')
    else:
        ax_c.plot(p_min, s_min, 'o', color='0.3', ms=3.5, zorder=6)
        ax_c.annotate(f'Maximum greenhouse\n$S_\\mathrm{{eff}}$ = {s_min:.3f}'
                      f' ({d_min:.2f} AU)',
                      xy=(p_min, s_min), xytext=(10, 12),
                      textcoords='offset points', color='0.25', fontsize=8,
                      ha='center', va='bottom')
    ax_c.set_ylabel('$S_\\mathrm{eff}$')
    ax_c.set_xlabel('CO₂ partial pressure (bar)')
    ax_c.set_xscale('log')
    ax_c.set_xlim(1., 35.)
    ax_c.set_ylim(0.30, 0.60)
    ax_c.set_xticks([1, 2, 5, 10, 20, 35])
    ax_c.set_xticklabels(['1', '2', '5', '10', '20', '35'])
    ax_c.set_facecolor('white')

    fig.tight_layout()
    for ext, dpi in [('pdf', 300), ('png', 150)]:
        path = os.path.join(HERE, f'hz_outer.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        print(f"Saved: {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
