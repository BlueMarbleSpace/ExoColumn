# Outer HZ (maximum greenhouse) reference case

ExoColumn version of Kopparapu et al. (2013) Section 3.3 / Figure 5: the
outer-edge (maximum greenhouse) habitable-zone limit for a G2V star.

## Method (Kopparapu §3.3, inverse climate calculation — no RCE)

Earth-like planet, **1 bar N2**, surface fixed at **Ts = 273 K**, CO2 partial
pressure swept **1 → 34.7 bar** (= psat_CO2(273 K), the same physical endpoint
as Kopparapu's "1 to 35 bar, the saturation vapor pressure at that
temperature"; layers deeper than 10 bar use clamped k-table pressure
broadening — see below).  For each pCO2 the cold start builds the prescribed
column — H2O-saturated moist adiabat from the surface, **pinned to the CO2
saturation curve wherever the ascent supersaturates in CO2** (Kasting 1991;
`co2_condense`), with the **CO2 share of cp evaluated at the local temperature**
(`cp_co2_tdep`, Kopparapu's Shomate update), capped by an **isothermal 154 K
stratosphere** — then calls ExoRT once (`flux_only`).  Surface albedo 0.32,
6-point Gauss–Legendre hemispheric zenith average, present solar constant.
Seff = F_IR/F_SOL; the maximum-greenhouse limit is the Seff minimum;
d = 1/√Seff.

New physics added for this case (both namelist-gated, default **off**; Earth RCE
and the inner-HZ reference verified bit-identical):

- `src/exocol_co2.F90` — CO2 saturation curve (Span & Wagner 1996 auxiliary
  equations, sublimation + liquid branches; the modern replacement for Kasting's
  Fanale et al. 1982 fits), frost-point inversion, and statistical-mechanics
  cp_CO2(T) (7/2 R + Einstein vibrational terms; matches NIST-JANAF to <0.5%
  over 150–650 K).  Validated by `test/test_co2.F90` (triple point, 1-atm
  dry-ice point, critical point, JANAF cp anchors — all pass).
- `&exocol_nml` switches `co2_condense` and `cp_co2_tdep` (see CLAUDE.md /
  `exocol_config.F90`), applied in the standard cold-start profile builder.

## Files

- `hz_outer.py` — sweep + figure generator.  Full sweep ~5 min
  (`python3 reference/max_greenhouse/hz_outer.py`, binary at PVER≥200, run with
  the Intel runtime sourced); `HZ_REPLOT=1` re-plots instantly from the cache.
- `hz_outer.npz` — cached sweep results (also caches T(p) profiles at
  pCO2 = 1, 2, 4, 8 bar for diagnostics).
- `hz_outer.pdf` / `.png` — the three-panel figure (Fig. 5 format, house style),
  with the digitized Kopparapu curves overlaid.
- `kopparapu2013_fig5.npz` — Kopparapu Fig. 5 curves, pixel-digitized from the
  paper PDF by `tools/digitize_kopparapu_fig5.py` (validated: F_IR(1 bar) =
  111.6 ≈ 112 W/m², albedo(1 bar) = 0.278, Seff min = 0.337 at 6.7 bar ≈ the
  caption values).

## Results (2026-06-10 sweep, PVER=200, full 1–34.7 bar range)

> **Cold-trap phase note (2026-06-11).**  This case briefly switched to
> `cold_trap_phase='liquid'` (Seff 0.395 → 0.384, d 1.59 → 1.61 AU) on the
> belief that supercooled-liquid saturation below 273.16 K was CLIMA's
> convention.  Inspection of the actual CLIMA source (atmos repo: `SATRAT`
> uses the sublimation latent heat below the triple point; `convec.f` label 13
> is the ice-saturated sub-freezing pseudoadiabat) showed CLIMA is fully
> **ice-based** below 273.16 K — i.e. ExoColumn's `'ice'` default *is* the
> Kopparapu-consistent choice — so the liquid run was reverted and the
> ice-based results below stand as the reference.

| quantity | ExoColumn | Kopparapu+2013 |
|---|---|---|
| F_IR at pCO2 = 1 bar | 123.1 W/m² | 111.6 W/m² |
| F_IR asymptote (≥ 15 bar) | 72.3–72.5 W/m² | ~65 W/m² |
| albedo at pCO2 = 1 bar | 0.328 | 0.278 |
| albedo at pCO2 = 34.7 bar | 0.538 | ~0.54 |
| Seff at pCO2 = 1 bar | 0.539 | 0.455 |
| **maximum greenhouse** | **Seff = 0.395 at 8.9 bar → d = 1.59 AU** | Seff = 0.337 at ~7 bar → d = 1.70 AU |

The Seff minimum is now resolved interior to the sweep.  The structure mirrors
Kopparapu's panel by panel: F_IR falls to a flat asymptote once the atmosphere
is LW-opaque (ours from ~15 bar, theirs from ~10), the albedo climbs with
Rayleigh scattering, and their competition produces the Seff minimum.
Notably, **F_SOL and the planetary albedo converge onto Kopparapu's curves at
high pCO2** (the dense-CO2 SW budget becomes Rayleigh-dominated, which both
models compute from the same Vardavas & Carver data — the near-IR absorption
difference stops mattering), while F_IR keeps a parallel ~+7 W/m² offset (the
LW window/far-wing gap below).  Both offsets push Seff high, so our 1.59 AU is
conservative relative to their 1.70 AU.

## Diagnosed radiation offsets (ExoRT n68equiv vs CLIMA, dense CO2)

**LW (~+11 W/m² at 1 bar):** reproduced Kopparapu's Figure 1 benchmark (early
Mars: 2 bar, 95% CO2/5% N2, Ts = 250 K, 167 K isothermal stratosphere):
ExoRT gives OLR = 97.4 (pure adiabat) / 100.0 (with CO2-condensation pinning)
W/m² where CLIMA = 86, SMART (line-by-line) = 88.4, and Wordsworth+2010 = 88.2.
Band-resolved comparison against the digitized Fig. 1 spectra localizes the leak
to the 400–850 cm⁻¹ 15-µm far wings (+9 W/m² vs SMART) and the 850–1100 cm⁻¹
window (+6 vs SMART).  n68equiv *does* carry CO2–CO2 CIA (GBB + dimer, LW & SW)
and sub-Lorentzian χ-factors (HITRAN2016 k-table, 500 cm⁻¹ line cut — same
conventions as Kopparapu/Wordsworth), but has no CLIMA-style empirical CO2
window continuum (ExoRT ships a `KCO2CONT` file only for its n28 band set).

**SW (albedo +0.05 at 1 bar, narrowing to +0.02 at 9 bar):** n68equiv absorbs
less near-IR sunlight in dense CO2 than CLIMA — the same family as the
documented inner-HZ near-IR H2O offset.  Rayleigh is *not* the cause (both
models use Vardavas & Carver 1984 CO2/N2 coefficients).

Both biases push Seff high (planet must sit closer in), so our d ≥ 1.59 AU is
*conservative* relative to Kopparapu's 1.70 AU.  These are ExoRT spectral-data
limitations (ExoRT is read-only), not methodology differences.

## Layers deeper than 10 bar: clamped pressure broadening

The n68equiv k-coefficient pressure grid ends at 10 bar (`radgrid.F90`).  Above
it ExoRT linearly **extrapolates** k in log10(p) (`rad_interp_mod`: "at the
tops of grids, extrapolate") — at 35 bar the extrapolation factor is ~6×, which
is unphysical and can even produce negative k.  The build therefore generates
`src/calc_opd_mod.F90` from ExoRT's source with one line patched (the same
PVER-style local-copy pattern; see `build/Makefile`, which fails loudly if the
anchor line ever changes upstream):

```fortran
pressure = min(log10(pmid(ik)), log10pgrid(kc_npress))
```

i.e. the k-table lookup is clamped at the 10-bar table edge — pressure
broadening frozen at its 10-bar value for deeper layers — while the
CIA/continuum paths keep the true `pmid` partial pressures and amagats.  This
is defensible because the >10-bar layers are LW-opaque and barely sunlit, and
the high-pCO2 albedo rise is Rayleigh-driven (path mass, not k-tables; confirmed
by our albedo converging onto Kopparapu's at 35 bar).  Columns that stay below
10 bar are **bit-identical** (verified: the Earth/IHZ flux anchor and the
8.99-bar sweep point reproduce exactly after the clamp build).
