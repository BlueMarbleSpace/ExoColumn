# Outer HZ (maximum greenhouse) reference case

> **ALBEDO UPDATE (2026-08-21).** Surface albedo changed from 0.32 to ExoColumn's
> own Earth calibration **α_s = 0.2736** (see `reference/albedo_sensitivity/`);
> the previous set is archived as `hz_outer_a032.{npz,pdf,png}`.
> **Primary maximum greenhouse: Seff = 0.385 at pCO₂ = 8.68 bar → d = 1.611 AU**
> (was 0.395 at 8.92 bar → 1.591 AU). The outer edge is nearly albedo-insensitive
> — dSeff/dα_s = 0.21, versus 0.65 at the inner edge — because the dense-CO₂
> Rayleigh layer screens the surface; the two α_p curves converge as pCO₂ rises.
> Comparison numbers below that are not restated here refer to the archived
> α_s = 0.32 set.


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
stratosphere** — then calls ExoRT once (`flux_only`).  Surface albedo 0.2736,
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
  paper PDF by `tools/digitize_kopparapu_fig5.py`.  The digitization is faithful
  to the figure (verified panel-by-panel against the rendered Fig 5): F_IR(1 bar)
  = 111.6 ≈ 112 W/m², F_IR asymptote ≈ 65, F_SOL(1 bar) ≈ 245, albedo 0.278 → 0.54,
  and the Seff curve bottoms at ~0.337 just below the 0.34 gridline at ~7–8 bar —
  exactly where the drawn blue curve sits.  **Published headline value (what we
  label in the figure): d = 1.70 AU, Seff = 0.343 (Fig 5 caption + Table 1).**
  Note Kopparapu's paper is internally inconsistent at the last digit: the §3.3
  text misprints "Seff = 0.325" (= 1.75 AU, below the drawn curve), the drawn
  curve minimum reads ~0.337 (= 1.72 AU), and the 2014 parametric fit (ApJL 787,
  L29, Eq. 4) lists Seff_sun = 0.356 (= 1.68 AU).  The 1.70 AU / 0.343 caption
  value is self-consistent (1/√0.343 = 1.707) and the one universally cited; the
  2014 mass-dependence paper confirms the maximum-greenhouse coefficients are
  unchanged from 2013.

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
| **maximum greenhouse** | **Seff = 0.395 at 8.9 bar → d = 1.59 AU** | Seff = 0.343 at ~8 bar → d = 1.70 AU |

(Kopparapu's maximum-greenhouse entry is his published caption/Table-1 headline,
Seff = 0.343 / 1.70 AU; the digitized curve we overlay bottoms at ~0.337 — see
the `kopparapu2013_fig5.npz` note above for the paper's internal 0.325/0.337/0.356
spread.)

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

**LW:**

> **GRAVITY CORRECTION (2026-08-24).**  The Figure-1 reproduction described in
> this paragraph was run at **Earth** gravity.  Kopparapu's Figure-1 case is a
> **Mars-mass** planet (their `input_clima.dat` uses `G = 373` cm/s²; a
> hydrostatic fit to the P–T–z profile they handed to SMART independently gives
> g₀ = 3.70 m/s² for a 95 % CO2 / 5 % N2 column), so the correct column mass is
> 2.6× larger than was used.  Redone properly — Mars gravity, and Kopparapu's own
> SMART profile rather than an ExoColumn cold start — ExoRT n68 gives
> **OLR = 82.3 W/m²**, versus 94.5 W/m² for the identical column at Earth gravity
> (which is what produced the 97.4 quoted below).  At Mars gravity ExoRT sits
> **between** Clima and SMART, not above both: Clima 2013-era k = 77.0,
> Clima Wolf-2016 k = 80.9, ExoRT n68 = 82.3, SMART = 88.5 W/m².  The
> "+11 W/m² dense-CO2 LW leak" inferred from the numbers below is therefore
> **not supported**; see "Kopparapu Figure-1 benchmark" below for the corrected
> case, which is now a manuscript figure panel.  The independent LBL benchmark
> on our *own* max-greenhouse column (next section) is unaffected — it was always
> run at Earth gravity, which is correct for that Earth-mass case.

*(superseded, retained for the record)* reproduced Kopparapu's Figure 1 benchmark (early
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

### Independent line-by-line benchmark (this column) — `tools/lbl_co2_benchmark.py`

The Fig.-1 comparison above uses *literature* LBL values (SMART, Wordsworth) on
the 2-bar early-Mars column.  `tools/lbl_co2_benchmark.py` adds a from-scratch
line-by-line OLR for **our own maximum-greenhouse limit column** (pCO2 = 8.87 bar
+ 1 bar N2, Ts = 273 K — the Seff-minimum), the dense-CO2 twin of the inner-HZ
`tools/lbl_olr_benchmark.py`.  Physics RADIS does not supply, added here:
Perrin & Hartmann (1989) sub-Lorentzian χ-factor on the CO2 line wings (500 cm⁻¹
cut, ClearSky.jl/Wordsworth coefficients), HITRAN-2024 CO2-CO2 CIA (carries the
Gruszka-Borysow far-IR + Baranov 7-µm bands), trace H2O lines + MT_CKD; the CO2
line sum reproduces RADIS Voigt to ~1 % at χ = 1 (`tools/check_co2_lbl.py`).

**Result (10–2000 cm⁻¹):**

| treatment | OLR (10–2000 cm⁻¹) | 8–12 µm window |
|---|---|---|
| LBL pure-Lorentz wings (opaque bound) | 44.1 W/m² | 3.5 W/m² |
| Clima 2013-era (Kopparapu) | 69.4 W/m² | 8.2 W/m² |
| Clima Wolf-HITRAN2016 | 72.3 W/m² | 8.0 W/m² |
| **ExoRT n68 (this work)** | **75.5 W/m²** | **9.0 W/m²** |
| LBL PH89 χ sub-Lorentzian (transparent bound) | 78.5 W/m² | 14.7 W/m² |

**ExoRT n68 sits inside the wing-treatment envelope, 3.0 W/m² (3.9 %) below the
PH89 LBL and clustered with Clima.**

> **Weak-line pruning fix (2026-08-21).**  The PH89 LBL previously read
> **91.7** W/m² (window 17.6), and the apparent ~16 W/m² window disagreement with
> every band/k model was an artefact of `co2_line_tau`'s own line-selection
> cutoff, not physics.  It pruned lines weaker than `1e-6 × S.max()`; on this
> column (N_col ≈ 1.3e26 molec/cm²) a line *at* that threshold still had a peak
> optical depth of ~215, and in the 8–12 µm window — which has no strong band to
> carry the opacity — **17 540 of 17 610 lines were being discarded**.  The cutoff
> is now on peak *column optical depth* (`tau_min = 1e-3`, converged: the total
> moves 0.017 W/m² between `tau_min` 1e-2 and 1e-5), which is also immune to the
> old criterion's silent dependence on the requested spectral range.  Effect is
> 96 % confined to the window: total 91.7 → 78.5, window 28.2 → 15.6, while
> 10–400 cm⁻¹ is unchanged to 0.01 W/m² and the pure-Lorentz bound barely moves
> (44.3 → 44.1 — with unclipped Lorentz wings the strong lines already saturate
> the window, so the weak lines add nothing there).  `tools/check_co2_lbl.py`
> check 4 now tests pruning convergence *in the window*; the pre-existing checks
> ran only at 580–720 cm⁻¹, inside the 15 µm band, where any sane cutoff is
> harmless — which is why this went unnoticed.

After the fix the band/k models and the LBL agree to |Δ| ≤ 0.25 W/m² per n68 band
across 720–1108 cm⁻¹, and **essentially the whole remaining residual is one band,
1108–1200 cm⁻¹ (LBL 5.34 vs ExoRT 2.14 W/m²)**.  That band straddles the edge of
the HITRAN-2024 CO2–CO2 CIA tabulation, which covers 1–750 and 1150–1850 cm⁻¹ with
**nothing in between** (directly verified; ExoRT's own GB/Baranov-derived CIA is
likewise ~0 across 875–1150).  So the residual is plausibly missing continuum/CIA
opacity on our side of the 1150 cm⁻¹ edge rather than a k-distribution error — but
that is a hypothesis, not a demonstration.  The three band/k models (ExoRT n68,
Clima 2013-era, Clima Wolf2016) still agree to within 6 W/m² overall and have
comparably opaque windows (8–9 W/m²), sharing the HITRAN-lineage CO2
k-distribution.  The 2013↔Wolf2016 shift is only ~3 W/m² here (CO2 k-coeffs), vs
~16 W/m² in the H2O-dominated IHZ.

The figure (`lbl_olr_benchmark_ohz.{png,pdf}`, `tools/plot_lbl_olr_figs.py`)
puts the max-greenhouse panel (left) alongside the early-Mars panel (right) in
the same style as the inner-edge twin in `reference/moist_runaway/`
(`lbl_olr_benchmark_ihz.{png,pdf}`, rendered by the same script): grey
line-by-line LBL, black LBL n68-band averages, red ExoRT n68, and
Clima in **both** k-coefficient generations — green for the 2013-era set used by
Kopparapu et al. (2013), blue for the post-2014 Wolf HITRAN-2016 set adopted by
the atmos repo in 2021 — with a `model − LBL` residual beneath (bars = ExoRT,
step curves = the two Clima generations).  The post-2014 curve is what shows that
the ExoRT–Clima gap is a k-coefficient-generation effect and not a model
difference: on the IHZ column Clima moves 251.8 → 267.5 W/m² against ExoRT's
269.9 and the LBL's 272.6, and on the OHZ column 69.4 → 72.3 against ExoRT's 75.5.  The reference LBL on the
right is the **PH89-χ** case (the sub-Lorentzian convention Kopparapu/CLIMA use);
the pure-Lorentz bound is the table's opaque end, not drawn.  After the
weak-line pruning fix the models track the LBL closely across the window; the
visible exception is the 1108–1200 cm⁻¹ band discussed above.
`lbl_olr_co2_maxgh.npz` caches both LBL wing bounds,
`clima_band_olr_maxgh.txt` the two Clima generations (regenerated by
running `/models/atmos` Clima inverse at TG0=273/PGO=9.87/fCO2≈0.9 through the
band-OLR `ir.f` patch, last 55 rows ×10⁻³; swap `ClimaMain.f` units 15–18 for the
k-coeff generation), `exocol_maxgh_8.87bar.nc` the benchmarked column.

### Kopparapu Figure-1 benchmark (early Mars) — the second OHZ panel

The right-hand panel of `lbl_olr_benchmark_ohz` is the **Kopparapu et al. (2013)
Figure-1 configuration** run through the same four radiation models: a Mars-mass
planet (**g = 3.73 m/s²** = their `G = 373` cm/s²), **2 bar of 95 % CO2 / 5 % N2**,
Ts = 250 K, 167 K isothermal stratosphere, saturated H2O.  Unlike the
max-greenhouse column on the left — which is ExoColumn's own — every model here
integrates **Kopparapu's own atmosphere**: the P–T profile is taken verbatim from
the file they handed to SMART (`smart_earlymars_thermal_newcia.hrt`, surface
2.001 bar / 250 K, top 5.65e-5 bar / 167 K), the H2O profile from the matched
Clima run, and Clima's own inverse solution on this deck reproduces that P–T to
< 0.25 K through the stratosphere and upper troposphere (5.4 K max, at the
CO2-condensation kink near 1 bar).

| model (10–2000 cm⁻¹) | OLR |
|---|---|
| LBL pure-Lorentz wings (opaque bound) | 46.3 W/m² |
| Clima 2013-era k (this reproduction) | 77.0 W/m² |
| Clima Wolf-HITRAN2016 k | 80.9 W/m² |
| **ExoRT n68 (this work)** | **82.3 W/m²** |
| LBL PH89 χ sub-Lorentzian (reference) | 82.9 W/m² |
| **SMART** (Kopparapu, updated CIA) | **88.5 W/m²** |
| published in Kopparapu et al. (2013) Fig. 1 | 86 (Clima), 88.4 (SMART) |

**ExoRT n68 reproduces the line-by-line reference to 0.6 % on this column** —
the closest agreement of any panel in either benchmark figure — and the
band-resolved residuals are ≤ 0.62 W/m² everywhere.  SMART is 5.6 W/m² above our
LBL, and essentially all of that sits in the 220–500 cm⁻¹ CO2 far wings, the same
region Kopparapu identified as the source of their Clima–SMART difference.  Our
reproduction of their 2013-era Clima lands 9 W/m² below their published 86; the
identical ~9 W/m² shortfall appears in the dense-H2O Figure-2 reproduction
(`reference/moist_runaway`), so it is a property of the current public `atmos`
Clima relative to the 2013 version rather than of this column.

Files: `lbl_olr_co2_earlymars.npz` (LBL, both wing bounds; regenerate with
`tools/lbl_co2_benchmark.py --g 3.73`), `clima_band_olr_earlymars.txt` (both
Clima generations), `smart_earlymars_olr.txt` (Kopparapu's SMART spectrum,
converted to W/m²/cm⁻¹), `exocol_earlymars_2bar.nc` (the benchmarked column;
requires a `make PVER=200 EXO_G=3.73` binary — gravity is compile-time).

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
