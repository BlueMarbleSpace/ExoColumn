# Outer HZ (maximum greenhouse) reference case

ExoColumn version of Kopparapu et al. (2013) Section 3.3 / Figure 5: the
outer-edge (maximum greenhouse) habitable-zone limit for a G2V star.

## Method (Kopparapu §3.3, inverse climate calculation — no RCE)

Earth-like planet, **1 bar N2**, surface fixed at **Ts = 273 K**, CO2 partial
pressure swept upward from 1 bar (Kopparapu: 1–35 bar; here **1–8.99 bar**, see
the 10-bar cap below).  For each pCO2 the cold start builds the prescribed
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

## Results (2026-06-10 sweep, PVER=200)

| quantity | ExoColumn | Kopparapu+2013 |
|---|---|---|
| F_IR at pCO2 = 1 bar | 123.1 W/m² | 111.6 W/m² |
| albedo at pCO2 = 1 bar | 0.328 | 0.278 |
| Seff at pCO2 = 1 bar | 0.539 | 0.455 |
| Seff at 8.99 bar (cap) | 0.3949, still falling | 0.341 (their curve) |
| maximum greenhouse | **Seff ≤ 0.395 → d ≥ 1.59 AU** (min beyond cap) | Seff = 0.337 → d = 1.70 AU |

The curve *shapes* track Kopparapu well (F_IR falling toward an asymptote,
albedo rising with Rayleigh scattering, Seff minimum from their competition);
the offsets are systematic and diagnosed below.  Our Seff is still decreasing
(by <0.001 per step) at the 8.99-bar sweep edge, so our maximum-greenhouse
minimum lies just beyond the current 10-bar cap.

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

## The 10-bar k-table cap

The n68equiv k-coefficient pressure grid ends at 10 bar (`radgrid.F90`); above
it ExoRT's reference-pressure search runs off the table (out-of-bounds index in
`rad_interp_mod`).  The sweep therefore stops at total ps = 10 bar
(pCO2 = 8.99 bar).  Extending to Kopparapu's full 35 bar needs a PVER-style
build-time clamp of the interpolation at the table edge (pressure broadening
frozen at its 10-bar value — defensible: deeper layers are LW-opaque and barely
sunlit, and the albedo rise is Rayleigh-driven) — deliberately deferred until
these results are reviewed.
