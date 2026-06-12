# Reference case: moist- / runaway-greenhouse inner HZ (non-ideal water EOS)

ExoColumn analogue of **Kopparapu et al. (2013), Figure 3** — the inner edge of
the habitable zone for a G2V (Sun-like) star, swept over surface temperature
`Ts`. This is the **non-ideal-EOS** result we validate against; the earlier
ideal-gas variant was a stepping stone and has been retired.

For each `Ts` the model cold-starts a fully-saturated moist column (isothermal
stratosphere capped at `t_strato = 200 K`), calls ExoRT once in `flux_only` mode,
and records OLR, absorbed SW, planetary albedo, `Seff = OLR/ASR`, and the H₂O
profile. The figure has four panels (Kopparapu Fig 3):

- **(a)** OLR and absorbed SW vs `Ts`
- **(b)** Planetary albedo vs `Ts`
- **(c)** Effective stellar flux `Seff` vs `Ts` (with Kopparapu + ExoColumn moist-/runaway-greenhouse limits)
- **(d)** H₂O volume-mixing-ratio profiles vs altitude for `Ts ∈ {280…380 K}`, clipped at 100 km

## Files

| File | What it is |
|------|------------|
| `hz_inner.py` | Self-contained sweep + plot generator. The full per-run namelist is embedded in `NML_TEMPLATE`, so the case is reproducible from source. |
| `hz_inner_nonideal.npz` | Cached sweep results (461 `Ts` points + 6 H₂O profiles). Enables instant re-plot. |
| `hz_inner_nonideal.pdf` / `.png` | The figure. **The PDF is the publication artifact** — hand-edit it as needed. |
| `waterloss_IHZ_present.dat` | Kopparapu et al. (2013) inner-HZ sweep vs surface temperature (provided directly by R. Kopparapu). Columns: `TGO  SEFF  PALB  FH2O  FTIR(1)`[OLR]`  FTSO(1)`[absorbed SW]. Overlaid on panels (a)–(c). |
| `clima_last.tab` | Kopparapu CLIMA water-loss vertical profiles — one ALT/P/T/FH2O/… block per surface temperature (220, 240, … K). Col 1 = altitude [km], col 4 = H₂O VMR. The six blocks matching the panel-(d) ExoColumn profiles (280–380 K) are overlaid on panel (d). |

The Kopparapu reference data are drawn as thin **dashed** curves, **colour-matched** to the
corresponding ExoColumn quantity (solid = ExoColumn, dashed = Kopparapu) — in panel (d)
that means each CLIMA profile shares the colour of the ExoColumn profile at the same
surface temperature. A style legend in panel (a) states the solid/dashed convention. If
either `.dat`/`.tab` file is absent the overlay is silently skipped.

## Configuration (fixed in `hz_inner.py`)

- **Water EOS:** `h2o_eos = 'nonideal'` — Kasting (1988) Appendix-A non-ideal moist
  pseudoadiabat on the native IAPWS-95 EOS (forces Wagner–Pruß steam saturation).
- **Composition:** pure N₂ (1 bar) + CO₂ = 3.3e-4 + prognostic H₂O; no O₂/Ar/O₃/CH₄ —
  matching Kopparapu (2013)'s IHZ "Earth" model (Fig 3: N₂ background, FCO₂ = 3.3e-4).
  (Earlier runs carried O₂ = 0.21 + Ar = 0.01; dropping them shifts Seff by only
  ~+0.003 but is the apples-to-apples choice and is consistent with the OHZ sweep.)
- **Surface albedo:** 0.32 (Kopparapu's IHZ value, Fig 3b, mimicking cloud reflection).
- **Solar zenith:** 6-point Gauss–Legendre **hemispheric** quadrature (`sw_zenith_quad`,
  `sw_nquad=6`) for the flux-weighted Bond average — same *scheme* as Kopparapu's
  6-angle Gaussian average. (Our nodes are GL-in-μ → zenith 14.9–88.1°; Kopparapu's
  are 11.0–82.8°; the integrand is the same flux-weighted average, so this is a small
  quadrature-node difference, and Kopparapu's weights are not published for an exact match.)
- **`variable_ps`:** `ps = p_N2(1 bar) + esat(Ts)`; **`msdist = 1 AU`**.
- **CO₂ mixing convention:** `co2_vmr_total = .true.` — CO₂ is 3.3×10⁻⁴ of the
  **total (moist)** air at every layer, Kopparapu/CLIMA's convention (their
  FCO₂ = const at all levels, verified in `clima_last.tab`). For a steam-rich
  column this means the CO₂ amount **grows** with the vapour inventory
  (×~200 in absolute terms at Ts → Tc where p ≈ 220 bar) — unphysical for a
  fixed CO₂ inventory (the model default holds CO₂/N₂ constant) but required
  for apples-to-apples: it is radiatively negligible at the moist-GH end
  (−0.02 W/m² OLR at 300 K) yet lowers the steam-plateau OLR by ~7 W/m² and
  the runaway-limit Seff from 1.088 to **1.060**. Yang et al. (2016) flag the
  same convention caveat for their own intercomparison.
- **Model top:** `p_top = 0.002 Pa` — so even the coldest profiled column (280 K)
  reaches ≥100 km in panel (d) (at 0.01 Pa it stopped at ~96 km). SW is converged
  wrt the top (mass above 0.002 Pa is ~1e-8 of the column).
- **Cold-trap saturation phase:** `cold_trap_phase = 'ice'` (the model default).
  Verified against the actual CLIMA source (atmos repo, 2026-06-11): below the
  273.16 K triple point CLIMA saturates over **ice** everywhere (`SATRAT`: CC
  with the sublimation latent heat) and its sub-freezing moist adiabat
  (`convec.f` label 13) is the ice-saturated two-component pseudoadiabat —
  which ExoColumn reproduces via `exocol_convadj::twocomp_dlnTdlnP` when
  `h2o_eos='nonideal'` (sub-273 K T(P) agrees with `clima_last.tab` to
  **±0.07 K** at Ts = 380 K). An earlier revision of this figure used
  `'liquid'` under the belief that supercooled liquid was CLIMA's convention;
  its apparent agreement at the hot cases came from two compensating errors
  (wetter liquid esat × the then too-steep textbook-`malr` sub-freezing
  fallback) and it over-watered the 300 K stratosphere ~1.5×.
- **`Ts` grid:** 200–2500 K, 5 K step (env-overridable; see below).

## Reproducing

The sweep needs an ExoColumn binary **built at `PVER ≥ 200`** (the wide
log-pressure span from `p_top = 0.002 Pa` would otherwise coarsen the
troposphere). From the project root:

```bash
source /opt/intel/oneapi/setvars.sh
cd /hugespace/models/ExoColumn/build && make clean && make PVER=200
```

### Full re-sweep (~14 min, regenerates the cache + figure)

```bash
source /opt/intel/oneapi/setvars.sh        # exe runtime libs
python reference/moist_runaway/hz_inner.py
```

The script finds the binary (`run/exocol.exe`) and scratch I/O
(`iofiles/`, `exocol_config.nml`) via paths relative to its own location, and
writes the figure + `.npz` cache back into this directory. It temporarily
writes `exocol_config.nml` per run and restores the original on exit.

### Re-plot only (instant — for figure/label/style tweaks)

```bash
HZ_REPLOT=1 python reference/moist_runaway/hz_inner.py
```

Loads `hz_inner_nonideal.npz` and regenerates the figure without any model runs.

### Useful environment overrides

| Var | Default | Effect |
|-----|---------|--------|
| `HZ_REPLOT` | (unset) | If set, re-plot from the caches; skip the sweep. |
| `HZ_BPS` | 1 | Overlay the BPS-continuum sweep (dotted) on panels (a)–(c) for a direct Kopparapu-continuum comparison; set 0 for the MT_CKD-only figure. Doubles the sweep time (a second radiation pass per `Ts`). |
| `HZ_COLD_TRAP` | ice | Sub-freezing saturation phase. `ice` (this figure) is both the model default and CLIMA's actual convention; `liquid` is a sensitivity variant (~2× wetter cold stratosphere). |
| `HZ_TRAP_EMU_SWEEP` | 0 | Verification mode: apply the Kopparapu cold-trap sampling to **every** sweep run (δ(Ts) interpolated between the six measured anchors). Verified 2026-06-12: the moist-GH and runaway limits are **identical to quoted precision** (1.087/1.060 MT_CKD, 1.074/1.072 BPS); max\|ΔSeff\| = 0.0011 at 300–400 K (0.004 worst-case at 235 K, where ΔOLR ≤ 0.8 W/m²). Pair with `HZ_TAG_SUFFIX` to keep the published caches untouched. |
| `HZ_TRAP_EMU` | 1 | Re-run the six panel-(d) profiles with `coldtrap_dT_offset` set to Kopparapu's measured per-case cold-trap grid offsets (T*−200 = 0.13–5.82 K from `clima_last.tab`), reproducing their tabulated stratospheric H₂O to 0–2% (sampling-faithful overlay). Panels (a)–(c) and all limits always use the model's own interpolated-200 K cold trap. Set 0 to plot the model's own cold trap in panel (d) (their curves then read ×1.3–2.1 high at Ts ≤ 340 K). |
| `HZ_TS_MIN` / `HZ_TS_MAX` / `HZ_TS_STEP` | 200 / 2500 / 5 | Restrict/refine the `Ts` grid (e.g. a quick validation subset). |
| `HZ_SMOOTH` | 5 | Running-median window (points) suppressing the resolution sawtooth; set 1 for the raw curve. |

## Key results

The figure overlays three curves per panel (a)–(c): ExoColumn **MT_CKD** (solid),
ExoColumn **BPS** continuum (dotted), and **Kopparapu/CLIMA** (dashed); panel (d)
overlays the ExoColumn and CLIMA H₂O profiles. The figure uses the `'ice'`
cold-trap convention — the model default *and* CLIMA's actual convention (see
Configuration).

- **Greenhouse-limit `Seff`** (this figure: `'ice'` cold trap + sub-273 K
  two-component adiabat + CLIMA CO₂-of-total convention):
  | | moist GH | runaway GH (Tc=647 K) |
  |---|---|---|
  | ExoColumn MT_CKD | 1.087 (Ts=345 K) | 1.060 |
  | ExoColumn BPS | 1.074 (Ts=345 K) | 1.072 |
  | Kopparapu/CLIMA | 1.016 | 1.060 |
  (Pure-N₂ background. With the **physical** fixed-CO₂-inventory convention
  the runaway limit is 1.088/1.101 — the CO₂-of-total convention's growing
  CO₂ inventory lowers the steam-plateau OLR by ~7 W/m². The MT_CKD runaway
  value landing exactly on Kopparapu's 1.060 is partly compensating offsets:
  our plateau OLR is ~7 W/m² **below** theirs under the matched convention
  while our absorbed SW is ~2.5% below theirs (albedo offset). Kopparapu's
  moist-GH limit derives from their grid-inflated stratospheric H₂O — see the
  panel-(d) bullet — so their water-loss threshold is crossed at a cooler Ts
  than a grid-converged CLIMA would give.)
- **Continuum (BPS vs MT_CKD):** BPS lowers the planetary albedo toward CLIMA,
  closing ~half the panel-(b) albedo offset and ~20 % of the moist-GH `Seff`
  gap; **OLR is continuum-insensitive**, so the residual `Seff`/OLR offset is
  near-IR H₂O **line** data (ExoRT n68/HITRAN-2016 k-distribution vs CLIMA),
  not the continuum.
- **SW line-by-line verdict (2026-06-12, `tools/lbl_sw_benchmark.py`):** the
  SW twin of the OLR benchmark — RADIS/HITRAN H₂O+CO₂ (0.5–5 µm) + the MT_CKD
  continuum, ExoRT's own Rayleigh formulas at band midpoints, a Toon-89
  **quadrature** two-stream (ExoRT's solar scheme; validated to ±0.0006 band
  albedo in the gas-free vis/UV bands and ±10⁻⁵ energy conservation), the
  same 6-node zenith quadrature, and ExoRT's own per-band incident fluxes as
  stellar weights.  Result at Ts = 300 K: **planetary albedo LBL = 0.2730 vs
  ExoRT n68 = 0.2733 (+0.0002)**; absorbed SW within 0.1 W/m².  So the n68
  k-distribution reproduces line-by-line SW absorption essentially exactly
  *given its ingredients* (HITRAN lines + MT_CKD continuum), and the old
  "n68 under-absorbs near-IR H₂O" interpretation is **revised**: Kopparapu's
  albedo (0.257 at 300 K, ~0.016 below the LBL) reflects their *ingredient*
  choices — the BPS continuum absorbs more near-IR than MT_CKD (our band
  BPS toggle closes about half the gap, consistently) — plus residual
  solar-spectrum/resolution differences, not an ExoRT k-table error.
  Combined with the LW verdict below: **both radiation cores of this model
  are LBL-grade at the bell peak; the remaining Seff offsets vs Kopparapu
  (2013) are attributable to their 2013-era opacity ingredients.**
- **LW (F_IR) offset vs Kopparapu (diagnosed 2026-06-11):** ΔOLR(Ts) is a bell:
  +4 W/m² at 220 K, peaking at **~+20 W/m² at 300–320 K**, collapsing through
  zero near 400 K — i.e. our F_IR *rises faster* toward the plateau. Under the
  matched CO₂-of-total convention our steam-plateau OLR then sits ~7 W/m²
  **below** theirs (284–285 vs 291.8 W/m² at 500–1300 K; with the physical
  fixed-CO₂-inventory convention the plateaus coincide at ~291.7 — a
  coincidence of our LW transparency against their growing CO₂ opacity). With T(P)
  and the H₂O profile now verified to match CLIMA, this is a pure
  radiation-core (gas-opacity) difference in the **semi-transparent** regime.
  Spectral home (band-resolved, Ts = 300 K): the 8–12.5 µm window is 67–82 %
  transparent in n68 (~117 W/m² of OLR) and the 380–800 cm⁻¹
  H₂O-rotation/CO₂-wing complex 43–54 %; the Ts = 220 K limit (an essentially
  H₂O-free N₂+CO₂ column, still +3.7 W/m²) bounds the CO₂-side share at
  ~4 W/m², leaving ~16 W/m² on the H₂O side (window + rotation-band wings).
  Ruled out: the continuum choice (BPS vs MT_CKD ≤ 0.4 W/m² per band), the
  T/q state (matched), and the cold-trap sampling (≤ 0.7 W/m²).
  **Line-by-line verdict (2026-06-12, `tools/lbl_olr_benchmark.py`):** an
  independent LBL calculation on the *same* Ts = 300 K column (RADIS/HITRAN
  lines for H₂O+CO₂ at 0.01 cm⁻¹ + a faithful port of the AER MT_CKD
  continuum, diffusivity-1.66 Schwarzschild) gives **OLR = 272.4 W/m²
  (10–3000 cm⁻¹) vs ExoRT n68 = 269.7 — agreement to 1 %**, with per-band
  residuals ≤ 0.35 W/m² everywhere except the CO₂-wing/rotation-band overlap
  (−1.7/−1.5/−0.5 W/m² in 720–800/545–617/667–720, ExoRT slightly *more*
  opaque than LBL).  Kopparapu's 250.2 W/m² is **22 W/m² below the LBL** — the
  size of an entire second MT_CKD continuum (the whole continuum removes
  23.7 W/m² from the lines-only 296.4) and far beyond plausible
  continuum-model spread.  Conclusion: the F_IR bell is **CLIMA-2013 being
  too opaque relative to modern line-by-line radiative transfer, not an
  ExoRT deficiency** — consistent with Yang et al. (2016), where the ExoRT
  lineage (CAM4_Wolf) tracks the SMART LBL while band models spread
  10–25 W/m², and with Kopparapu's ≈1420 W/m² runaway threshold vs
  Goldblatt/SMART's ≈1340.  LBL caveats (all ≲ 1 W/m² here): no N₂–N₂ CIA,
  air- (not N₂-) broadened widths, no CO₂ χ-factors at 330 ppm.  This
  reframes the moist-GH `Seff` gap: ~60 % of it (the OLR ratio 1.047 at
  340 K) is CLIMA's LW opacity bias, the rest the SW/albedo offset (ASR
  ratio 1.028) — i.e. our limit is the better-grounded one by the LBL
  standard, and "agreement with Kopparapu" should not be pursued further on
  the LW side.
- **H₂O profiles (panel d):** fully apples-to-apples as of 2026-06-11. Both
  models integrate the same physics: two-component moist pseudoadiabat (steam
  tables / IAPWS-95 above 273.16 K, matching <0.1 %; ice-saturated two-component
  ideal-gas form below it, `twocomp_dlnTdlnP` ≡ CLIMA `convec.f` label 13,
  T(P) matching ±0.07 K), capped by an isothermal 200 K stratosphere and
  freeze-dried at the cold trap over **ice**. Surface VMRs match CLIMA to
  0.3–1 % at every profiled Ts (the dry partial pressure is 1.0005 bar in both).
  The plotted VMR uses the exact two-component mmr→VMR conversion (the dilute
  `q·mw/18` form read ~25 % high at the 380 K surface) and includes the true
  z = 0 surface point. **Cold-trap sampling:** Kopparapu's *tabulated* profiles
  freeze-dry at their cold trap's grid level, one coarse layer (ΔlnP ≈ 0.10–0.19)
  below the 200 K cap — i.e. at T* = 200.1–205.8 K — which inflates their
  stratospheric H₂O by es_ice(T*)/es_ice(200 K) ≈ 1.0–2.1× over the true
  cold-trap value. The default figure (`HZ_TRAP_EMU=1`) therefore re-runs the
  six profile cases with `coldtrap_dT_offset` = their measured per-case T*−200,
  sampling our continuous adiabat at the same temperature: the overlay then
  matches their tabulated stratospheric H₂O to **0–2 % at every profiled Ts**
  (and the surface to 0.1–1 %). ExoColumn's own cold trap remains the
  *interpolated* 200 K crossing (grid-snapped traps caused the panel (a)–(c)
  staircase), and panels (a)–(c) + all quoted limits use it exclusively — the
  emulation's radiative effect is < 0.7 W/m² in OLR (at 280 K) and < 0.003 in
  Seff, i.e. the IHZ limits are **insensitive to the cold-trap sampling
  convention**; the moist-GH/runaway gaps vs Kopparapu are entirely the
  radiation offsets (LW: see below; SW: near-IR H₂O absorption/albedo).

## Caveats

- **k-table ceiling:** ExoRT's n68 k-tables have a hardcoded 500 K temperature
  ceiling (read-only `radgrid.F90`). The high-`Ts` branch (≳1600 K, the
  post-runaway OLR rise) uses clamped opacity and is qualitatively right but
  quantitatively extrapolated — do not over-interpret absolute values there.
- **Sawtooth:** the small high-frequency ripple in panels (a)–(c) is a
  vertical-resolution sampling artifact (the tropopause kink migrating across
  fixed log-pressure layers as `ps` grows), not physics; it converges as ~1/N
  with `PVER` and the median smoother removes the residual.
