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
| `HZ_TS_MIN` / `HZ_TS_MAX` / `HZ_TS_STEP` | 200 / 2500 / 5 | Restrict/refine the `Ts` grid (e.g. a quick validation subset). |
| `HZ_SMOOTH` | 5 | Running-median window (points) suppressing the resolution sawtooth; set 1 for the raw curve. |

## Key results

The figure overlays three curves per panel (a)–(c): ExoColumn **MT_CKD** (solid),
ExoColumn **BPS** continuum (dotted), and **Kopparapu/CLIMA** (dashed); panel (d)
overlays the ExoColumn and CLIMA H₂O profiles. The figure uses the `'ice'`
cold-trap convention — the model default *and* CLIMA's actual convention (see
Configuration).

- **Greenhouse-limit `Seff`** (this figure, `'ice'` cold trap + sub-273 K
  two-component adiabat):
  | | moist GH | runaway GH (Tc=647 K) |
  |---|---|---|
  | ExoColumn MT_CKD | 1.088 (Ts=345 K) | 1.088 |
  | ExoColumn BPS | 1.075 (Ts=345 K) | 1.101 |
  | Kopparapu/CLIMA | 1.018 | 1.060 |
  (Pure-N₂ background. Note Kopparapu's moist-GH limit derives from their
  grid-inflated stratospheric H₂O — see the panel-(d) bullet below — so their
  water-loss threshold is crossed at a cooler Ts than a grid-converged CLIMA
  would give.)
- **Continuum (BPS vs MT_CKD):** BPS lowers the planetary albedo toward CLIMA,
  closing ~half the panel-(b) albedo offset and ~20 % of the moist-GH `Seff`
  gap; **OLR is continuum-insensitive**, so the residual `Seff`/OLR offset is
  near-IR H₂O **line** data (ExoRT n68/HITRAN-2016 k-distribution vs CLIMA),
  not the continuum.
- **H₂O profiles (panel d):** fully apples-to-apples as of 2026-06-11. Both
  models integrate the same physics: two-component moist pseudoadiabat (steam
  tables / IAPWS-95 above 273.16 K, matching <0.1 %; ice-saturated two-component
  ideal-gas form below it, `twocomp_dlnTdlnP` ≡ CLIMA `convec.f` label 13,
  T(P) matching ±0.07 K), capped by an isothermal 200 K stratosphere and
  freeze-dried at the cold trap over **ice**. Surface VMRs match CLIMA to
  0.3–1 % at every profiled Ts (the dry partial pressure is 1.0005 bar in both).
  The plotted VMR uses the exact two-component mmr→VMR conversion (the dilute
  `q·mw/18` form read ~25 % high at the 380 K surface) and includes the true
  z = 0 surface point. **Known residual:** Kopparapu's *tabulated* profiles
  freeze-dry at their cold trap's grid level, one coarse layer (ΔlnP ≈ 0.10–0.19)
  below the 200 K cap — i.e. at T* = 202–206 K for `Ts ≤ 340 K` — which inflates
  their stratospheric H₂O by es_ice(T*)/es_ice(200 K) ≈ 1.3–2.1×. ExoColumn
  freeze-dries at the *interpolated* 200 K crossing (grid-snapped cold traps are
  what caused the panel (a)–(c) staircase). Accounting for that factor, the two
  models agree to 2–3 % at **every** profiled Ts (280–380 K); the raw tabulated
  ratios are ~0.45–0.75 at Ts ≤ 340 K and 0.92–0.98 at 360/380 K. This grid
  artifact also means Kopparapu's water-loss (moist-GH) limit was diagnosed from
  slightly inflated stratospheric H₂O.

## Caveats

- **k-table ceiling:** ExoRT's n68 k-tables have a hardcoded 500 K temperature
  ceiling (read-only `radgrid.F90`). The high-`Ts` branch (≳1600 K, the
  post-runaway OLR rise) uses clamped opacity and is qualitatively right but
  quantitatively extrapolated — do not over-interpret absolute values there.
- **Sawtooth:** the small high-frequency ripple in panels (a)–(c) is a
  vertical-resolution sampling artifact (the tropopause kink migrating across
  fixed log-pressure layers as `ps` grows), not physics; it converges as ~1/N
  with `PVER` and the median smoother removes the residual.
