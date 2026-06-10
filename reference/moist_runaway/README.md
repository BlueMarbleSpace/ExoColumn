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
- **Composition:** N₂ = 0.78, O₂ = 0.21, Ar = 0.01, CO₂ = 3.3e-4, prognostic H₂O; no O₃/CH₄.
- **Surface albedo:** 0.32 (Kopparapu's cloud-free value, mimicking cloud reflection).
- **Solar zenith:** 6-point Gauss–Legendre **hemispheric** quadrature (`sw_zenith_quad`, `sw_nquad=6`), matching Kopparapu's 6-angle averaging.
- **`variable_ps`:** `ps = p_N2(1 bar) + esat(Ts)`; **`msdist = 1 AU`**.
- **Model top:** `p_top = 0.002 Pa` — so even the coldest profiled column (280 K)
  reaches ≥100 km in panel (d) (at 0.01 Pa it stopped at ~96 km). SW is converged
  wrt the top (mass above 0.002 Pa is ~1e-8 of the column).
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
| `HZ_REPLOT` | (unset) | If set, re-plot from the cache; skip the sweep. |
| `HZ_TS_MIN` / `HZ_TS_MAX` / `HZ_TS_STEP` | 200 / 2500 / 5 | Restrict/refine the `Ts` grid (e.g. a quick validation subset). |
| `HZ_SMOOTH` | 5 | Running-median window (points) suppressing the resolution sawtooth; set 1 for the raw curve. |

## Key results

- **ExoColumn moist-greenhouse `Seff` = 1.080** (at `Ts ≈ 350 K`, where the
  stratospheric H₂O VMR reaches the 3e-3 water-loss threshold).
- **ExoColumn runaway-greenhouse `Seff` = 1.099** (the plateau over 600 K … 1600 K).
- Kopparapu (2013) references: runaway 1.06, moist GH 1.015. The residual ~0.04
  in `Seff` is near-IR H₂O **shortwave** absorption (ExoRT n68 under-absorbs vs
  HITEMP-2010 + BPS continuum), seen as the panel-(b) albedo offset.

## Caveats

- **k-table ceiling:** ExoRT's n68 k-tables have a hardcoded 500 K temperature
  ceiling (read-only `radgrid.F90`). The high-`Ts` branch (≳1600 K, the
  post-runaway OLR rise) uses clamped opacity and is qualitatively right but
  quantitatively extrapolated — do not over-interpret absolute values there.
- **Sawtooth:** the small high-frequency ripple in panels (a)–(c) is a
  vertical-resolution sampling artifact (the tropopause kink migrating across
  fixed log-pressure layers as `ps` grows), not physics; it converges as ~1/N
  with `PVER` and the median smoother removes the residual.
