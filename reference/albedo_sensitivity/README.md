# Supplementary case — sensitivity of the HZ limits to the surface albedo

Raised in coauthor review (2026-08-21): the HZ figures originally fixed the
surface albedo at **α_s = 0.32**, the value \citet{kopparapu2013} tuned inside
*Clima* so their 1-D Earth reaches 288 K. ExoColumn's own Earth calibration gives
**α_s = 0.2736**. Does using our tuned value move the derived limits toward the
*Clima* values?

**Outcome: yes — and on the strength of this result α_s = 0.2736 was adopted as
the primary albedo for every HZ calculation in the paper** (matching Kopparapu's
tuning *procedure* rather than his number). All HZ reference cases now carry the
α_s = 0.2736 results under their canonical filenames; the α_s = 0.32 set is
archived beside them under the `_a032` tag. This directory documents the
comparison that motivated the change.

## Why this is a clean experiment

The HZ limits are computed in **inverse (`flux_only`) mode**: the `T(p)` and
`H₂O(p)` columns are prescribed by the pseudoadiabat, so they cannot respond to
the surface albedo. Surface albedo therefore acts **only on the shortwave**.
The generator asserts this on the cached sweeps and it holds exactly:

```
check [mtckd]: max|dOLR| = 0.000e+00 W/m2 , max|d(strat H2O VMR)| = 0.000e+00
check [bps]  : max|dOLR| = 0.000e+00 W/m2 , max|d(strat H2O VMR)| = 0.000e+00
check [OHZ]  : max|dOLR| = 0.000e+00 W/m2
```

So `F_IR`, the H₂O profiles, and the moist-greenhouse trigger temperature
(`Ts = 345 K`) are **bit-identical** between the two albedos; the entire effect is
a rescaling of `F_SOL`, `α_p` and `S_eff = F_IR/F_SOL`. That also makes the
result portable: a single slope `dS_eff/dα_s` (reported below) lets any other
tuned albedo be scaled without re-running the model.

## Key results

| limit | α_s = 0.32 | α_s = 0.2736 | Δ | dS_eff/dα_s | Clima |
|---|---|---|---|---|---|
| moist greenhouse (MT_CKD) | 1.087 | 1.063 | −0.024 | 0.51 | 1.016 |
| moist greenhouse (BPS)    | 1.074 | 1.052 | −0.022 | 0.48 | 1.016 |
| runaway greenhouse (MT_CKD) | 1.101 | 1.071 | −0.030 | 0.65 | 1.060 |
| runaway greenhouse (BPS)    | 1.093 | 1.062 | −0.032 | 0.68 | 1.060 |
| maximum greenhouse | 0.395 | 0.385 | −0.010 | 0.21 | 0.343 |

Corresponding distances `d = 1/√S_eff` [AU]:

| limit | α_s = 0.32 | α_s = 0.2736 | Clima |
|---|---|---|---|
| moist greenhouse (MT_CKD / BPS) | 0.959 / 0.965 | 0.970 / 0.975 | 0.992 |
| runaway greenhouse (MT_CKD / BPS) | 0.953 / 0.956 | 0.966 / 0.970 | 0.971 |
| maximum greenhouse | 1.591 | 1.611 | 1.707 |

Max-greenhouse pCO₂ shifts only 8.92 → 8.68 bar (Clima ~8 bar).

1. **Lowering α_s moves every limit toward *Clima***, as expected — more absorbed
   SW at fixed `F_IR` means a lower `S_eff`, i.e. both edges move outward.
2. **The runaway (inner) limit is where it matters.** With *both* conventions
   matched to Kopparapu — our tuned albedo **and** the BPS continuum they used —
   ExoColumn gives 1.062 vs *Clima*'s 1.060 (0.970 vs 0.971 AU): agreement to
   0.002 in `S_eff`, 0.001 AU. The albedo convention alone closes ~73 % of the
   MT_CKD gap; albedo + continuum together close ~94 %.
3. **The moist-greenhouse limit does not follow.** Its gap narrows from 0.071 to
   0.047 (MT_CKD) — the residual is the stratospheric-H₂O/line-data difference
   already diagnosed in `reference/moist_runaway/README.md`, not the albedo.
4. **The outer edge is nearly albedo-insensitive** (`dS_eff/dα_s` = 0.21, a third
   of the inner-edge value): at the 8.7-bar maximum-greenhouse point the CO₂
   Rayleigh layer hides the surface, so the two `α_p` curves converge as pCO₂
   rises (bottom-left panel). It closes only ~19 % of the outer-edge gap.

**Interpretation.** Neither albedo is "more correct". 0.32 is *Clima*'s cloud
proxy, tuned inside *Clima*; 0.2736 is ours, tuned inside ExoColumn — and tuned
in a *different configuration* than the HZ sweeps use (RCEMIP O₃, 400 ppm CO₂,
N₂/O₂/Ar Earth air, full RCE, vs pure N₂ + 330 ppm CO₂, no O₃, a saturated
inverse column and 6-point zenith quadrature). Holding α_s at Kopparapu's value
in the main text isolates the **radiative-transfer** differences between the two
models; this supplement quantifies how much of the residual `S_eff` offset is
instead attributable to the **albedo convention** — most of it, at the inner edge.

## Files

| file | contents |
|---|---|
| `plot_albedo_sensitivity.py` | Figure generator. Runs **no** model — reads the six sweep caches, asserts the SW-only invariance above, prints the tables, draws the figure. |
| `albedo_sensitivity.{pdf,png}` | The supplementary figure (2×2: inner-edge `α_p` and `S_eff`; outer-edge `α_p` and `S_eff`). |

Sweep caches live with their own reference cases and are **not** duplicated here.
Primary (α_s = 0.2736): `reference/moist_runaway/hz_inner_nonideal[_bps].npz`,
`reference/max_greenhouse/hz_outer.npz`.  Archived (α_s = 0.32): the same names
with the `_a032` tag.

## Reproducing

The α_s = 0.2736 sweeps (~30 min inner + ~5 min outer, binary at `PVER=200`):

The primary (α_s = 0.2736) sweeps are just the scripts' defaults.  To regenerate
the archived α_s = 0.32 comparison set (~30 min inner + ~5 min outer, binary at
`PVER=200`):

```bash
source /opt/intel/oneapi/setvars.sh
HZ_ALBEDO=0.32  HZ_TAG_SUFFIX=_a032  python reference/moist_runaway/hz_inner.py
OHZ_ALBEDO=0.32 OHZ_TAG_SUFFIX=_a032 python reference/max_greenhouse/hz_outer.py
python reference/albedo_sensitivity/plot_albedo_sensitivity.py
```

Always pair `HZ_ALBEDO`/`OHZ_ALBEDO` with a tag suffix so a variant writes its own
files instead of overwriting the primary ones. Run the two sweeps
**sequentially**: both drive the same `exocol_config.nml` and
`iofiles/exocol_out.nc` scratch files.

## Caveat for anyone re-deriving the limits

`hz_inner.py` takes the runaway limit as the Simpson–Nakajima peak of `S_eff`
over `RUNAWAY_TS_LO/HI = 280–700 K` for MT_CKD, and reads the **BPS** limits off
the BPS curve at the *same* two `Ts`. That second step is not just a
convenience. Taking an independent peak for the BPS curve is **not robust**: at
α_s = 0.2736 the BPS plateau is shallow enough that the 280–700 K maximum
migrates to 690 K — the supercritical branch climbing past `Tc = 647.1 K`, not
the Simpson–Nakajima peak — and would report 1.076 instead of 1.062. The
published α_s = 0.32 curves are window-insensitive (identical over 280–700 /
500 / 450 K), so the **main-text numbers are unaffected**; the fragility only
appears in variants with a less-depressed plateau.
