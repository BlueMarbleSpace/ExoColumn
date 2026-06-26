# Planetary-mass HZ reference case

ExoColumn version of Kopparapu et al. (2014, ApJL 787, L29): the dependence of
the habitable-zone (HZ) limits on **planetary mass** (0.1, 1, 5 M⊕).  A more
massive planet has higher surface gravity, which compresses the H₂O (inner
edge) and CO₂ (outer edge) column depths, so the inner edge moves *inward*
(higher Seff) and the outer edge changes little.  The target is their **Figure
3**: HZ limits vs stellar effective temperature for the three masses.

## What sets the mass dependence

Two quantities, both following Kopparapu+2014's prescription, anchored so that
**1 M⊕ = Earth exactly** (recovering ExoColumn's validated inner/outer HZ
reference and the multi-stellar Figure-6/7 Sun point):

1. **Surface gravity `g`.**  ExoRT bakes gravity in at *compile* time
   (`exo_g → SHR_CONST_G`; the radiation core's column amount is `pdel/g`, and
   the cold-start adiabat / hydrostatic heights use `exo_g`).  So **each mass
   needs its own binary**, built with the Makefile's `EXO_G` override
   (`make PVER=200 EXO_G=<g>`; default `9.80616` is bit-identical to the stock
   Earth build).
2. **Background N₂ pressure `p_N2`.**  Scaled with mass per their Eq. (3)
   (case 3, "N₂ scaled with planet radius").  This is a *runtime* namelist
   value, threaded into `hz_inner`/`hz_outer.run_one(n2_bar=…)`.

Mass–radius relation (Kopparapu+2014, fit to exoplanets.org):

    M/M⊕ = 0.968 (R/R⊕)^3.2          (M < 5 M⊕)
  ⟹  g(M)    = g⊕ · M^(1 − 2/3.2) = g⊕ · M^0.375
      p_N2(M) = 1 bar · M^(2.40/3.2) = 1 bar · M^0.75

| M [M⊕] | R/R⊕ | g [m/s²] | p_N2 [bar] |
|:------:|:----:|:--------:|:----------:|
|  0.1   | 0.49 |   4.135  |   0.178    |
|  1.0   | 1.01 |   9.806  |   1.000    |
|  5.0   | 1.67 |  17.931  |   3.344    |

(Kopparapu's empirical fit coefficients 0.968 / 0.937 put `g(1 M⊕)` ~2% below
Earth; anchoring 1 M⊕ to Earth removes that offset while keeping their mass
*trend* — the exact powers M^0.375, M^0.75 are independent of the
normalisation.)

## HZ limits

- **Inner edge = runaway greenhouse** — the Simpson–Nakajima peak of the
  inner-edge `Seff(Ts)` curve (the same definition standardised in
  `reference/moist_runaway` / `hz_figure7`).  Kopparapu+2014 adopt the runaway
  rather than the moist-greenhouse limit (the two differ by <2%, and the
  runaway limit is less sensitive to the assumed 200 K tropopause).
- **Outer edge = maximum greenhouse** — the `Seff` minimum of the dense-CO₂
  sweep (1 → 34.7 bar at Ts = 273 K); nearly mass-independent.

Each point reuses the validated single-column configs of
`reference/moist_runaway/hz_inner.run_one` (inner) and
`reference/max_greenhouse/hz_outer.run_one` (outer), only changing the host-star
SED, the binary's gravity, and `n2_bar`.

## Files

- `hz_mass.py` — mass parameterization + build/sweep driver + Figure-3 plot.
- `sweep_all_masses.sh` — thin wrapper that builds & sweeps all three masses,
  restores the Earth binary, and plots.
- `hz_mass_m{0.1,1,5}.npz` — per-mass cached sweep results (written by the sweep).
- `hz_mass.{pdf,png}` — the figure (PDF is the manuscript copy).

## Reproduce

The build steps need the Intel OneAPI runtime; source it first:

```bash
source /opt/intel/oneapi/setvars.sh
bash reference/planet_mass/sweep_all_masses.sh      # ≈ build+sweep ×3 + restore + plot
```

Equivalent / piecewise (run from the project root):

```bash
python reference/planet_mass/hz_mass.py all          # n68 M->G sweep ×3 masses
python reference/planet_mass/hz_mass.py addf          # + n84 F 7200 K endpoint ×3 masses
python reference/planet_mass/hz_mass.py gravity 5     # -> 17.93134 (EXO_G for 5 M⊕)
python reference/planet_mass/hz_mass.py sweep 5       # n68 sweep the CURRENT binary as 5 M⊕
python reference/planet_mass/hz_mass.py addf1 5       # merge F 7200 K into the 5 M⊕ cache
python reference/planet_mass/hz_mass.py plot          # re-plot from caches (or HZ_REPLOT=1)
```

The driver leaves the default **Earth-gravity** binary in place (`run/exocol.exe`
rebuilt at `EXO_G=9.80616`) so the repo's validated reference is restored.

## Host-star set & scope

The sweep uses the n68 cool/solar ladder (M 2600 K → G Sun) shared with
`reference/habitablezone/hz_figure6.py`, **plus the hot F 7200 K endpoint on the
n84 core** (its BT-Settl SED exists only at n84, which also resolves the strong
F-star shortwave ~2% better; a separate n84 build per mass, as in
`hz_add_n84_stars.py`).  So each mass curve combines n68 (M→G) points with the
n84 F point, spanning Kopparapu's full 2600–7200 K range.  The `all` mode runs
the n68 ladder; `addf` then appends the n84 F point to each cache.  (A/B hotter
stars give unphysical albedo > 1 in the 1-D two-stream shortwave and are
excluded — see `hz_add_n84_stars.py`.)  The Kopparapu+2014 Table-1 parametric
fit (their Eq. 4) is drawn as the dashed reference across the same range.

## Expected result (at the Sun)

Kopparapu+2014's headline inner-edge (runaway) numbers are `Seff = 0.99 / 1.107
/ 1.188` for 0.1 / 1 / 5 M⊕ (~10% lower / ~7% higher flux than Earth).
ExoColumn's 1 M⊕ inner edge already sits ~5% above Kopparapu's at the Sun (the
documented near-IR shortwave H₂O albedo offset; see
`reference/moist_runaway/README.md` and the LBL benchmarks), so the three
ExoColumn curves are expected to lie a few percent inward of the dashed
Kopparapu curves while reproducing the **mass ordering and spacing**.
