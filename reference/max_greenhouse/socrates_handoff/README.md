# Outer-HZ maximum-greenhouse column — profile handoff

This is the **exact atmospheric column** used in the right-hand panel of Figure 2
of the ExoColumn paper (the clear-sky radiative-transfer benchmark: ExoRT
`n68equiv` vs `Clima` vs a from-scratch line-by-line calculation).

It is the S_eff-minimum ("maximum greenhouse") point of the outer-HZ sweep,
i.e. the ExoColumn analogue of Kopparapu et al. (2013) Figure 5.

## The case in one line

**Ts = 273.0 K (fixed), 8.87 bar CO2 over 1 bar N2 (ps = 9.875 bar), trace H2O,
154 K isothermal stratosphere, clear sky, Earth gravity.**

This is an **inverse** (prescribed-profile) calculation, not a radiative–convective
equilibrium solution: the T(p) profile is an *input*, and only the fluxes are
diagnosed. So SOCRATES should be run in single-call flux mode on exactly the
profile in these files — no relaxation, no adjustment.

## Files

| file | contents |
|---|---|
| `ohz_maxgh_profile.txt` | Human-readable. Full setup header, then Table 1 = 201 level (interface) p/T/z, Table 2 = 200 layer (midpoint) p/T/vmr, Table 3 = layer column amounts as a cross-check. |
| `ohz_maxgh_profile.csv` | Same 200 layers, machine-readable, one row per layer with both bounding interfaces. |
| `ohz_maxgh_band_olr.txt` | Band-resolved TOA OLR on the ExoRT n68 band grid: ExoRT, and the LBL reference under both far-wing conventions. This is what the figure plots. |
| `clima_band_olr_maxgh.txt` | The same quantity from Clima, on its own 55-interval band grid (two k-coefficient generations). |

Index convention throughout: **1 = model top, last = surface.**

## Conventions worth stating explicitly

- `vmr` = mole fraction of **total (moist)** air, so CO2 + N2 + H2O = 1 exactly
  at every layer. H2O is trace here (<= 2.2e-4 by mass), so CO2 and N2 are
  essentially constant with height (x_CO2 = 0.8982 at the surface, 0.8987 at the top).
- Pressures in Pa, temperatures in K, heights in m above the surface.
- Gravity g = 9.80616 m/s^2, held constant (no variation with height).
- The 200-layer grid is log-spaced from 0.01 Pa to the surface. Feel free to
  coarsen — for the line-by-line benchmark we mass-weighted it down to 40 layers
  and the OLR moved by < 0.1 W/m^2.

## How the profile was constructed

- H2O-saturated (RH = 1) moist adiabat integrated up from Ts = 273 K.
- Pinned to the CO2 saturation curve wherever the ascent supersaturates in CO2
  (Kasting 1991; Span & Wagner 1996 saturation curve), with the CO2 share of cp
  evaluated at the local temperature.
- CO2 is held **well mixed** — CO2 condensation adjusts the temperature profile
  but the condensate is not removed and CO2 clouds are neglected, as in Clima.
- Capped by an isothermal 154 K stratosphere.

## Longwave result to compare against (0-2000 cm^-1 OLR)

| model | OLR (W/m^2) |
|---|---|
| LBL, pure-Lorentz CO2 wings (opaque bound) | 44.1 |
| Clima, Kopparapu-2013 k-coefficients | 69.4 |
| Clima, Wolf HITRAN-2016 k-coefficients | 72.3 |
| **ExoRT n68equiv (this work)** | **75.5** |
| LBL, Perrin & Hartmann (1989) sub-Lorentzian chi (transparent bound) | 78.5 |

The three band/k models agree to within 6 W/m^2, and ExoRT is 3.0 W/m^2 (3.9%)
below the line-by-line reference.

Per n68 band, ExoRT and the LBL agree to |delta| <= 0.25 W/m^2 across
720-1108 cm^-1. Essentially all of the remaining residual is a single band,
**1108-1200 cm^-1 (LBL 5.34 vs ExoRT 2.14 W/m^2)** — which is exactly where the
HITRAN-2024 CO2-CO2 CIA tabulation has a gap: it covers 1-750 and 1150-1850
cm^-1 with nothing in between, and ExoRT's own GB/Baranov CIA is likewise ~0
across 875-1150. So that band is the one place worth looking closely at in a
SOCRATES comparison.

Note (2026-08-21): an earlier version of this table quoted 91.7 W/m^2 for the
PH89 LBL. That was an artefact of a too-aggressive weak-line cutoff in our own
LBL code, which discarded 17540 of the 17610 CO2 lines in the 8-12 um window on
this dense column. It is fixed; the LBL total moved 91.7 -> 78.5 and the window
17.6 -> 14.7. The atmospheric profile itself is unaffected.

## Shortwave, if you also want the S_eff comparison

The figure itself is longwave-only, but for completeness the SW side used:

- Present-day solar spectrum (G2V), S0 = 1360 W/m^2 at 1 AU -> TOA down 339.99 W/m^2.
- Surface albedo 0.32, Lambertian and spectrally grey (direct and diffuse).
- 6-point Gauss-Legendre hemispheric solar-zenith average, insolation normalised to S0/4.

ExoRT gives TOA SW up = 148.83 W/m^2 (planetary albedo 0.438), hence
S_eff = OLR / (SW_dn - SW_up) = 0.395, i.e. d = 1/sqrt(S_eff) = 1.59 AU
(Kopparapu et al. 2013: 0.343, 1.70 AU).

## Regenerating

```
python3 reference/max_greenhouse/make_socrates_handoff.py
```

reads `reference/max_greenhouse/exocol_maxgh_8.87bar.nc` (the archived ExoColumn
output for this column) and rewrites everything here.
