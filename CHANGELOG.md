# Changelog

All notable changes to ExoColumn are documented here. Versions follow
[semantic versioning](https://semver.org): the major version changes when
results change, the minor version when capabilities are added.

## [1.0.0] — 2026-08-28

First public release, accompanying Haqq-Misra, Wolf & Kopparapu (submitted to
*The Astrophysical Journal*), *"Validation of Habitable Zone Limits with a New
Radiative-Convective Equilibrium Climate Model"*.

### Model

- One-dimensional radiative–convective equilibrium column coupled to the
  [ExoRT](https://github.com/storyofthewolf/ExoRT) correlated-*k* radiation core
  (`aerad_driver`), which is used without modification. Compile-time settings
  ExoColumn owns (vertical levels, surface gravity) are applied by generating
  local copies of a few ExoRT files at build time.
- Convection: simplified Betts–Miller (default), dry adjustment,
  Manabe–Wetherald, and Zhang–McFarlane schemes, plus a surface-coupled mixed
  layer and a non-deadlocking dry fallback.
- Surface: prognostic slab with a simplified Monin–Obukhov flux scheme
  (resolution-independent 10 m reference height) or a legacy bulk-aerodynamic
  scheme; selectable wet or dry surface.
- Moisture: prognostic water vapour with Newton saturation adjustment and a
  cold-point (Brewer–Dobson) stratospheric cold trap.
- Thermodynamics: native IAPWS-95 water equation of state, the Wagner–Pruß
  saturation curve, the Kasting (1988) non-ideal moist pseudoadiabat, and a
  two-component ideal-gas pseudoadiabat below the triple point.
- CO2: Span & Wagner saturation curve, CO2 condensation for the outer habitable
  zone, and a temperature-dependent CO2 specific heat.
- Runtime options for the H2O continuum (MT_CKD or BPS), host-star spectrum,
  multi-point solar-zenith quadrature, and cold-trap phase.

### Validation and reproducibility

- `reference/` contains eight self-contained cases, each with its data and the
  script that regenerates the corresponding published figure: present-day Earth
  against konrad and Clima; line-by-line OLR benchmarks at both habitable-zone
  edges; inner- and outer-edge limits; multi-stellar boundaries for twelve host
  stars; and habitable-zone boundaries for 0.1, 1, and 5 Earth-mass planets.
- `test/` contains unit tests for the IAPWS-95, steam-adiabat, and CO2 modules.

### Release engineering

- Added `CITATION.cff`, `.zenodo.json`, and this changelog.
- The build no longer hardcodes a machine-specific Intel compiler path: `ifx` is
  taken from `PATH` by default, and the compiler vendor is detected from the
  basename of `USER_FC`, so absolute or versioned `gfortran` paths are now
  classified correctly.
- The `EXORT_ROOT` environment variable overrides the ExoRT location in the
  Python tooling, matching the build's `config.mk` setting.
- `data/exort_extra/` ships the three BT-Settl stellar spectra (2000, 2200 and
  2400 K) that the multi-stellar figures need and that are not yet part of the
  public ExoRT distribution.
