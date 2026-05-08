# ExoColumn — Project Brief for Claude Code

## Overview

ExoColumn is a new 1-D radiative-convective equilibrium (RCE) climate model written in
Fortran. It wraps the ExoRT radiative transfer (RT) model developed by Eric T. Wolf
(University of Colorado), calling ExoRT's core computation routines directly rather than
through file-based I/O. The result is a self-contained, compiled RCE column model suitable
for planetary atmosphere research, including exoplanet and early Earth applications.

This document is the authoritative design brief for the initial implementation. Claude Code
should treat it as the source of truth for architecture, module responsibilities, and
development order.

---

## Background

### What is an RCE model?

A radiative-convective equilibrium model finds the vertical temperature profile T(z) at
which the atmosphere is in steady state under two competing processes:

- **Radiation**: shortwave (SW) heating and longwave (LW) cooling at each layer, computed
  by a radiative transfer scheme.
- **Convection**: when the radiative lapse rate exceeds the adiabatic lapse rate, convective
  adjustment redistributes heat upward to restore stability.

Equilibrium is reached when net heating at every level (radiative + convective) goes to
zero and the TOA energy budget closes.

### What is ExoRT?

ExoRT (https://github.com/storyofthewolf/ExoRT) is a two-stream, correlated-k radiative
transfer code written in Fortran (~93% of the codebase). It is designed for use with the
NCAR CESM/ExoCAM 3-D GCM but also supports offline 1-D column calculations.

Key ExoRT facts relevant to ExoColumn:

- The **computational kernel** is the subroutine `aerad_driver`. It takes a single
  atmospheric column (T, p, gas profiles, surface properties) and returns upwelling and
  downwelling LW and SW fluxes at every interface level, plus heating rates at every
  layer midpoint.
- In CESM, `aerad_driver` is called by `exo_radiation_tend`, which sets up column
  profiles from GCM data structures. ExoColumn needs its own equivalent wrapper.
- The recommended spectral version for terrestrial planets is `src.n68equiv` (68 spectral
  intervals, 8 Gauss points, HITRAN 2016 line data, correlated-k via HELIOS-K).
- ExoRT uses **NetCDF** for its standalone offline I/O (`RTprofile_in.nc` /
  `RTprofile_out.nc`), but ExoColumn bypasses this and calls `aerad_driver` directly in
  memory.

### ExoRT input variables (what aerad_driver needs per column)

These are the quantities ExoColumn must maintain in its column state and pass to the
radiation wrapper:

| Variable   | Size    | Description                                  |
|------------|---------|----------------------------------------------|
| `ts`       | scalar  | Surface temperature (K)                      |
| `ps`       | scalar  | Surface pressure (Pa)                        |
| `tmid`     | (pver)  | Temperature at layer midpoints (K)           |
| `tint`     | (pverp) | Temperature at layer interfaces (K)          |
| `pdel`     | (pver)  | Pressure thickness of each layer (Pa)        |
| `pint`     | (pverp) | Pressure at interface levels (Pa)            |
| `zint`     | (pverp) | Height at interfaces (m)                     |
| `asdir`    | scalar  | SW albedo, direct beam                       |
| `asdif`    | scalar  | SW albedo, diffuse                           |
| `aldir`    | scalar  | Near-IR albedo, direct beam                  |
| `aldif`    | scalar  | Near-IR albedo, diffuse                      |
| `coszrs`   | scalar  | Cosine of solar zenith angle                 |
| `mw`       | scalar  | Molecular weight of dry air (g/mol)          |
| `cp`       | scalar  | Specific heat of dry air (J/kg/K)            |
| `h2ommr`   | (pver)  | H₂O specific humidity, kg(wv)/kg(air)       |
| `co2mmr`   | (pver)  | CO₂ mass mixing ratio (dry, optional)        |
| `ch4mmr`   | (pver)  | CH₄ mass mixing ratio (dry, optional)        |
| *(others)* | (pver)  | NH₃, CO, O₂, O₃, H₂, N₂, C₂H₆ (optional)  |

`pver` = number of layer midpoint levels; `pverp = pver + 1` (interface levels).

### ExoRT output variables (what aerad_driver returns)

| Variable | Size    | Description                             |
|----------|---------|-----------------------------------------|
| `LWUP`   | (pverp) | LW upwelling flux (W m⁻²)              |
| `LWDN`   | (pverp) | LW downwelling flux (W m⁻²)            |
| `SWUP`   | (pverp) | SW upwelling flux (W m⁻²)              |
| `SWDN`   | (pverp) | SW downwelling flux (W m⁻²)            |
| `LWHR`   | (pver)  | LW heating rate (K/Earth day)           |
| `SWHR`   | (pver)  | SW heating rate (K/Earth day)           |

---

## ExoColumn Architecture

### Design principles

1. **Direct subroutine coupling**: ExoColumn links against ExoRT source and calls
   `aerad_driver` directly. No intermediate NetCDF file I/O occurs inside the RCE loop.
2. **Fortran throughout**: Fully compiled Fortran for performance and natural interop with
   ExoRT's existing Fortran modules.
3. **Modular, single-responsibility files**: Each Fortran source file has a clearly bounded
   job (state, radiation, convection, loop, I/O, driver).
4. **NetCDF for external I/O only**: Initial conditions are read from NetCDF (or a namelist)
   before the loop; final equilibrium output is written to NetCDF after the loop.

### Call hierarchy (analog to CESM)

```
CESM:      radiation_tend → exo_radiation_tend → aerad_driver
ExoColumn: exocol_rce_loop → exocol_radiation   → aerad_driver
```

---

## Proposed Directory Structure

```
ExoColumn/
├── src/
│   ├── exocol_mod.F90          ! Column state: T, p, gas profiles; allocate/deallocate
│   ├── exocol_radiation.F90    ! Wrapper around aerad_driver (analog of exo_radiation_tend)
│   ├── exocol_convadj.F90      ! Convective adjustment (dry adiabat first; moist later)
│   ├── exocol_rce_loop.F90     ! Main RCE iteration loop + convergence check
│   ├── exocol_io.F90           ! Read initial conditions; write equilibrium output
│   └── exocol_driver.F90       ! Top-level PROGRAM: init, run loop, finalize
├── build/
│   └── Makefile                ! Links ExoColumn src + ExoRT src + NetCDF
├── iofiles/
│   └── exocol_in.nc            ! Input: initial T profile, gas concentrations, etc.
└── README.md
```

---

## Module Descriptions

### `exocol_mod.F90` — Column state

Defines all arrays and scalars that describe the current state of the atmospheric column.
This module is `USE`d by every other module. Key responsibilities:

- Declare `pver` and `pverp` as module-level parameters (set at init from input).
- Allocate/deallocate all column arrays (`tmid`, `tint`, `pint`, `pdel`, `zint`,
  `h2ommr`, `co2mmr`, etc.).
- Provide `exocol_init(pver)` and `exocol_finalize()` subroutines.
- Store scalar surface/atmospheric parameters (`ts`, `ps`, `mw`, `cp`, `coszrs`,
  albedos).

### `exocol_radiation.F90` — Radiation wrapper

Mirrors the role of `exo_radiation_tend` in CESM. Key responsibilities:

- Accept the current column state from `exocol_mod` and package it into exactly what
  `aerad_driver` expects.
- Call `aerad_driver`.
- Return `LWHR(pver)` and `SWHR(pver)` (heating rates) and the full flux arrays
  (`LWUP`, `LWDN`, `SWUP`, `SWDN`) at interfaces.
- **First task for Claude Code**: read `exo_radiation_tend.F90` from ExoRT's
  `source/src.main/` to extract the exact `aerad_driver` argument list and replicate
  the setup logic here. This is the most critical interface to get right.

### `exocol_convadj.F90` — Convective adjustment

Implements convective adjustment to enforce atmospheric stability. Key responsibilities:

- **Phase 1 (implement first)**: Dry adiabatic adjustment. Working from the surface
  upward, compare the local lapse rate to the dry adiabatic lapse rate
  `Γ_dry = g / cp`. If `Γ_actual > Γ_dry` between adjacent levels, redistribute
  enthalpy between those layers to restore `Γ_dry` while conserving total enthalpy
  in the column.
- **Phase 2 (later)**: Moist adiabatic adjustment. Compute the local moist adiabatic
  lapse rate from T(z) and q(z) and use it as the adjustment threshold instead.
- Interface: `subroutine convadj_dry(tmid, tint, pint, pdel, cp, g, pver)`

### `exocol_rce_loop.F90` — RCE iteration loop

The main solver. Key responsibilities:

- Time-march the column state forward using a large virtual timestep `dt` (recommended
  starting value: 1–10 Earth days per step, expressed in seconds).
- Each iteration:
  1. Call `exocol_radiation` → get `LWHR`, `SWHR`
  2. Update `tmid(k) ← tmid(k) + dt * (LWHR(k) + SWHR(k))` for all k
  3. Update `ts` using net surface flux: `ts ← ts + dt * F_net_surf / (rho_surf * cp_surf * dz_surf)`
     (or a simpler slab-ocean formulation)
  4. Recompute `tint` by interpolation from `tmid` and `ts`
  5. Call `convadj_dry` to enforce stability
  6. Check convergence
- Convergence criterion (implement both; stop when either is met):
  - `max(|LWHR(k) + SWHR(k)|) < 0.01` K/day across all levels
  - `|TOA net flux| = |SWDN_TOA - SWUP_TOA - LWUP_TOA| < 0.1` W/m²
- Cap iterations at a user-specified `nmax` (e.g. 100,000) to prevent runaway.
- Print a convergence diagnostic every N steps (e.g. every 1000 steps).

### `exocol_io.F90` — Input/output

Key responsibilities:

- `read_initial_conditions(filename)`: read `RTprofile_in.nc` (ExoRT-format input) or
  an ExoColumn-specific namelist to populate `exocol_mod` arrays. Reuse ExoRT's existing
  `makeColumn.py` tool to generate the input file.
- `write_output(filename)`: write the equilibrium state to NetCDF in ExoRT's
  `RTprofile_out.nc` format (so existing ExoRT plotting tools work unchanged).

### `exocol_driver.F90` — Top-level program

```fortran
PROGRAM exocol_driver
  USE exocol_mod
  USE exocol_io
  USE exocol_rce_loop
  implicit none
  call exocol_init(...)
  call read_initial_conditions('exocol_in.nc')
  call run_rce_loop()
  call write_output('exocol_out.nc')
  call exocol_finalize()
END PROGRAM exocol_driver
```

---

## Build System

The Makefile must compile ExoColumn's `src/` together with ExoRT's source (specifically
`source/src.main/` and `source/src.n68equiv/` for the recommended spectral version), and
link against NetCDF-Fortran.

Recommended compiler: `ifort` on Linux/HPC; `gfortran` on macOS Apple Silicon (ifort was
discontinued for arm64). ExoRT's existing Makefile in `build/` is the template — extend
it rather than starting from scratch.

---

## Key ExoRT Files to Review First

Before writing any ExoColumn code, Claude Code should read these ExoRT source files to
understand the `aerad_driver` interface:

1. `source/src.main/exo_radiation_tend.F90` — the CESM wrapper; defines what `aerad_driver` expects
2. `source/src.main/aerad_driver.F90` (or equivalent) — the actual kernel subroutine signature
3. `source/src.main/exoplanet_mod.F90` — module-level constants, solar spectrum selector
4. `build/Makefile` — understand how ExoRT compiles, then extend for ExoColumn

---

## Development Order (Recommended)

1. **Read ExoRT source** — understand `aerad_driver` argument list (see above).
2. **Write `exocol_mod.F90`** — all other modules depend on this data structure.
3. **Write `exocol_radiation.F90`** — verify it compiles and can call `aerad_driver`
   with a hardcoded test column before the loop exists.
4. **Write `exocol_io.F90`** — read an existing `RTprofile_in.nc` to populate the
   column state; confirm values are reasonable.
5. **Write `exocol_convadj.F90`** — dry adiabat; test independently with a synthetic
   super-adiabatic profile before integrating.
6. **Write `exocol_rce_loop.F90`** — wire radiation + convection together; run to
   convergence on a simple test case (e.g. Earth-like atmosphere, TS273K profile as
   initial condition).
7. **Write `exocol_driver.F90`** — top-level glue.
8. **Validate**: compare equilibrium T profile against ExoRT's standalone offline output
   for the same input column.

---

## Open Design Questions (to be resolved during development)

- **Surface energy balance**: Should `ts` be prognostic (updated by net surface flux each
  step) or held fixed? Prognostic is preferred; a simple "slab" with a specified heat
  capacity is a good starting point.
- **Interface temperature interpolation**: `tint` must be consistent with `tmid` and `ts`
  at every step. Confirm what interpolation scheme ExoRT's offline driver uses.
- **Timestep `dt`**: Start large (1–10 days) for speed; reduce if instabilities appear
  near the surface or tropopause.
- **Moist convection**: Defer to Phase 2. When implemented, will require carrying `h2ommr`
  as a prognostic variable and computing the moist adiabatic lapse rate at each level.
- **Parallelism**: Eric notes that `aerad_driver` is the parallelizable kernel (one call
  per column). Single-column RCE doesn't need this, but ensemble runs (e.g., sweeping
  over CO₂ concentrations) could trivially parallelize over columns with OpenMP or MPI.

---

## References

- ExoRT GitHub: https://github.com/storyofthewolf/ExoRT
- Wolf et al. (2022), PSJ 3:7 — primary ExoRT reference for n68equiv
- Wolf & Toon (2013), Astrobiology 13(7) — Archean applications
- Manabe & Wetherald (1967) — foundational RCE methodology
