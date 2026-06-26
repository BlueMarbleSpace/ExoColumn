# ExoColumn

A one-dimensional radiative–convective equilibrium (RCE) model for planetary
atmospheres, built directly on the [ExoRT](https://github.com/storyofthewolf/ExoRT)
correlated-k radiative transfer code (the radiation core of ExoCAM).

ExoColumn time-marches a single atmospheric column to radiative–convective
equilibrium by calling ExoRT's `aerad_driver` each step, with a prognostic
surface slab, surface turbulent fluxes, moist convection, condensation, and a
stratospheric cold-trap closure. It is intended for habitable-zone and
terrestrial-atmosphere studies where the full 3-D GCM is unnecessary.

## Status

- **Earth validation**: surface temperature 287.8 K, agreeing with the
  independent RCE codes [konrad](https://github.com/atmtools/konrad) (288.0 K)
  and CLIMA/Kasting (289.5 K); energy and water budgets closed. See
  `reference/earth/`.
- **Habitable zone**: inner- and outer-edge calculations reproducing
  Kopparapu et al. (2013), including the moist/runaway greenhouse inner edge,
  the maximum-greenhouse outer edge, and the multi-stellar (F–G–K–M) boundaries.
  See `reference/moist_runaway/`, `reference/max_greenhouse/`,
  `reference/habitablezone/`.

## Requirements

- A copy of **ExoRT** (default location `/models/ExoRT`; override with
  `EXORT_ROOT`). ExoRT is treated as read-only.
- **NetCDF-Fortran** built with the same compiler as ExoColumn.
- A Fortran compiler: **Intel `ifx`** (OneAPI; the default) or **`gfortran`**.

> The system NetCDF on most Linux distributions is built with gfortran and is
> incompatible with `ifx`. For an `ifx` build you need a NetCDF-Fortran compiled
> with `ifx` (see the build notes in `CLAUDE.md`). On Apple Silicon/macOS the
> system NetCDF works directly with `gfortran`.

## Build

Copy `config.mk.example` to `config.mk` and set your local paths
(`NETCDF_ROOT`, `EXORT_ROOT`, optionally `USER_FC`/`PVER`). Then:

```bash
source /opt/intel/oneapi/setvars.sh      # Intel OneAPI (ifx build)
make -C build                            # → run/exocol.exe
make -C build clean
```

The number of vertical layers is a compile-time choice owned by ExoColumn
(not ExoRT). The default is 70; override on the command line:

```bash
make -C build PVER=200
```

A `make clean` is required after changing `PVER`. The build never modifies
ExoRT: it generates local copies of the handful of ExoRT files whose
compile-time settings ExoColumn controls.

## Run

Run from the project root so ExoRT can resolve its data-file paths:

```bash
./run/exocol.exe
```

The run is configured by `exocol_config.nml` (three optional namelist blocks:
`&exocol_nml` runtime physics, `&exocol_init` initial conditions, and
`&exocol_composition` dry-air composition). With an empty `input_file` (the
default) the column is built from a cold start; otherwise an ExoRT-format
NetCDF profile is read. Calibrated example configurations are in `presets/`.

Output is written to `iofiles/exocol_out.nc` in ExoRT's `RTprofile_out.nc`
format.

## Plot

```bash
python tools/plot_exocol.py              # 4-panel summary of iofiles/exocol_out.nc
```

Each `reference/<case>/` directory is self-contained and includes its own
plotting script that regenerates the published comparison figure against the
reference data (konrad, CLIMA, Kopparapu et al. 2013).

## Repository layout

| Path | Contents |
|------|----------|
| `src/` | ExoColumn Fortran source |
| `build/` | Makefile and build-time ExoRT patch tooling |
| `tools/` | Plotting, input-generation, and line-by-line benchmark scripts |
| `test/` | Unit tests for the IAPWS-95, steam-adiabat, and CO2 modules |
| `presets/` | Calibrated example namelists |
| `reference/` | Self-contained validation/HZ cases with data and figures |
| `CLAUDE.md` | Detailed architecture and developer notes |

## License

MIT License © 2026 Blue Marble Space. See `LICENSE`.
