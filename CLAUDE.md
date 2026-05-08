# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Guiding principles

1. **ExoRT is read-only.** Never modify any file under `/models/ExoRT`.
2. **No duplication.** Before writing a constant, utility, or subroutine, check whether it already exists in ExoRT (`grep -r` in `/models/ExoRT/source`). If it does, `USE` it.
3. **Forward-compatibility.** Couple only to ExoRT's public module interfaces (`aerad_driver`, `physconst_setgas`, etc.), not to implementation internals that may change in a future ExoRT update.

## Build

### Prerequisites — NetCDF-Fortran for ifx

The system NetCDF (`/usr/lib64/gfortran/modules/netcdf.mod`) was compiled with gfortran and is **incompatible with ifx**. You must build and install NetCDF-C and NetCDF-Fortran with ifx before the first build. NetCDF-Fortran is already installed at `/opt/netcdf`. If it ever needs to be rebuilt (suggested install prefix: `/opt/netcdf`):

```bash
source /opt/intel/oneapi/setvars.sh

# 1. NetCDF-C (required by NetCDF-Fortran)
wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz
tar xf v4.9.2.tar.gz && cd netcdf-c-4.9.2
CC=icx ./configure --prefix=/hugespace/local/netcdf-ifx --disable-dap --disable-byterange
make -j$(nproc) && make install && cd ..

# 2. NetCDF-Fortran
wget https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz
tar xf v4.6.1.tar.gz && cd netcdf-fortran-4.6.1
FC=ifx CC=icx \
  CPPFLAGS=-I/hugespace/local/netcdf-ifx/include \
  LDFLAGS=-L/hugespace/local/netcdf-ifx/lib \
  ./configure --prefix=/hugespace/local/netcdf-ifx
make -j$(nproc) && make install
```

### Compiling ExoColumn

Activate Intel OneAPI, then build from the `build/` directory:

```bash
source /opt/intel/oneapi/setvars.sh
cd build
make                                                 # NETCDF_ROOT read from config.mk
make clean
```

To use gfortran (e.g. Apple Silicon macOS), the system NetCDF works directly:

```bash
make USER_FC=gfortran
```

## Running

The executable must be invoked from the **project root** so ExoRT can resolve its data file paths:

```bash
cd /hugespace/models/ExoColumn
run/exocol.exe
```

Input: `iofiles/exocol_in.nc` (ExoRT `RTprofile_in.nc` format)  
Output: `iofiles/exocol_out.nc` (ExoRT `RTprofile_out.nc` format, compatible with ExoRT plotting tools)

Use ExoRT's `makeColumn.py` script to generate the input file.

## Plotting

After a run, inspect the output with:

```bash
python tools/plot_exocol.py                        # reads/writes iofiles/exocol_out.{nc,pdf}
python tools/plot_exocol.py my_in.nc my_out.pdf    # explicit paths
```

Produces a 4-panel PDF: temperature profile, radiative fluxes, heating rates, water vapour — all on a log-pressure axis.

## Architecture

ExoColumn is a 1-D radiative-convective equilibrium (RCE) model written in Fortran that directly calls ExoRT's `aerad_driver` subroutine. It does **not** use ExoRT's file-based I/O inside the RCE loop.

```
exocol_driver (PROGRAM)
  └── exocol_io          :: read_initial_conditions → populate exocol_mod
  └── exocol_mod         :: exocol_setgas, exocol_update_derived
  └── ExoRT init sequence:: initialize_kcoeff → initialize_solar → init_ref
                            → init_model_specific → init_planck → initialize_radbuffer
  └── exocol_rce_loop    :: run_rce_loop (main iteration)
        ├── exocol_radiation :: exocol_rad_tend → aerad_driver
        └── exocol_convadj   :: convadj_dry
  └── exocol_io          :: write_output
```

**Analogy to CESM/ExoCAM:**
```
CESM:      exo_radiation_tend → aerad_driver
ExoColumn: exocol_radiation   → aerad_driver
```

### Module responsibilities

- **`exocol_mod`** — Defines the entire column state (all arrays and scalars). `USE`d by every other ExoColumn module. `pver`/`pverp` come from ExoRT's compile-time `ppgrid` module (set by `exoplanet_mod::exo_pver`). Call order after init: `exocol_setgas()` → `exocol_update_derived()`.

- **`exocol_radiation`** — Wraps `aerad_driver`. Packages column state into the exact argument list `aerad_driver` expects. Converts heating rates from K/s (raw output) to K/day for the RCE loop.

- **`exocol_rce_loop`** — Time-marches the column with a virtual timestep (`dt_days = 5` Earth days). Each step: radiation → update `tmid`/`ts` → recompute `tint` (log-p interpolation) → dry convective adjustment → update `zint` (hypsometric). Two convergence paths: **Path A** (radiative equilibrium): `max|LWHR+SWHR| < 0.01 K/day` AND `|TOA net flux| < 0.1 W/m²`. **Path B** (frozen-state stability): all three diagnostics — `tmid`, `ts`, and TOA net flux — change by less than 0.001 (K or W/m²) over 100 consecutive steps. Path B does **not** require the TOA flux itself to be small; a dry-only column reaches a genuine equilibrium with a structural TOA imbalance (~79 W/m² for Earth-like inputs) that Path B detects and accepts, printing a warning. Path B is necessary because a convectively active layer always has a large instantaneous HR (balanced by `convadj_dry`), so Path A is never satisfied there.

- **`exocol_convadj`** — Dry adiabatic adjustment using potential-temperature stability criterion. Sweeps surface→TOA, adjusting unstable pairs while conserving column enthalpy. Repeats up to 30 passes per timestep. Phase 2 (moist) is not yet implemented.

- **`exocol_io`** — NetCDF I/O using the Fortran 90 interface. Validates `pver`/`pverp` dimensions against compile-time constants on read. Output format is ExoRT-compatible.

- **`exocol_driver`** — Top-level `PROGRAM`. Owns the ExoRT init sequence (mirrors `ExoRT/source/src.main/main.F90`); does **not** call ExoRT's `input_profile`.

### ExoRT dependency

ExoRT source lives at `/models/ExoRT`. The Makefile compiles:
- `source/src.misc/` — shared utilities (kinds, constants, shr_*)
- `source/src.main/` — core radiation routines including `aerad_driver`
- `source/src.n68equiv/` — 68-band correlated-k spectral tables (recommended for terrestrial planets)

The key ExoRT files for understanding the `aerad_driver` interface are:
- `source/src.main/exo_radiation_mod.F90` — defines `aerad_driver`
- `source/src.main/exo_radiation_tend.F90` — CESM wrapper (the model for `exocol_radiation`)
- `source/exoplanet_mod.F90` — sets `exo_pver`, `exo_g`, solar spectrum selector

### Array indexing convention

Index 1 = TOA (top of atmosphere), index `pver`/`pverp` = surface. Midpoint arrays have size `pver`; interface arrays have size `pverp = pver + 1`.

### Derived pressure arrays

`pdeldry` and `pintdry` must be recomputed via `exocol_update_derived()` any time `pdel`, `pint`, or `h2ommr` change. These are not read from the input file; they are derived from wet quantities using `pdeldry(k) = pdel(k) * (1 - h2ommr(k))`.

### Surface energy balance

`ts` is prognostic: updated each step by `dt * F_net_srf / H_slab` where `H_slab = rho_w * cp_w * dz_slab = 1026 * 4000 * 50 = 2.052×10⁸ J/m²/K`. After updating `ts`, `tint(pverp)` is pinned to `ts` before convective adjustment.
