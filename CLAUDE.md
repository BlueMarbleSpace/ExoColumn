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
  └── exocol_config      :: read_config (namelist: conv_scheme, moisture_scheme, wind_speed, C_D)
  └── exocol_mod         :: exocol_init, exocol_setgas, exocol_update_derived
  └── exocol_io          :: read_initial_conditions → populate exocol_mod
  └── ExoRT init sequence:: initialize_kcoeff → initialize_solar → init_ref
                            → init_model_specific → init_planck → initialize_radbuffer
  └── exocol_rce_loop    :: run_rce_loop (main iteration)
        ├── exocol_radiation :: exocol_rad_tend → aerad_driver
        ├── exocol_surface   :: compute_surface_fluxes (bulk aerodynamic LE, SH)
        └── exocol_convadj   :: convadj_dry | convadj_moist | convadj_manabe
  └── exocol_io          :: write_output (state + LE, SH, precip, cond_heating)
```

**Analogy to CESM/ExoCAM:**
```
CESM:      exo_radiation_tend → aerad_driver
ExoColumn: exocol_radiation   → aerad_driver
```

### Module responsibilities

- **`exocol_mod`** — Defines the entire column state (all arrays, scalars, and diagnostics `LE_diag`, `SH_diag`, `precip_diag`, `cond_heating`). `USE`d by every other ExoColumn module. `pver`/`pverp` come from ExoRT's compile-time `ppgrid` module (set by `exoplanet_mod::exo_pver`). Call order after init: `exocol_setgas()` → `exocol_update_derived()`.

- **`exocol_radiation`** — Wraps `aerad_driver`. Packages column state into the exact argument list `aerad_driver` expects. Converts heating rates from K/s (raw output) to K/day for the RCE loop.

- **`exocol_surface`** — Bulk-aerodynamic surface fluxes. `compute_surface_fluxes(ts, t_bot, q_bot, p_bot, mwdry, cpdry, U, C_D) → LE, SH` using `LE = ρ·L(Ts)·C_D·U·(qsat(Ts)−q_bot)` and `SH = ρ·cp·C_D·U·(Ts−T_bot)`. L is phase-aware: `L_v` for `Ts ≥ 273.16 K`, `L_sub` below. The rce loop applies an implicit-Euler damping factor `1/(1+dt/τ)` so the raw bulk formulas remain stable at the large virtual `dt` (τ ≈ 8.5 h vs `dt` = 5 d).

- **`exocol_config`** — Reads `exocol_config.nml` (namelist `&exocol_nml`). Exports `conv_scheme` ∈ {`'dry'`, `'moist'`, `'manabe'`}, `moisture_scheme` ∈ {`'prognostic'`, `'fixed_rh'`, `'off'`}, `wind_speed` (default 5 m/s), and `C_D` (default 1.5e-3). Silently uses defaults if the file is absent.

- **`exocol_rce_loop`** — Time-marches the column with a virtual timestep (`dt_days = 5` Earth days). Each step:
  1. radiation tendency on `tmid`
  2. bulk surface fluxes LE, SH (implicit-damped)
  3. slab budget: `ts += dt · (F_net_srf_rad − LE − SH) / H_slab`
  4. bottom-layer sources: `tmid(pver) += dt·SH/(cp·pdel/g)`, `h2ommr(pver) += dt·LE/(L·pdel/g)`
  5. `update_tint`
  6. convective adjustment (T + q mixed; see `exocol_convadj`)
  7. condensation cap: where `h2ommr(k) > qsat(T(k),p(k))`, set `h2ommr = qsat` and add `L(T)·q_excess/cp` to `tmid(k)` (phase-aware L)
  8. `update_derived`, `update_zint`

  Two convergence paths: **Path A** (radiative equilibrium): `max|LWHR+SWHR| < 0.01 K/day` AND `|TOA net flux| < 0.1 W/m²`. **Path B** (frozen-state stability): `max|Δtmid| < 0.001 K` AND `|ΔTs| < 0.001 K` over 100 steps, AND either `|TOA net flux| < 0.1 W/m²` OR `|ΔTOA flux| < 0.001 W/m²`. The `'fixed_rh'` and `'off'` moisture schemes are legacy code paths preserved for diagnostics (`fixed_rh` retains the historical RH relaxation closure with `tau_relax = 50 days`).

- **`exocol_convadj`** — Three schemes selectable via `conv_scheme`. All operate purely on atmosphere-atmosphere pairs (no surface-bottom pair adjustment — surface coupling is handled by the bulk SH/LE fluxes). All conserve `cp·T·Δp` in each adjusted pair. The `'moist'` and `'manabe'` schemes also mix `h2ommr` in adjusted pairs.
  - `'dry'`: potential-temperature stability criterion; q mixed fully (mass-weighted) when a pair is adjusted.
  - `'moist'`: rh-weighted local lapse rate `Γ_eff = rh·Γm + (1−rh)·Γd` where `rh = q/qsat` and `Γm = malr(T̄, p̄)` (saturated moist adiabat, phase-aware L). q-mixing is also rh-weighted: saturated pairs homogenize q; subsaturated pairs adjust only T.
  - `'manabe'`: fixed 6.5 K/km lapse rate (Manabe-Wetherald 1967); q mixed fully (mass-weighted).

  `esat_cc(T)` is phase-aware (Clausius-Clapeyron with L_v above 0 °C, L_sub below; continuous at 273.16 K). `Lvap_T(T)` returns the phase-appropriate latent heat and is reused by `exocol_surface` and the condensation step.

- **`exocol_io`** — NetCDF I/O using the Fortran 90 interface. Validates `pver`/`pverp` dimensions against compile-time constants on read. Output includes the column state plus diagnostics: `LE`, `SH` (W/m²), `precip` (mm/day), `cond_heating(pver)` (K/day from the final step's condensation).

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

`ts` is prognostic and the slab budget includes turbulent fluxes:

```
ts ← ts + dt · (F_net_srf_rad − LE − SH) / H_slab
H_slab = rho_w · cp_w · dz_slab = 1026 · 4000 · 50 = 2.052×10⁸ J/m²/K
```

`F_net_srf_rad = (SWDN − SWUP) + (LWDN − LWUP)` at the surface interface. `LE` and `SH` are bulk-aerodynamic fluxes from `exocol_surface` damped by `1/(1+dt/τ)` for stability (`τ = (pdel/g)/(ρ·C_D·U) ≈ 8.5 h` for Earth-like surface conditions; without damping the explicit Euler step overshoots by ~14×). LE and SH are also applied as bottom-layer sources for `tmid(pver)` and `h2ommr(pver)`, and the latent enthalpy enters the column as condensation heat where the vapor later saturates.

After updating `tmid` and `ts`, `tint(pverp)` is pinned to `ts` in `update_tint`. The convadj schemes operate purely on interior atmosphere-atmosphere pairs — there is no surface-bottom-layer pair adjustment (it would inject energy into the bottom layer without a matching slab debit, breaking column conservation).

### Moisture and energy conservation

With `moisture_scheme = 'prognostic'` the column conserves moist static energy in steady state: column-integrated `F_TOA − F_net_srf_rad + LE + SH = 0`. Mass-balanced steady state requires surface evaporation rate = column precipitation rate (verified diagnostically by comparing `LE/L_v` to `precip_diag`).
