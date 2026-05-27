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

### Setting vertical levels

`pver` (number of vertical layers) is compile-time, but owned by ExoColumn — not by ExoRT. At build time, the Makefile reads `/models/ExoRT/source/exoplanet_mod.F90`, substitutes our chosen `PVER`, and writes the result to `src/exoplanet_mod.F90`. The build compiles our local copy; ExoRT's source is never modified.

Default is **70 layers** (60 log-spaced upper layers + 10 near-surface geometric layers, giving a bottom midpoint at ~8.5 m altitude). Override via `config.mk` (`PVER = 80`) or on the command line:

```bash
make PVER=80
```

A `make clean && make PVER=N` is required after changing because every ExoRT static array depends on `pver`. The driver echoes the active value on startup: `pver = N layers (from exoplanet_mod::exo_pver)`.

## Running

The executable must be invoked from the **project root** so ExoRT can resolve its data file paths:

```bash
cd /hugespace/models/ExoColumn
run/exocol.exe
```

Initial conditions are controlled by `&exocol_init::input_file` in `exocol_config.nml`:
- **Empty** (default) → cold start. The column is built from the `&exocol_init` and `&exocol_composition` namelists: log-spaced pressure grid from `p_top` to `ps`, moist adiabat from `ts` capped at `t_strato`, `h2ommr = rh_init · qsat` in the troposphere, well-mixed dry gases. No input file required.
- **Non-empty** (e.g. `'iofiles/exocol_in.nc'`) → read the ExoRT-format NetCDF file (`RTprofile_in.nc` format).

Output: `iofiles/exocol_out.nc` (ExoRT `RTprofile_out.nc` format, compatible with ExoRT plotting tools)

Use ExoRT's `makeColumn.py` script to generate the input file.

## Plotting

After a run, inspect the output with:

```bash
python tools/plot_exocol.py                        # reads/writes iofiles/exocol_out.{nc,pdf}
python tools/plot_exocol.py my_in.nc my_out.pdf    # explicit paths
```

Produces a 4-panel PDF: temperature profile, radiative fluxes, heating rates, water vapour — all on a log-pressure axis. The temperature panel includes the US Standard Atmosphere 1976 (USSA-1976) as a dotted reference line.

## Architecture

ExoColumn is a 1-D radiative-convective equilibrium (RCE) model written in Fortran that directly calls ExoRT's `aerad_driver` subroutine. It does **not** use ExoRT's file-based I/O inside the RCE loop.

```
exocol_driver (PROGRAM)
  └── exocol_config      :: read_config (3 namelists: &exocol_nml, &exocol_init, &exocol_composition)
  └── exocol_mod         :: exocol_init, exocol_setgas, exocol_update_derived
  └── EITHER:
        exocol_coldstart :: cold_start_init  (when input_file = '')
      OR:
        exocol_io        :: read_initial_conditions  +
        exocol_config    :: apply_composition_overrides
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

- **`exocol_surface`** — Bulk-aerodynamic surface fluxes. `compute_surface_fluxes(ts, t_bot, q_bot, p_bot, mwdry, cpdry, U, C_D) → LE, SH` using `LE = ρ·L(Ts)·C_D·U·(qsat(Ts)−q_bot)` and `SH = ρ·cp·C_D·U·(Ts−T_bot)`. L is phase-aware: `L_v` for `Ts ≥ 273.16 K`, `L_sub` below. The rce loop applies an implicit-Euler damping factor `1/(1+dt/τ)` so the raw bulk formulas remain stable at the large virtual `dt` (τ ≈ 8.5 h vs `dt` = 1 d).

- **`exocol_coldstart`** — Self-contained initial-condition builder for cold-start runs (when `&exocol_init::input_file = ''`). Builds a log-spaced pressure grid from `p_top` to `ps`, integrates a moist adiabat upward from `ts` (`exocol_convadj::malr`) capped at `t_strato`, sets `h2ommr = rh_init·qsat` in the troposphere (zero above), broadcasts well-mixed dry-gas MMRs (Earth-like defaults overridden by `&exocol_composition`), and computes interface heights from hydrostatic balance. O3 is handled separately via `o3_profile`: `'uniform'` broadcasts the scalar, `'earth'` calls `exocol_ozone::set_earth_o3_profile`, `'none'` zeros it. Sets all surface scalars (`ts`, `coszrs`, albedos, `mwdry_col`, `msdist`) from the namelist. **`cpdry_col` is auto-computed** as the mass-weighted mean `Σ mmr_i·cp_i` from the dry-air composition (using CP_* constants in `exocol_config`); the namelist `cpdry` field is ignored. This ensures correct specific heat for non-Earth compositions (e.g., H₂-dominated atmospheres where cp differs by a factor of ~14). Calls `exocol_update_derived()` internally so the column is fully consistent on return.

- **`exocol_config`** — Reads `exocol_config.nml`. Three optional namelist blocks:
  - **`&exocol_nml`** (runtime physics): `conv_scheme` ∈ {`'dry'`, `'moist'`, `'manabe'`}, `moisture_scheme` ∈ {`'prognostic'`, `'fixed_rh'`, `'off'`}, `o3_profile` ∈ {`'uniform'`, `'earth'`, `'none'`} (default `'uniform'`; see `exocol_ozone`), `wind_speed` (default 5 m/s), `C_D` (default 1.5e-3), `msdist` (planet-star distance in AU, default 1.0; TOA stellar flux scales as 1/msdist²), `n_sfc_layers` (default 10; near-surface geometric layers), `dp_sfc_bot` (default 200 Pa; bottom layer pressure thickness, ~8.5 m midpoint altitude), `sfc_stretch` (default 1.5; geometric ratio between consecutive surface layers). Set `n_sfc_layers = 0` for pure log-spacing (legacy). `exocol_io::read_initial_conditions` copies `cfg_msdist` into `exocol_mod::msdist` after the input file is read (file mode); `exocol_coldstart::cold_start_init` does the equivalent for cold-start mode.
  - **`&exocol_init`** (initial conditions): `input_file` (default `''` → cold start), `ts`, `t_strato`, `p_top`, `rh_init`, `coszrs`, `cpdry` (ignored in cold start — see below), `asdir`, `asdif`, `aldir`, `aldif`. Selects between cold start (when `input_file` is empty) and reading an ExoRT NetCDF (when non-empty). All other fields supply the cold-start initial state and are ignored in file mode. Note: the Fortran namelist variables for `ts`, `coszrs`, `cpdry`, and the four albedos collide with identically-named `exocol_mod` state variables; consumers that need to USE both modules must rename one side (see `exocol_coldstart`).
  - **`&exocol_composition`** (dry-air composition + surface pressure): `ps` (Pa) plus per-gas volume mixing ratios `co2_vmr`, `ch4_vmr`, `c2h6_vmr`, `h2_vmr`, `n2_vmr`, `o3_vmr`, `o2_vmr`. **H2O is intentionally excluded** — water vapor is prognostic. Any field left unset (sentinel = -1) is taken from the input file in file mode, or from built-in Earth-like defaults in cold-start mode (N2=0.78, O2=0.21, CO2=4e-4, others=0). In file mode, applied by `apply_composition_overrides()` (called between `read_initial_conditions` and `exocol_setgas`). In cold-start mode, applied directly inside `cold_start_init`. Semantics (file mode):
    - **Per-gas VMR**: when ≥ 0, broadcast as a well-mixed scalar to every layer; when < 0, the input layer profile is kept and uniformly rescaled by `mwdry_in/mwdry_new` to preserve mass.
    - **`mwdry`**: recomputed as `Σ_dry VMR_i · Mw_i`. **`cpdry`** is also recomputed as `Σ mmr_i · cp_i` using the CP_* constants in `exocol_config` (CO2=844, CH4=2220, C2H6=1729, H2=14310, N2=1039, O3=820, O2=919 J/kg/K at ~300 K).
    - **`ps`**: when ≥ 0, `pmid`/`pint`/`pdel` are multiplied by `ps_new/ps_old`, preserving σ-structure.
    - A warning is printed if the resulting dry-air VMRs don't sum to ≈ 1.
  Silently uses defaults if `exocol_config.nml` is absent or any block is missing. Active settings are reported in the startup banner.

- **`exocol_rce_loop`** — Time-marches the column with an adaptive CFL-limited timestep (`dt = cfl_safety · dT_target / max|HR|`, clamped to `[dt_min, dt_max] = [1e-4, 1.0]` days; `dT_target = 2 K`). Each step:
  1. radiation tendency on `tmid`
  2. bulk surface fluxes LE, SH (implicit-damped)
  3. slab budget: `ts += dt · (F_net_srf_rad − LE − SH) / (H_slab + dt·4σTs³)` (semi-implicit Planck damping)
  4. bottom-layer sources: `tmid(pver) += dt·SH/(cp·pdel/g)`, `h2ommr(pver) += dt·LE/(L·pdel/g)`
  5. `update_tint`
  6. **iterated saturation adjustment** — `convadj` + `condense` are run in a fixed-point loop until no layer's T changes more than `sat_T_tol = 1e-4 K` between iterations (capped at `max_sat_iter = 50`). One `convadj` pass equilibrates the lapse rate; `condense` then removes excess vapor via implicit-Euler relaxation with timescale `τ_cond = 3600 s`: `q_new = (q + (dt/τ)·qsat) / (1 + dt/τ)`, releasing `Lvap_T(ts)·Δq/cp` into `tmid(k)`. `L` is evaluated at `ts` so the latent-heat ledger balances the surface evap debit; phase-aware `esat_cc` controls when condensation triggers. Latent release can re-destabilize adjacent pairs in convadj, so iterating to fixed point is required. `precip_total` and `cond_heating` accumulate over inner iterations.
  7. `update_derived`, `update_zint`

  Two convergence paths: **Path A** (radiative equilibrium): `max|LWHR+SWHR| < 0.01 K/day` AND `|TOA net flux| < 0.1 W/m²`. **Path B** (profile stability): time-mean quantities over 100-day windows must satisfy `max|Δtmid| < 1.0 K` AND `|ΔTs| < 0.02 K` AND `max|Δh2ommr| < 5e-4 kg/kg` AND `|⟨TOA⟩| < 0.1 W/m²`. The 1.0 K profile tolerance accounts for the irreducible noise floor (~0.3–0.6 K) from discrete convective bursts (τ_corr ~3–5 days, N_eff ~20–25 per 100-day window). The 5e-4 q tolerance accounts for condensation-relaxation flicker in 100-day means. The `'fixed_rh'` and `'off'` moisture schemes are legacy code paths preserved for diagnostics (`fixed_rh` retains the historical RH relaxation closure with `tau_relax = 50 days`).

- **`exocol_ozone`** — Provides `set_earth_o3_profile(pmid, mwdry, o3mmr)`, which interpolates a 15-point tabulated mid-latitude climatological O3 profile (Anderson et al. 1986 / SPARC, 1000–0.1 hPa) onto the model pressure grid and returns mass mixing ratios. Interpolation is linear in log(p) with clamping at table boundaries. Called by `exocol_coldstart` and `exocol_config::apply_composition_overrides` when `o3_profile = 'earth'`. Has no dependency on `exocol_config` (avoids circular USE). `aldir`/`aldif` in ExoRT are near-IR **solar** band albedos (0.7–5 μm), not thermal LW emissivity controls — the surface emits as a near-blackbody regardless of these values.

- **`exocol_convadj`** — Three schemes selectable via `conv_scheme`. All operate purely on atmosphere-atmosphere pairs (no surface-bottom pair adjustment — surface coupling is handled by the bulk SH/LE fluxes). All conserve `cp·T·Δp` in each adjusted pair. The `'moist'` and `'manabe'` schemes also mix `h2ommr` in adjusted pairs.
  - `'dry'`: potential-temperature stability criterion; q mixed fully (mass-weighted) when a pair is adjusted.
  - `'moist'`: rh-weighted local lapse rate `Γ_eff = rh·Γm + (1−rh)·Γd` where `rh = q/qsat` and `Γm = malr(T̄, p̄)` (saturated moist adiabat, phase-aware L). q-mixing is also rh-weighted: saturated pairs homogenize q; subsaturated pairs adjust only T.
  - `'manabe'`: fixed 6.5 K/km lapse rate (Manabe-Wetherald 1967); q mixed fully (mass-weighted).

  `esat_cc(T)` is phase-aware (Clausius-Clapeyron with L_v above 0 °C, L_sub below; continuous at 273.16 K). `Lvap_T(T)` returns the phase-appropriate latent heat and is reused by `exocol_surface` and the condensation step. `compute_tint_interp(tmid, pint, nv, tint)` fills `tint(1:nv)` by log-pressure interpolation/extrapolation and is the single canonical implementation used by all three schemes and by `exocol_rce_loop::update_tint`.

- **`exocol_io`** — NetCDF I/O using the Fortran 90 interface. Validates `pver`/`pverp` dimensions against compile-time constants on read. Output includes the column state plus diagnostics: `LE`, `SH` (W/m²), `precip` (mm/day), `cond_heating(pver)` (K/day from the final step's condensation).

- **`exocol_driver`** — Top-level `PROGRAM`. Reads namelists, branches on `cfg_input_file` to either cold-start (`exocol_coldstart::cold_start_init`) or file-read (`exocol_io::read_initial_conditions` + `exocol_config::apply_composition_overrides`), owns the ExoRT init sequence (mirrors `ExoRT/source/src.main/main.F90`), runs the RCE loop, writes output.

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

`F_net_srf_rad = (SWDN − SWUP) + (LWDN − LWUP)` at the surface interface. `LE` and `SH` are bulk-aerodynamic fluxes from `exocol_surface` damped by `1/(1+dt/τ)` for stability (`τ = (pdel/g)/(ρ·C_D·U) ≈ 8.5 h` for Earth-like surface conditions; at `dt = 1 day`, `dt/τ ≈ 2.8` so the raw explicit step would still overshoot ~3×). LE and SH are also applied as bottom-layer sources for `tmid(pver)` and `h2ommr(pver)`, and the latent enthalpy enters the column as condensation heat where the vapor later saturates.

After updating `tmid` and `ts`, `tint(pverp)` is pinned to `ts` in `update_tint`. The convadj schemes operate purely on interior atmosphere-atmosphere pairs — there is no surface-bottom-layer pair adjustment (it would inject energy into the bottom layer without a matching slab debit, breaking column conservation).

### Moisture and energy conservation

With `moisture_scheme = 'prognostic'` the column conserves moist static energy in steady state: column-integrated `F_TOA − F_net_srf_rad + LE + SH = 0`. Mass-balanced steady state requires surface evaporation rate = column precipitation rate (verified diagnostically by comparing `LE/L_v` to `precip_diag`).
