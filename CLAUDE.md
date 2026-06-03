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

Default is **70 layers**, **log-spaced** from `p_top` to `ps` (lowest level ~400 m). The near-surface coupling is handled by the surface-coupled mixed layer (`convadj_surface`), so no fine near-surface grid is needed. Override layer count via `config.mk` (`PVER = 80`) or on the command line:

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
        ├── exocol_surface   :: compute_surface_fluxes (Monin-Obukhov 'mos' | 'bulk' LE, SH)
        └── exocol_convadj   :: convadj_surface (surface mixed layer) → convadj_sbm
                                → convadj_dry (non-deadlocking fallback)
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

- **`exocol_surface`** — Surface turbulent fluxes, scheme chosen by `surface_flux`. **`'mos'` (default)**: simplified Monin-Obukhov (Frierson 2006) — potential-temperature/dry-static-energy differences with a neutral drag coefficient `C = κ²/ln²(z_ref/z0)` evaluated at a **fixed 10 m reference height** (`z_flux_ref`), not the lowest model level. Referencing C at 10 m makes the exchange coefficient resolution-independent (it no longer depends on where the bottom level falls) and Earth-magnitude (C ≈ 1e-3); paired with the surface-coupled mixed layer it gives an Earth-like climate (Ts ≈ 288, Bowen ≈ 0.24) on the plain log grid with **no boundary-layer mixing scheme**. **`'bulk'` (legacy)**: fixed-`C_D` bulk aerodynamic in actual temperature (`SH = ρ·cp·C_D·U·(Ts−T_bot)`, `LE = ρ·L(Ts)·C_D·U·(qsat(Ts)−q_bot)`); resolution-dependent, retained for comparison/regression (reproduces the calibrated Ts=288.01). L is phase-aware (`L_v` for `Ts ≥ 273.16 K`, `L_sub` below). Returns the exchange coefficient C. (Historical note: the Frierson K-profile boundary layer + hybrid fine-surface grid were prototyped and **removed** — the troposphere collapse they exhibited was actually the SBM gate, fixed by the dry-convective fallback; see `project_mos_no_bl_result`.)

- **`exocol_coldstart`** — Self-contained initial-condition builder for cold-start runs (when `&exocol_init::input_file = ''`). Builds a log-spaced pressure grid from `p_top` to `ps`, integrates a moist adiabat upward from `ts` (`exocol_convadj::malr`) capped at `t_strato`, sets `h2ommr = rh_init·qsat` in the troposphere (zero above), broadcasts well-mixed dry-gas MMRs (Earth-like defaults overridden by `&exocol_composition`), and computes interface heights from hydrostatic balance. O3 is handled separately via `o3_profile`: `'uniform'` broadcasts the scalar, `'earth'` calls `exocol_ozone::set_earth_o3_profile`, `'none'` zeros it. Sets all surface scalars (`ts`, `coszrs`, albedos, `mwdry_col`, `msdist`) from the namelist. **`cpdry_col` is auto-computed** as the mass-weighted mean `Σ mmr_i·cp_i` from the dry-air composition (using CP_* constants in `exocol_config`); the namelist `cpdry` field is ignored. This ensures correct specific heat for non-Earth compositions (e.g., H₂-dominated atmospheres where cp differs by a factor of ~14). Calls `exocol_update_derived()` internally so the column is fully consistent on return.

- **`exocol_config`** — Reads `exocol_config.nml`. Three optional namelist blocks:
  - **`&exocol_nml`** (runtime physics): `conv_scheme` ∈ {`'dry'`, `'moist'`, `'manabe'`, `'zm'`, `'sbm'`} (default `'sbm'`, the recommended simplified Betts-Miller scheme), `moisture_scheme` ∈ {`'prognostic'`, `'fixed_rh'`, `'off'`}, `o3_profile` ∈ {`'uniform'`, `'earth'`, `'none'`, `'rcemip'`} (default `'uniform'`; see `exocol_ozone`), `wind_speed` (default 5 m/s), `C_D` (default 1.5e-3; legacy `'bulk'` exchange coefficient), `msdist` (planet-star distance in AU, default 1.0; TOA stellar flux scales as 1/msdist²), `rh_sbm` (default 0.7; SBM reference relative humidity), `surface_flux` ∈ {`'mos'`, `'bulk'`} (default `'mos'`; see `exocol_surface`), `z0_rough` (default 3.21e-5 m; surface roughness used by `'mos'`), `latent_heat_mode` ∈ {`'phase_aware'`, `'fixed_vap'`}. `exocol_io::read_initial_conditions` copies `cfg_msdist` into `exocol_mod::msdist` after the input file is read (file mode); `exocol_coldstart::cold_start_init` does the equivalent for cold-start mode.
  - **`&exocol_init`** (initial conditions): `input_file` (default `''` → cold start), `ts`, `t_strato`, `p_top`, `rh_init`, `coszrs`, `cpdry` (ignored in cold start — see below), `asdir`, `asdif`, `aldir`, `aldif`. Selects between cold start (when `input_file` is empty) and reading an ExoRT NetCDF (when non-empty). All other fields supply the cold-start initial state and are ignored in file mode. Note: the Fortran namelist variables for `ts`, `coszrs`, `cpdry`, and the four albedos collide with identically-named `exocol_mod` state variables; consumers that need to USE both modules must rename one side (see `exocol_coldstart`).
  - **`&exocol_composition`** (dry-air composition + surface pressure): `ps` (Pa) plus per-gas volume mixing ratios `co2_vmr`, `ch4_vmr`, `c2h6_vmr`, `h2_vmr`, `n2_vmr`, `o3_vmr`, `o2_vmr`. **H2O is intentionally excluded** — water vapor is prognostic. Any field left unset (sentinel = -1) is taken from the input file in file mode, or from built-in Earth-like defaults in cold-start mode (N2=0.78, O2=0.21, CO2=4e-4, others=0). In file mode, applied by `apply_composition_overrides()` (called between `read_initial_conditions` and `exocol_setgas`). In cold-start mode, applied directly inside `cold_start_init`. Semantics (file mode):
    - **Per-gas VMR**: when ≥ 0, broadcast as a well-mixed scalar to every layer; when < 0, the input layer profile is kept and uniformly rescaled by `mwdry_in/mwdry_new` to preserve mass.
    - **`mwdry`**: recomputed as `Σ_dry VMR_i · Mw_i`. **`cpdry`** is also recomputed as `Σ mmr_i · cp_i` using the CP_* constants in `exocol_config` (CO2=844, CH4=2220, C2H6=1729, H2=14310, N2=1039, O3=820, O2=919 J/kg/K at ~300 K).
    - **`ps`**: when ≥ 0, `pmid`/`pint`/`pdel` are multiplied by `ps_new/ps_old`, preserving σ-structure.
    - A warning is printed if the resulting dry-air VMRs don't sum to ≈ 1.
  Silently uses defaults if `exocol_config.nml` is absent or any block is missing. Active settings are reported in the startup banner.

- **`exocol_rce_loop`** — Time-marches the column with an adaptive CFL-limited timestep (`dt = cfl_safety · dT_target / max|HR|`, clamped to `[dt_min, dt_max] = [1e-4, 1.0]` days; `dT_target = 1 K`). Each step:
  1. radiation tendency on `tmid`
  2. surface fluxes LE, SH (`exocol_surface`; Monin-Obukhov `'mos'` by default)
  3. slab budget: `ts += dt · (F_net_srf_rad − LE − SH) / (H_slab + dt·4σTs³)` (semi-implicit Planck damping)
  4. bottom-layer sources: `tmid(pver) += dt·SH/(cp·pdel/g)`, `h2ommr(pver) += dt·LE/(L·pdel/g)`
  5. `update_tint`
  6. **convective adjustment.** For `'sbm'`: `convadj_surface` (surface mixed layer) → `convadj_sbm` → `convadj_dry` (non-deadlocking fallback) → `condense` (Newton satadj). For the hard schemes (`dry`/`moist`/`manabe`): **thermodynamic consistency loop** — `convadj` + `condense` iterated together up to `max_inner_phys = 20` times until `condense` removes no vapour (`precip_iter == 0`), then one final `convadj` pass. Either way the column ends up simultaneously lapse-rate stable AND `q ≤ qsat` before each radiation step.
     - **`condense`** is a Newton satadj: `Δq = (q − qsat) / (1 + (L/cp)·dqsat/dT)` where `dqsat/dT = qsat·(1+qsat/ε)·L/(Rv·T²)`. CC convexity guarantees `q_new ≤ qsat(T_new)` in one step per layer. `L = Lvap_T(ts)` (surface temperature, not layer T) balances the surface evap ledger and prevents free `L_fusion` injection for vapor cycling through the ice phase. **No `tau_cond`** — the implicit-Euler relaxation scheme was removed.
     - **`sat_iter_convadj`** iterates the chosen `conv_scheme` until `max|ΔT| < sat_T_tol = 1e-4 K` (capped at `max_sat_iter = 50`).
  6b. **stratospheric cold-point cold trap** (`apply_stratospheric_coldtrap`, prognostic moisture only) — sets stratospheric water vapour by the Brewer-Dobson freeze-drying mechanism: air entering the stratosphere is dehydrated to the saturation mixing ratio at the cold-point tropopause (the coldest model level), and that value is conserved on ascent. Implemented as a **column-wide humidity floor** `h2ommr(k) = max(h2ommr(k), qsat(T_cp))`: freeze-dried air at the cold point is the driest air in the column, so the floor is a no-op in the moist troposphere but fills the otherwise-`q=0` stratosphere (both above the cold point and in the gap between the convective top and the cold point). Without it the stratosphere stays bone-dry (ExoColumn has no vertical moisture transport above the convective top), removing the principal stratospheric LW coolant and producing an unphysically warm stratosphere with a sharp single-layer tropopause where the cold convective top abuts the warm dry layer above. The floor is recomputed each step, so it is self-consistent: a colder cold point → drier stratosphere → less H₂O cooling. Same closure as konrad's `ColdPointCoupling`. After the fix the upper stratosphere matches konrad within ~1.4 K; a residual warm-cold-point bias (~204 K vs konrad's ~197 K) remains and is attributed to ExoRT-n68equiv vs RRTMG core radiation differences.
  7. `update_derived`, `update_zint`

  **Instantaneous TOA noise (~0.35 W/m² mean, ~1.3 W/m² peak near equilibrium) is expected and not a bug.** The mean comes from per-step radiation perturbation (`dt × max|HR| ≈ 0.06 d × 13 K/d ≈ 0.8 K`, OLR sensitivity ~0.5 W/m²/K → ~0.4 W/m²). The peaks come from discrete convective adjustment events, which set an irreducible noise floor of ~0.5–1.3 W/m² independent of `dT_target`. Only the 100-day mean must be < 0.1 W/m² for Path B.

  Two convergence paths: **Path A** (radiative equilibrium): `max|LWHR+SWHR| < 0.01 K/day` AND `|TOA net flux| < 0.1 W/m²`. **Path B** (profile stability): time-mean quantities over 100-day windows must satisfy `max|Δtmid| < 1.0 K` AND `|ΔTs| < 0.02 K` AND `max|Δh2ommr| < 5e-4 kg/kg` AND `|⟨TOA⟩| < 0.1 W/m²`. With `dT_target = 1 K`, the observed ΔTs_win at convergence is typically 0.003–0.007 K and ΔTpro < 0.04 K — well inside the tolerances. The `'fixed_rh'` and `'off'` moisture schemes are legacy code paths preserved for diagnostics (`fixed_rh` retains the historical RH relaxation closure with `tau_relax = 50 days`).

- **`exocol_ozone`** — Provides `set_earth_o3_profile(pmid, mwdry, o3mmr)`, which interpolates a 15-point tabulated mid-latitude climatological O3 profile (Anderson et al. 1986 / SPARC, 1000–0.1 hPa) onto the model pressure grid and returns mass mixing ratios. Interpolation is linear in log(p) with clamping at table boundaries. Called by `exocol_coldstart` and `exocol_config::apply_composition_overrides` when `o3_profile = 'earth'`. Has no dependency on `exocol_config` (avoids circular USE). `aldir`/`aldif` in ExoRT are near-IR **solar** band albedos (0.7–5 μm), not thermal LW emissivity controls — the surface emits as a near-blackbody regardless of these values.

- **`exocol_convadj`** — Five convection schemes selectable via `conv_scheme`, plus two always-on helpers used with `'sbm'` (`convadj_surface`, the surface-coupled mixed layer; and `convadj_dry` reused as a non-deadlocking fallback). The pairwise schemes (`dry`, `moist`, `manabe`, `zm`) operate purely on atmosphere-atmosphere pairs and conserve `cp·T·Δp` in each adjusted pair. The `'moist'` and `'manabe'` schemes also mix `h2ommr` in adjusted pairs. `sbm` is a whole-column scheme (see below).
  - **`convadj_surface`** (surface mixed layer): a **slab-rooted** dry convective adjustment — the slab (Ts, heat capacity `H_slab`) is the bottom node of a dry adjustment over a **fixed pressure depth** (`dp_surf_mix ≈ 150 hPa`). Super-adiabatic near-surface layers are mixed toward the dry adiabat rooted at Ts; the slab↔bottom pair conserves `H_slab·Ts + Σcp·T·Δp/g` (the slab is debited — this **is** the convective surface sensible flux, the RCE-appropriate replacement for a boundary-layer scheme). It is a **no-op on a stable surface layer** (e.g. the `'bulk'` reference, where θ_bottom > θ_surf), so it does not disturb that equilibrium. Run before `convadj_sbm` so SBM lifts its moist adiabat from a bottom layer coupled to Ts. Resolution-independent (fixed pressure depth). Gated by `use_surf_couple` in `exocol_rce_loop`.
  - **Non-deadlocking dry fallback**: SBM gates off when the column is net-subsaturated (`Σ(q−q_ref) ≤ 0`) — in a 1-D column at a long timestep that gate is a *trap* (the free troposphere dries, radiatively cools, and cannot recover → cold-dry collapse). The RCE loop therefore calls `convadj_dry` after `convadj_sbm`: it fires only on dry-superadiabatic layers (a **no-op on the moist-adiabatic reference**, so `'bulk'` Ts=288.02 is preserved bit-for-bit) and keeps the column convectively coupled when SBM is gated off, preventing the collapse. Gated by `use_dry_fallback`. This is the fix that made the prototype Frierson boundary layer unnecessary (the collapse was the gate, not the absence of BL mixing).
  - `'dry'`: potential-temperature stability criterion; q mixed fully (mass-weighted) when a pair is adjusted.
  - `'moist'`: rh-weighted local lapse rate `Γ_eff = rh·Γm + (1−rh)·Γd` where `rh = q/qsat` and `Γm = malr(T̄, p̄)` (saturated moist adiabat, phase-aware L). q-mixing is also rh-weighted: saturated pairs homogenize q; subsaturated pairs adjust only T.
  - `'manabe'`: fixed 6.5 K/km lapse rate (Manabe-Wetherald 1967); q mixed fully (mass-weighted).
  - `'zm'`: Zhang-McFarlane soft adjustment (Zhang & McFarlane 1995, consistent with ExoCAM). For each unstable pair, the full moist-adiabatic correction is computed but only a fraction `f = 1 − exp(−dt/τ_conv)` is applied. With the CAM default `τ_conv = 7200 s`, each step removes ~50% of the instability near equilibrium (dt ≈ 0.06 d) vs ~100% during fast transients (dt >> τ_conv). This smooths discrete convective temperature jumps and reduces TOA flux noise by ~5–10×. Namelist parameters `tau_conv` [s] (default 7200) and `cape_trigger` [J/kg] (default 0, disabled) are set in `&exocol_nml`. The RCE loop applies one soft pass then condensation then one hard cleanup pass (f=1) per outer step — no inner iteration for the soft pass. `compute_cape` (public from `exocol_convadj`) returns surface-parcel CAPE [J/kg] and is used internally for the trigger check.
  - `'sbm'`: simplified Betts-Miller (Frierson 2007), **recommended for moist Earth-like columns**. A whole-column scheme: relaxes T toward the moist adiabat lifted from the lowest level (parcel base = `tmid(pver)`), and — when `moisture_scheme='prognostic'` — q toward `rh_sbm·qsat(T_ref)` (`rh_sbm` default 0.7), over relaxation fraction `α = min(dt/τ_conv, 1)`. Unlike `moist`/`zm`, the **temperature target is the pure moist adiabat regardless of environmental RH** (RH enters only via the moisture target) — this fixes the too-steep lower-troposphere lapse rate and dry-aloft problem of the rh-weighted schemes (q no longer collapses to 0 above the boundary layer). Energy is conserved exactly by a single reference-temperature shift `dT_shift = (Σcp(T_ref−T)Δp/g − L·Σ(q−q_ref)Δp/g)/(cp·ΣΔp/g)`, so column enthalpy gain equals latent heat of precipitated water; water is conserved (`precip = α·Σ(q−q_ref)Δp/g`). `L = Lvap_T(ts)` matches the surface evap ledger. Cloud top = highest contiguous buoyant level (`T_ref ≥ T_env`); above it the column is left to radiative equilibrium. When the column is net subsaturated vs the reference (`Σ(q−q_ref) ≤ 0`) the SBM step itself makes no change, but the **non-deadlocking dry fallback** (`convadj_dry`, see above) then keeps the column convectively coupled so it cannot collapse cold-and-dry. At equilibrium `dt ≥ τ_conv` so `α→1` (hard adjustment to the moist adiabat — equivalent to konrad's `HardAdjustment`+`MoistLapseRate`); `α<1` during transients smooths TOA noise. Per outer step the RCE loop runs `convadj_surface` → `convadj_sbm` → `convadj_dry` → a `condense` (Newton satadj) mop-up for any residual stratiform supersaturation.

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

`ts` is prognostic and the slab budget includes turbulent fluxes (semi-implicit Planck damping):

```
ts ← ts + dt · (F_net_srf_rad − LE − SH) / (H_slab + dt·4σTs³)
H_slab = rho_w · cp_w · dz_slab = 1026 · 4000 · dz_slab   (dz_slab default 10 m → 4.1×10⁷ J/m²/K)
```

`F_net_srf_rad = (SWDN − SWUP) + (LWDN − LWUP)` at the surface interface. `LE` and `SH` are the `exocol_surface` fluxes (Monin-Obukhov `'mos'` by default), applied as bottom-layer sources for `tmid(pver)` and `h2ommr(pver)`; the latent enthalpy enters the column as condensation heat where the vapor later saturates. In addition, `convadj_surface` carries a **convective surface sensible flux** (it warms the sub-cloud layer toward Ts, debiting the slab by exactly the enthalpy gained — energy-conserving); that flux is folded into the `SH` diagnostic so the reported budget `F_srf_rad − LE − SH ≈ 0` closes.

After updating `tmid` and `ts`, `tint(pverp)` is pinned to `ts` in `update_tint`. The standalone convadj schemes (`dry`/`moist`/`manabe`/`zm`/`sbm`) operate purely on interior atmosphere-atmosphere pairs; the slab↔bottom coupling is done **only** by `convadj_surface`, which debits the slab for the enthalpy it adds to the bottom layer so column energy is conserved.

### Moisture and energy conservation

With `moisture_scheme = 'prognostic'` the column conserves moist static energy in steady state: column-integrated `F_TOA − F_net_srf_rad + LE + SH = 0`. Mass-balanced steady state requires surface evaporation rate = column precipitation rate (verified diagnostically by comparing `LE/L_v` to `precip_diag`).
