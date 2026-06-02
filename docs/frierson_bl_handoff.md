# Frierson Boundary-Layer Scheme — Handoff Note

**Status:** WORK IN PROGRESS. **Date:** 2026-06-02. Start here.

## Goal

Make ExoColumn's equilibrium climate **resolution-independent** (same Ts at any
`PVER`). The root cause was diagnosed earlier (`docs/resolution_independence.md`
is now partly obsolete): it is **not** the radiation timestep — it is the
**surface/boundary-layer coupling**. Without BL mixing the near-surface T,q
(which set the bulk surface fluxes) depend on the bottom grid-layer thickness,
so refining the grid changes LE/SH (a ~4× swing in SH between PVER=70 and 140)
and hence the climate (PVER=140 reference froze at Ts=279.5, not 288).

The fix being built is the **Frierson, Held & Zurita-Gotor (2006)** surface +
boundary-layer scheme (the lineage of ExoColumn's SBM convection). Exact
formulas are in the `reference_frierson_bl` memory; paper at
https://www.gfdl.noaa.gov/bibliography/related_files/dmwf0601.pdf.

## Three pieces (all implemented; toggled by namelist)

1. **MOS surface fluxes** (`exocol_surface.F90`, `surface_flux='mos'`):
   potential-temperature / dry-static-energy flux differences with a
   height-dependent drag `C = κ²/ln²(z_a/z0)` (stable branch reduces it). The
   `z_a` dependence makes the flux resolution-independent in the surface layer.
   Legacy `'bulk'` (fixed C_D, actual-T) retained and is the default.
2. **Hybrid fine-surface grid** (`n_sfc_layers=12`): puts the lowest level at
   ~8.5 m, **inside** the surface layer where MOS is valid (the default 400 m
   log-grid bottom is too high — that is *why* a gradient BL scheme could not
   engage at coarse resolution).
3. **Diffusive K-profile BL** (`exocol_pbl.F90`, `pbl_scheme='kprofile'`):
   mechanical eddy diffusivity `K = κ·U·√C·z` (surface layer) + parabolic taper
   to zero at the BL top `h`; `h` from a bulk-Richardson criterion (Ric=1) on
   virtual dry static energy. Conservative implicit (backward-Euler) tridiagonal
   diffusion of dry static energy `s = cp·T + g·z` and `q`.

### Two architectural fixes already made (both essential, both in)

- **Surface flux injected as a BL bottom *source*, not deposited into the bottom
  layer.** `exocol_rce_loop.F90` no longer deposits LE/SH into `tmid(pver)`,
  `h2ommr(pver)` (when a BL scheme is active); instead `surf_sh` and
  `surf_e=LE/L` are passed to `pbl_diffuse` and added to the bottom-layer RHS of
  the implicit solve (`tridiag_diffuse` `src_bot` arg), so the flux is spread
  through the BL in one conservative implicit step — no thin-layer spike.
  Energy-conserving: slab debits LE+SH, BL injects the same. The legacy explicit
  deposit is kept only for `pbl_scheme='none'`.
- **Physics sub-stepping.** Frierson's BL+SBM assume short (~minute) steps;
  ExoColumn's radiative step is 0.1–0.4 d. The local physics (surface, BL, SBM,
  condense, cold-trap) is now sub-stepped at `dt_sub ≤ dt_phys_max = 0.01 d`
  (~14 min) inside each frozen-radiation outer step. This reuses the existing
  `isub_rad`/`N_sub` loop (which already runs all local physics per sub-step):
  `N_sub = max(radiative N_sub, ceil(dt_outer/dt_phys_max))`, gated to
  `pbl_scheme/='none'`. This **fixed the timestepping instability** — gave the
  best result so far: a stable, correct BL (h≈454 m, SH≈15, LE≈63, Ts 288→289,
  TOA falling) for ~110 model days at PVER=70.

## THE REMAINING PROBLEM — BL depth diagnosis is metastable

With physics sub-stepping the *timestep* is no longer the issue, but the
**bulk-Richardson BL depth itself does not have a stable fixed point at the
right (~450 m) depth in single-column RCE**:
- From the good 454 m state it slowly drifts (over thousands of 0.008 d
  sub-steps, so NOT a timestep artifact). **Without** any cap: `h` deepens to
  ~1.7–2.7 km and SH collapses to ~0.6 (deep mixing equilibrates Θ_BL→Ts,
  killing the surface sensible gradient; LE stays ~83, near reference). **With**
  an LCL cap (tried, then removed): q_bot slowly saturates → LCL collapses → BL
  switches off → thin bottom layer goes radiatively stiff (max|HR|→124) →
  diverges.
- Why: Frierson diagnoses `h` from bulk-Ri **every step** and the paper notes
  the BL *"extends to the tropopause on occasion"* — i.e. it is noisy in his GCM
  too, but **GCM variability + horizontal averaging keep the time-mean
  sensible.** A single-column RCE has no such variability, so the depth gets
  *stuck* in the metastable deep state.

## RECOMMENDED NEXT STEP — prognostic, entrainment-limited BL depth

Make `h` **prognostic**: store it, and each step relax it toward the bulk-Ri
target at a **bounded entrainment rate** (e.g. `dh/dt` capped, or a relaxation
`h += (h_target − h)·min(dt/τ_h, 1)`) instead of snapping to a fresh, jumpy
bulk-Ri value. This damps the drift/oscillation and should hold `h` at the
stable ~450 m. The 454 m metastable state IS the physically-correct answer
(SH=15, LE=63) — the depth scheme just needs to STAY there.

`h` would need to live in module state (`exocol_mod` or a `save` in
`exocol_pbl`), initialised on first call, and reset when the loop re-enters.
Alternatives if (1) proves fiddly: a physical depth cap (~ps−150 hPa); a smaller
Ric; or fixed-pressure-depth mixing (crude but robust).

## How to run / toggle

- **Reference (default, Ts≈288.01):** namelist `surface_flux='bulk'`,
  `pbl_scheme='none'`, `n_sfc_layers=0`. This is what is committed as the
  default and what the regression must keep reproducing.
- **Frierson scheme (WIP):** set `surface_flux='mos'`, `pbl_scheme='kprofile'`,
  `n_sfc_layers=12` in `exocol_config.nml`.
- Build/run: `cd build && source /opt/intel/oneapi/setvars.sh && make PVER=70`,
  then from project root `./run/exocol.exe`. Sweep `make clean && make PVER=N`.
- Watch stdout: the loop prints `PBL: h[m]= … ktop= …` each report. Healthy =
  `h` steady ~400–600 m, SH ~15, LE ~60–90, max|HR| modest, TOA→0. Failure =
  `h` drifting to km-scale with SH→0, or `h`→bottom with max|HR| exploding.

## Test plan once depth is stable

1. PVER=70 with the Frierson scheme converges to a sensible Ts (it will NOT be
   288.01 — adding a real BL changes the climate; **re-tune albedo** to recover
   Ts≈288 and re-validate vs konrad).
2. Sweep PVER=70/100/140/200 (hybrid grid) → confirm the converged Ts is
   resolution-INDEPENDENT (the whole point). This is the success criterion.

## Key files

- `src/exocol_surface.F90` — MOS + legacy bulk surface fluxes (`mos_drag_coef`).
- `src/exocol_pbl.F90` — K-profile BL, implicit tridiagonal diffusion with
  surface-source bottom BC. **The depth diagnosis to fix is here** (bulk-Ri
  loop computing `ktop`).
- `src/exocol_rce_loop.F90` — surface block (no deposit), physics sub-stepping
  (`dt_phys_max`, N_sub), BL call after SBM. Throwaway `diag_jacobian` flag +
  `diagnose_radiative_stiffness` still present (guarded `.false.`; can delete).
- `src/exocol_config.F90` / `exocol_config.nml` — `surface_flux`, `pbl_scheme`,
  `z0_rough`, `n_sfc_layers`.

## Relevant memories

`project_frierson_bl_impl` (detailed status), `reference_frierson_bl` (exact
formulas), `project_resolution_root_cause` (why it's physical not numerical),
`project_resolution_dt_stiffness` (superseded), `project_sbm`,
`project_stratospheric_coldtrap`.
