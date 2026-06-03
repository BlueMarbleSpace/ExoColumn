# Frierson Boundary-Layer Scheme — Handoff Note

**Status:** WORK IN PROGRESS, **PARKED behind namelist flags.** **Date:** 2026-06-02
(updated). The default config is the validated bulk scheme (Ts=288.01); the
Frierson BL is enabled by flags but currently produces an unphysical climate
(see "THE BLOCKER" below). Start here.

## Goal

Make ExoColumn's equilibrium climate **resolution-independent** (same Ts at any
`PVER`) via the **Frierson, Held & Zurita-Gotor (2006)** surface + boundary-layer
scheme. Root cause of the resolution dependence is the surface/BL coupling, not
the radiation timestep (see `project_resolution_root_cause`).

## What was done this session (2026-06-02)

1. **Prognostic, fixed-depth-anchored BL depth** (`exocol_pbl.F90`) — REPLACES
   the runaway/collapse-prone bulk-Richardson depth. The depth is carried as
   prognostic state, relaxed toward `min(bulk-Ri height, h_fixed)` floored at
   `h_floor`, where `h_fixed`/`h_floor`/`h_cap` are FIXED pressure depths
   (dp_mix=6000 Pa ≈ 500 m, dp_floor=2000 Pa, dp_cap=3e4 Pa) → same physical
   height at any resolution. **This part works: h is stable at ~510 m, SH/LE
   sensible, no runaway or collapse of the depth.** `pbl_reset()` clears it per
   run. (The original handoff's bulk-Ri-only depth is gone — it had no stable
   fixed point: deep runaway without a cap, LCL collapse with one.)

2. **Cold-trap energy-leak FIX** (`exocol_rce_loop.F90::apply_stratospheric_coldtrap`)
   — a REAL BUG, fixed and verified. The cold trap was CREATING water (assign/
   floor q to q_cp each step with no source); that spurious vapour precipitated
   and released **~9 W/m² of latent heat from nothing** (precip − evap ≈ +0.3
   mm/d → a structural TOA leak that blocked convergence and faked a slab-coupled
   oscillation). Now water-conservative: it sources its stratospheric water from
   the troposphere (Brewer-Dobson; vapour moved, not created). After the fix:
   precip≈evap (leak <0.01 W/m²), ⟨TOA⟩→0 monotonically, no oscillation. Added
   `BUDGET`/`WATER` per-window diagnostics to the RCE loop (kept — they caught
   this) and a `use_coldtrap` toggle. See `project_coldtrap_leak`.

## THE BLOCKER — the BL config collapses the troposphere (bistable)

With the BL enabled the column "converges" (Path B) but to an **unphysical
profile**: deep convection is confined to the bottom ~6 layers (the BL, moist),
and the **entire free troposphere goes cold and bone-dry** (q = cold-point value,
500 hPa at ~220 K, "tropopause" at ~947 hPa). Verified at every albedo (the
Ts=277/285.6/287.8 values are all collapsed states). The **bulk config never
does this** — its thick 400 m bottom layer keeps deep convection alive.

**Mechanism (verified):** SBM has a gate — if the convecting column is
net-subsaturated (Σ(q−q_ref) ≤ 0) it does nothing (can't sustain a precipitating
moist adiabat). During the thin-slab cold-start **warming overshoot**,
q_ref=rh_sbm·qsat(T_ref) rises faster than the BL delivers moisture → column goes
net-dry → gate trips → convection shuts off above the BL → the free troposphere
radiatively cools and dries → it **can't recover** (a cold layer can't hold
moisture without latent warming, and can't warm without deep convection, which
needs the moisture). A genuine **bistability**: warm-moist deep-convecting vs
cold-dry radiative attractors; the BL (cooling/drying the lowest layer + spreading
surface moisture over 500 m) plus the overshoot tip it into the dry one.

**Tried and REVERTED:** anchoring the SBM reference adiabat at `Ts` (surface
parcel) instead of `tmid(pver)` — konrad-consistent and arguably more correct,
but does NOT fix the collapse (the *gate*, not the parcel temperature, is the
blocker; a warmer adiabat demands MORE moisture → worse). `convadj_sbm` is back
to lifting from `tmid(pver)`.

## RECOMMENDED NEXT STEPS (need a design decision — touches convection/spin-up)

1. **Diagnostic first (lightest):** start the BL run from a warm-moist
   deep-convecting state (or suppress the SBM net-subsaturation gate during
   spin-up) and see if it STAYS warm-moist. If yes → it's a spin-up artifact
   (fix via initialization / gentler approach / no overshoot). If it collapses
   anyway → the BL fundamentally can't sustain deep convection here.
2. **Non-deadlocking convective closure:** a shallow / mass-flux branch that
   transports moist static energy up when the column is net-dry. NOTE a *pure
   moisture-redistribution* shallow branch does NOT work — the cold free trop just
   rains the moisture back out (it needs simultaneous latent warming). A real
   mass-flux scheme (or at least a combined T+q transport) is likely required.
   Big change to validated convection code.
3. **Reconsider the BL↔SBM coupling** (the BL spreading moisture over 500 m keeps
   the BL/column under q_ref; a shallower/gentler BL, or convection sourced
   before BL mixing, may help).

## How to run / toggle

- **Default (validated, Ts≈288.01):** `surface_flux='bulk'`, `pbl_scheme='none'`,
  `n_sfc_layers=0`, albedo 0.2674. This is what is committed.
- **Frierson BL (WIP, collapses):** `surface_flux='mos'`, `pbl_scheme='kprofile'`,
  `n_sfc_layers=12`. Watch the `WATER`/`BUDGET` lines (energy/water closure) and
  the profile via `python tools/plot_exocol.py` — the collapse shows as a cold,
  dry free troposphere with a ~947 hPa "tropopause".
- Build/run: `cd build && source /opt/intel/oneapi/setvars.sh && make PVER=70`,
  then from project root `./run/exocol.exe`. (NB: source setvars in the run shell
  too, or the exe fails on libimf.so.)

## Key files

- `src/exocol_pbl.F90` — K-profile BL + **prognostic fixed-depth-anchored depth**
  (works). Implicit tridiagonal diffusion with surface-flux bottom source.
- `src/exocol_surface.F90` — MOS + legacy bulk surface fluxes (`mos_drag_coef`).
- `src/exocol_rce_loop.F90` — surface block (no deposit when BL active), physics
  sub-stepping (`dt_phys_max`), BL call after SBM, **water-conservative cold trap**,
  `use_coldtrap` toggle, `BUDGET`/`WATER` diagnostics. Throwaway `diag_jacobian`
  flag + `diagnose_radiative_stiffness` still present (guarded `.false.`).
- `src/exocol_convadj.F90` — `convadj_sbm` (the net-subsaturation GATE that
  deadlocks is here: `if (Wvap <= 0._r8 ...) return`).

## Relevant memories

`project_bl_trop_collapse` (THE BLOCKER — read first), `project_coldtrap_leak`
(the leak fix), `project_frierson_bl_impl` (BL infrastructure + depth history),
`reference_frierson_bl` (exact formulas), `project_resolution_root_cause`,
`project_sbm`, `project_konrad_ref`.
