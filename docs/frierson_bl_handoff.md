# Frierson Boundary-Layer Scheme — RESOLVED (BL removed)

**Status:** CLOSED, 2026-06-03. The Frierson K-profile boundary layer, the hybrid
fine-surface grid, and the physics sub-stepping described in earlier versions of
this note have been **removed**. They are preserved in git history
(`exocol_pbl.F90` was deleted in the surface-scheme refactor commit).

## What happened

The BL was added to fix the model's resolution-dependent climate (Ts and surface
fluxes changed with vertical resolution). While bringing it up it produced an
unphysical **cold-dry troposphere collapse**, which was blamed on the BL.

Experiment 1 (2026-06-03) showed the collapse is reproduced with **no boundary
layer at all** — its real cause is the **SBM convective gate** (`if net-subsaturated:
return` in `convadj_sbm`) deadlocking in a 1-D column at a long timestep. The BL
was unnecessary *and* was masking that bug.

## The scheme that replaced it

Three focused mechanisms, all composition-general, on the plain log grid:

1. **MOS surface flux** (`exocol_surface`, `surface_flux='mos'`) with the drag
   coefficient referenced at the conventional **10 m** height (`z_flux_ref`) →
   resolution-independent and Earth-magnitude evaporation.
2. **Non-deadlocking dry-convective fallback** (`convadj_dry` after `convadj_sbm`,
   `use_dry_fallback`) → prevents the collapse; a no-op on the moist reference
   (bulk@70 = 288.02 preserved bit-for-bit).
3. **Slab-rooted surface mixed layer** (`convadj_surface`, `use_surf_couple`) →
   couples the lowest layers to Ts, restoring the moist-adiabatic profile; the
   slab is debited for the convective sensible flux (energy-conserving).

**Result (PVER=70):** Ts ≈ 288, LE/SH ≈ 68/16 (Bowen ≈ 0.24), P ≈ 2.4 mm/d,
energy and water budgets close, moist-adiabatic profile with realistic RH.
Resolution spread Ts(70→140) ≈ 2.2 K (vs 8.5 K for the old bulk scheme, which
additionally collapsed at high resolution).

## Open / future

- Residual ~2.2 K Ts resolution spread (the surface coupling still references the
  grid-dependent bottom layer for part of the flux). Accepted for now.
- Residual ~4 K cooler-aloft lapse rate vs the (warm-biased) bulk reference, from
  the dry-fallback/gate-margin interaction giving a slightly-steeper-than-moist
  lapse rate. Plausibly closer to konrad than bulk was; worth a konrad comparison.

See memory `project_mos_no_bl_result`, `project_bl_trop_collapse`,
`reference_frierson_bl`, and CLAUDE.md (`exocol_surface`, `exocol_convadj`).
