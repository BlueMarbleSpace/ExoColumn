# Resolution-Independent Convergence — Handoff Note

**Status:** OPEN — required work. **Date:** 2026-06-01.

ExoColumn must converge for an **arbitrary number of vertical levels** (`PVER`). Today it
converges cleanly only at `PVER=70`; at `PVER≥100` it fails to converge within `nmax`. This
note documents the root cause, two failed fixes (so they are not repeated), and the
recommended path. Start here.

---

## 1. The problem

`PVER=70` converges: Ts=288.01 K, cold-point 199.27 K, OLR=258.75 W/m², 2312 steps / 1001
model days. `PVER=100` does **not** converge within `nmax=200000`: a thin layer develops a
persistent radiative heating `|LWHR+SWHR| ≈ 16 K/day` (vs ≈1.8 at 70 levels) that convection
cancels almost exactly — the profile freezes (net ΔT≈0) but the slab cannot equilibrate, Ts
crawls to a non-equilibrium ~280 K, TOA stalls ~0.24 W/m². **The Ts at 100 levels is a stalled
crawl, not a real equilibrium — do not trust it.**

### Root cause (the binding constraint)

`src/exocol_rce_loop.F90`, main loop:
- line ~260: `max_hr = maxval(abs(LWHR + SWHR))` — radiative heating rate.
- line ~269: `dt_days = cfl_safety * dT_target / max_hr` — adaptive timestep.
- line ~279: `tmid(k) = tmid(k) + dt_days*(LWHR(k)+SWHR(k))` — **explicit** radiation, then
  convection/condensation/surface adjust.

The explicit radiation substep requires `dt·max|radiative HR| ≤ dT_target` for its own
stability. At fine resolution a thin layer (sharp flux divergence — suspect the H₂O cold-trap
moisture step at the cold point / convective top, smeared at 70 levels, spiking when thinned)
has large radiative HR. That clamps `dt`, even though the layer is convectively balanced and
contributes ~zero *net* tendency. With `dt` clamped, the slab (10 m, large heat capacity) needs
more than `nmax` steps to equilibrate. **The grid is not the limiter — the explicit-radiation
CFL is.** Note: convergence in a convecting column is via **Path B** (profile stability +
`|⟨TOA⟩|<0.1`), never Path A (`max|radiative HR|<0.01`) — radiation is balanced by convection,
not zero, so Path A is the wrong test for these columns.

---

## 2. What was tried and FAILED (do not repeat as-is)

### (a) Key `dt` on the NET per-step ΔT — FAILED (instability)
Sized `dt` by the realized `max|Δtmid_net|/dt` from the previous step (lagged, growth capped
1.5×/step), keeping explicit radiation. **Broke even PVER=70.** Once `dt` grew to `dt_max=1 d`
on quiet steps, the explicit radiation substep `tmid += dt·HR` overshot a stiff layer (~17 K in
one step) before convection/surface could react; SH exploded 36→232 W/m²; the column derailed
to a spurious cold attractor (Ts=279.6) and falsely "converged" there via Path B.
**Lesson: keying `dt` on net tendency is incompatible with applying radiation in one EXPLICIT
step. The lagged limiter reacts one step too late to the overshoot.**

### (b) Diagonal semi-implicit radiation — FAILED (energy non-conservation)
Linearized backward-Euler. `κ(k)=∂HR/∂T` from one extra radiation call (perturb `tmid` by
+1 K uniformly, restore; clamp `κ≤0`). Apply `ΔT_rad = dt·HR/(1−dt·κ)`; size
`dt = cfl·dT_target / max_k(0, |HR|−dT_target·|κ|)`. The `dt` math is correct (bounds applied
`|ΔT_rad|≤dT_target` exactly; a layer with saturated displacement `|HR|/|κ|≤dT_target` imposes
no limit) and it did **not** blow up. **But it leaks energy → slow secular drift.** At PVER=70
it passed through the correct equilibrium (Ts=288.06, TOA=0.043 at step 308) then drifted
monotonically to Ts=289.84 / TOA=2.9 and stalled. The damped increment `dt·HR/(1+dt|κ|)`
applies *less* cooling than the true flux divergence `dt·HR`, while OLR/surface fluxes report
the full radiative loss → the atmosphere under-cools → energy accumulates.
**Lesson: a DIAGONAL semi-implicit is inherently non-conserving — it drops the off-diagonal
flux coupling, and that mismatch is the leak. Cannot be patched into conservation.**

---

## 3. Recommended path: radiation sub-cycling (frozen HR)

The two failures point straight at the fix. The net-tendency `dt` idea (a) was right; it only
failed because radiation was applied in **one** explicit shot at large `dt`. Apply the radiation
tendency in **N explicit sub-steps with convective adjustment between them**, keeping HR frozen
over the outer step:

```
Outer step (one expensive radiation call):
  call exocol_rad_tend → HR(k), fluxes (OLR, surface)        ! once
  choose dt_outer  (see below)
  N_sub = ceil( dt_outer * max|HR| / dT_target )             ! explicit-CFL sub-steps
  dt_sub = dt_outer / N_sub
  do s = 1, N_sub
    tmid(k) += dt_sub * HR(k)         ! frozen HR, small stable explicit step
    surface fluxes + slab + bottom-layer deposit  (existing dt_sfc subcycle)
    convadj + condense + cold trap
  end do
```

- **Energy-conserving exactly:** the full `dt_outer·HR` is applied; convection/condensation/
  surface all conserve. No leak (unlike (b)).
- **Stable:** each sub-step satisfies the explicit CFL (unlike (a)).
- **`dt_outer` keyed to the SLOW (net) tendency** — now SAFE because the stiff layer is
  handled by the sub-steps, not one big explicit jump. Use the net-tendency limiter from (a),
  or limit `dt_outer` by the slab/surface evolution rate.
- **Cost win:** radiation (68-band correlated-k) is the expensive call and happens **once** per
  `dt_outer`; the `N_sub` sub-steps are cheap (convection only). At high res `max|HR|~16`,
  `dt_outer~1 d` ⇒ `N_sub~16` cheap steps per radiation call — far fewer radiation calls per
  model-day than today's clamped `dt~0.05 d` (one radiation call every 0.05 d).

### Caveats to validate
- **Frozen-HR accuracy:** HR is held constant while `tmid` evolves over `dt_outer`. Exact at the
  stiff layer (it returns to the convective reference each sub-step), but slowly-evolving free
  layers (stratosphere) drift under frozen HR — bound `dt_outer` so this O(`dt_outer·dHR/dt`)
  error is small. This is the standard "infrequent radiation" approximation.
- **Integration with the existing surface/slab subcycle** (`dt_sfc ≤ τ_surf`, already in the
  loop) — the radiation sub-cycle wraps it; check the nesting.
- **TOA/flux diagnostics** are from the start-of-outer-step radiation call; fine for the
  budget since energy is conserved, but confirm the window-mean TOA used by Path B is sensible.

### Alternative (heavier, if sub-cycling accuracy proves insufficient)
Full/banded-Jacobian implicit radiation: solve `(I − dt·J)ΔT = dt·HR` with `J=∂HR/∂T`. Correct
and unconditionally stable, but `J` needs ≈`pver` radiation calls per step (or a banded
approximation) — likely too expensive. Sub-cycling is the better first attempt.

---

## 4. How to test

1. **Regression at PVER=70 first.** Any new scheme must reproduce Ts=288.01 K, cold-point
   199.27 K, OLR=258.75 W/m² (current explicit baseline). If 70 changes, the scheme is wrong.
2. **Then sweep PVER=70/100/140/200** (uniform log grid). Build: `make clean && make PVER=N`
   (compile-time; `src/exoplanet_mod.F90` is regenerated each build — see CLAUDE.md). Run from
   project root. Each run writes `iofiles/exocol_out.nc`; copy to `iofiles/sweep_N.nc`.
3. **Watch live convergence** via stdout — the loop now has `flush(6)` after the `step=` and
   `WINDOW@t=` prints, so a redirected run (`stdbuf -oL run/exocol.exe > log 2>&1`) streams
   progress. Healthy: `max|HR|` and `TOA` fall, `dt` rises, Ts steadies. Failure signatures
   seen here: `max|HR|` pinned ~16–23 while Ts drifts (stall), or Ts derailing to ~280/~290 with
   SH exploding (instability/leak).
4. **Convergence criteria** (`exocol_rce_loop.F90`): Path B = window-mean profile stable +
   `|⟨TOA⟩|<0.1 W/m²`. Confirm convergence is in a **reasonable step count** (70-level reference:
   2312 steps), Ts is resolution-converged (compare 70/100/140/200), and energy is conserved
   (no secular Ts drift — the (b) failure mode).

### Why this matters beyond convergence
Once arbitrary `PVER` converges, the original motivating question can finally be answered:
does finer resolution (a) round the single-layer **SBM convective-top kink at ~200 hPa** (which
is a real scheme feature — see `project_sbm_coldnotch` memory — not a bug), and (b) reduce the
~2.6 K cold-point warm bias vs konrad (currently charged to ExoRT-vs-RRTMG radiation, but partly
plausibly numerical diffusion)? Neither is known yet — the high-res runs never reached a
trustworthy equilibrium.

---

## 5. Repo state at handoff
- Working tree: only `flush(6)` added to `src/exocol_rce_loop.F90` (after the `step=` and
  `WINDOW@t=` prints). **Keep it** — it is what made all three failures diagnosable in minutes.
- `run/exocol.exe` and `iofiles/exocol_out.nc`: clean PVER=70 baseline (Ts=288.01).
- All failed-fix code has been reverted. Nothing committed.
- Relevant memories: `project_resolution_dt_stiffness`, `project_sbm_coldnotch`,
  `project_convergence_speed`, `project_stratospheric_coldtrap`, `project_sbm`.
