# Handoff: sub-273 K moist pseudoadiabat (panel-(d) hot-edge fix)

**Goal.** Close the remaining inner-HZ panel-(d) H₂O-profile residual at `Ts = 360–380 K`
by making the moist pseudoadiabat **below the 273.16 K triple point** carry the
water-vapour heat capacity (a full two-component treatment), instead of the
textbook `malr` fallback that uses `g/cp_dry` only.

**Success criterion.** In `tools/diag_lapse.py`, the lapse rate Γ = dlnT/dlnP below
273 K should track CLIMA (currently ~18 % too steep at 380 K); and the strat-H₂O
ratio (ExoColumn / CLIMA) at `Ts = 360/380 K` should move from **0.70 / 0.54**
toward ~1.0 (it is already ~1.0 at `Ts ≤ 340 K`). Verify with the snippet at the
bottom of this file.

---

## Diagnosis (verified 2026-06-11)

Lapse-rate comparison (`/tmp/diag_lapse.py`, reproduced below) at the pure-N₂,
liquid-cold-trap, BPS config of `hz_inner.py`:

| Ts=380 K | T | Γ_exo | Γ_CLIMA | ΔΓ |
|---|---|---|---|---|
| P = 10⁴ Pa | 306 K | 0.061 | 0.061 | +0.1 % |
| P = 3000 Pa | 285 K | 0.056 | 0.056 | +0.0 % |
| P = 1000 Pa | 269 K | 0.054 | 0.046 | **+18.7 %** |
| P = 100 Pa | 239 K | 0.049 | 0.042 | +17 % |

- **Above 273 K**: ExoColumn matches CLIMA to <0.1 % — our full two-component
  IAPWS-95 steam adiabat `exocol_steam::steam_dlnTdlnP_sat` (carries the heat
  capacity of *both* components).
- **Below 273 K**: the cold start falls back to `exocol_convadj::malr`
  (`exocol_coldstart.F90`: `if (use_nonideal_adiabat .and. T_lev >= IAPWS_TT
  .and. T_lev < IAPWS_TC) <steam> else <malr>`). `malr` is the textbook ideal
  moist adiabat Γ_d = g/cp_**dry**, which **drops the water-vapour heat
  capacity**. At the hot-case upper-troposphere mixing ratios (ws ~ 0.3) that
  term is large, so the adiabat is too steep → cools to the 200 K cap at higher
  pressure → freeze-dries to a drier strat. The error grows with Ts (more upper-
  trop vapour), matching the observed Ts-dependence (340 K ratio 1.02 ✓, 380 K 0.54).

The divergence **onset is exactly at the triple point** — the 273 K steam→malr
switch is the artifact.

## Failed approach — do NOT repeat

Making `malr` saturate over supercooled **liquid** (`esat_cc_liq` + L_v) below
273 K when `cold_trap_phase='liquid'` made it **worse** (strat VMR 340: 1.02→0.78,
360: 0.70→0.42, 380: 0.54→0.32; Γ +34 % vs +18 %). At high ws, more vapour +
lower L_vap net to a *steeper* textbook adiabat. **Phase (ice/liquid) is a
secondary, wrong-way lever; the water-vapour heat capacity is the real one.**
This change was implemented and reverted (commit history); `malr` is back to its
original phase-aware-ice form.

## Physical constraint (read before coding)

Below ~235 K there is **no supercooled liquid–vapour saturation** — real water is
ice. So `exocol_iapws95::iapws95_sat` (the Maxwell construction) is **not valid at
the 200 K cold point**; you cannot simply lower the `IAPWS_TT` bound and call
`steam_dlnTdlnP_sat` down to 200 K (the solver will fail or return garbage). CLIMA
sidesteps this with a *parameterized* supercooled-liquid saturation extrapolation,
not a real EOS.

## Candidate approaches

1. **Two-component ice pseudoadiabat (recommended).** Derive/implement a moist
   adiabat below 273 K that carries `cp_v` (and the H₂O latent heat of
   sublimation) over **ice** saturation — i.e. the ice analogue of
   `steam_dlnTdlnP_sat`'s dilute limit. This keeps the physically-correct ice
   phase (arguably better than CLIMA) while fixing the dropped heat capacity.
2. **Extend the steam adiabat with a supercooled-liquid parameterization**
   (CLIMA-style) down to the cold point — matches CLIMA but adopts the less
   physical supercooled-liquid convention; needs a metastable `es`/`L` fit since
   IAPWS-95 saturation is unavailable there.
3. **Minimal patch**: add the vapour heat-capacity term to `malr` (e.g. effective
   cp = cp_dry + ws·cp_v with the correct moist-adiabat algebra). Cheapest, but
   re-derive carefully — the naive cp_dry→cp_total substitution is *not* the
   correct formula, and the first attempt above shows reasoning errors are easy.

Whichever path: keep the default (current ice `malr`) bit-identical, and gate any
new behaviour so the Earth RCE and existing reference cases are unchanged.

## Relevant code

- `src/exocol_coldstart.F90` — the steam/`malr` switch at lines ~314/383/436
  (`if … T_lev >= IAPWS_TT … else Gm = malr(...)`). The below-273 K branch is what
  to replace.
- `src/exocol_steam.F90` — `steam_dlnTdlnP_sat` (the full two-component adiabat to
  emulate below freezing); header documents the A4/A5 Kasting-88 formulation.
- `src/exocol_convadj.F90` — `malr` (textbook, pure), `esat_cc`/`esat_cc_liq`
  (ice / supercooled-liquid CC saturation), `Lvap_T` (phase-aware latent heat).
- `src/exocol_iapws95.F90` — `iapws95_sat` (NOT valid below the triple point),
  `iapws95_psat_aux` (WP analytic Psat, *is* evaluable below 273 K but liquid-branch).

## Diagnostics

- `tools/diag_ptz_clima.py` — overlays ExoColumn vs CLIMA T(P), H₂O VMR(P), z(P)
  at Ts = 340/360/380 K → `tools/diag_ptz_clima.png`. The localizer.
- `/tmp/diag_lapse.py` (reproduce: it's a ~40-line script; the Γ table above) —
  prints Γ_exo vs Γ_CLIMA at matched pressures. The decisive test.
- After any change: rebuild `make PVER=200`, run the two diagnostics, confirm the
  Earth RCE is still `Ts = 287.845 K / 2107 steps` (bit-identical), then regenerate
  the figure (`python reference/moist_runaway/hz_inner.py`, ~28 min dual sweep).

## Context

This is the last open item from the 2026-06 Kopparapu-(2013) apples-to-apples
audit (memory `project_bps_continuum`, `project_ihz_kasting_fidelity`,
`project_hz_roadmap`). Everything else is matched: composition (pure N₂ + 330 ppm
CO₂), albedo 0.32, S0 = 1360 G2V, isothermal 200 K stratosphere, BPS continuum,
liquid cold trap. The HZ *limits* (≤ 340 K, incl. the moist-greenhouse boundary)
already match CLIMA — this fix is for the 360–380 K profiles only, beyond the HZ
edge, so it is a fidelity nicety, not a result-changing correction.
