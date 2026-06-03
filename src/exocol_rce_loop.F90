module exocol_rce_loop
! RCE time-marching loop for ExoColumn.
!
! run_rce_loop time-steps the column forward with a CFL-limited adaptive
! timestep and the following physics per step, in order:
!
!   1. Radiation tendency LWHR, SWHR on tmid (ExoRT aerad_driver).
!   2. Adaptive timestep: dt chosen so that |dt · max|HR|| ≤ dT_target.
!        dt = clamp( cfl_safety · dT_target / max|HR| , dt_min , dt_max )
!      Standard radiative-CFL criterion used in research-grade RCE codes
!      (Frierson, Manabe-Wetherald): keep the per-step temperature change
!      in the most active layer below ~0.5 K so explicit-Euler radiation
!      stays in its linear regime.
!   3. Apply radiation tendency: tmid ← tmid + dt · (LWHR + SWHR).
!   4–6. Subcycled surface turbulent exchange:
!      τ_surf = (pdel/g)/(ρ·C_D·U) — bottom-layer moisture relaxation time.
!      The old implicit-Euler damping 1/(1+dt/τ) collapsed LE → 0 for thin
!      layers (τ << dt); subcycling at dt_sfc < τ keeps explicit Euler stable
!      for any bottom-layer thickness.  n_sub = max(1, floor(dt/τ)+1) so that
!      for thick layers (τ >> dt) n_sub = 1 and cost is identical to before.
!      Each subcycle updates ts, tmid(pver), and h2ommr(pver) together.
!   7. Saturation adjustment — convadj + satadj (Newton) + convadj.
!      Pass A: convadj equilibrates the lapse rate.
!      Pass B: satadj removes all supersaturation in one Newton step, releasing
!              Lvap_T(ts)·Δq into tmid; q ≤ qsat(T) is guaranteed after one call.
!      Pass C: convadj cleans up lapse-rate instability from the latent release.
!   8. Update derived (pdeldry, pintdry), tint, zint.
!
! CFL / stability summary (Earth-like column at dt_max = 1 d):
!   Radiative ΔT/step          → dt_max · max|HR|        bounded by dT_target.
!   Surface turbulent flux     → τ_surf ≈ 8 h            implicit damping.
!   Slab radiative equilibration → τ_slab ≈ 350 d        well below dt_max.
!   Condensation                → saturation adjustment  one Newton step/layer.
!
! Time accounting: model_time_days accumulates dt_days each step.  Stability
! snapshots and console prints are model-time based (not step-count based)
! so output remains readable when dt shrinks during transients.
!
! Moisture is fully prognostic by default.  Setting moisture_scheme='fixed_rh'
! restores the legacy CC-relaxation closure; moisture_scheme='off' holds q
! frozen at the input values.  When the prognostic scheme is active there is
! no fixed-RH or rh_init state — q evolves via evaporation source, convective
! transport, and condensation sink.
!
! Convergence is declared when either path is satisfied:
!
!   Path A — radiative equilibrium (quiescent atmosphere):
!     max |LWHR(k) + SWHR(k)| < hr_tol   AND   |TOA net flux| < toa_tol
!
!   Path B — profile stability (convectively active atmosphere):
!     max |Δtmid| < prof_stab_tol   AND  |ΔTs| < ts_stab_tol  AND
!     max |Δh2ommr| < q_stab_tol    AND  |⟨TOA⟩| < toa_tol
!
!   Path B handles cases where convective adjustment continuously balances a
!   large instantaneous radiative tendency.  The h2ommr term guards against
!   declaring victory while T/Ts appear frozen but moisture is still drifting
!   (which would leave a residual TOA imbalance from the latent heat sink).

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use shr_const_mod,   only: SHR_CONST_CSEC, SHR_CONST_RGAS, SHR_CONST_MWWV, &
                             SHR_CONST_STEBOL, SHR_CONST_RWV
  use ppgrid,          only: pver, pverp
  use exocol_mod
  use exocol_radiation, only: exocol_rad_tend
  use exocol_config,    only: conv_scheme, moisture_scheme, wind_speed, C_D, &
                              cfg_dz_slab => dz_slab, tau_conv, cape_trigger, &
                              rh_sbm, surface_flux, z0_rough
  use exocol_convadj,   only: convadj_dry, convadj_surface, convadj_moist, &
                              convadj_manabe, convadj_zm, convadj_sbm, &
                              esat_cc, Lvap_T
  use exocol_surface,   only: compute_surface_fluxes

  implicit none
  private

  public :: run_rce_loop

  ! ---- Adaptive timestep parameters ----
  ! Each step, dt is chosen so the radiation-induced ΔT in the most actively
  ! heated/cooled layer stays below dT_target.  dt is clamped to [dt_min, dt_max].
  ! Near radiative equilibrium dt → dt_max and the loop reverts to the classic
  ! "1 day per step" behaviour; during cold starts where max|HR| ~ 10 K/day, dt
  ! shrinks to ~0.04 day so the explicit Euler radiation step remains in its
  ! linear regime (Manabe-Wetherald 1967; CliMT; RRTM-RCE).
  real(r8), parameter :: dt_max      = 1.0_r8       ! cap [days]
  real(r8), parameter :: dt_min      = 1.0e-4_r8    ! floor [days] (~10 s)
  real(r8), parameter :: dT_target   = 1.0_r8       ! target |dt·HR| per sub-step [K]
  real(r8), parameter :: cfl_safety  = 0.8_r8       ! safety factor [-]

  ! ---- Radiation sub-cycling (resolution-independent convergence) ----
  ! The expensive 68-band radiation call is made ONCE per outer step; the
  ! resulting heating rate HR(k) is held frozen and applied in N_sub cheap
  ! explicit sub-steps, with convection / condensation / surface exchange
  ! between each.  This decouples two limits the old single-explicit step
  ! conflated:
  !   * the explicit radiative CFL (dt_sub·max|HR| ≤ dT_target) — a STABILITY
  !     limit set by the single stiffest layer.  At fine resolution a thin
  !     cold-point / cold-trap layer concentrates a sharp flux divergence into
  !     little mass and develops max|HR| ~ 16 K/day; keying the whole timestep
  !     on it throttles the column to ~0.05 d and the slab never equilibrates
  !     (the original resolution-stiffness bug).  Handled by N_sub.
  !   * the explicit radiative CFL of the BULK column — the rate at which the
  !     climatically-relevant (mass-bearing) layers evolve.  This sets dt_outer.
  ! dt_outer is keyed to hr_for_dt = max|HR| EXCLUDING the stiffest layers that
  ! together hold less than f_excl_mass of the column mass (the thin cold-point
  ! spikes).  At PVER=70 there is no stiff outlier so hr_for_dt = max|HR| and
  ! dt_outer reproduces the reference scheme; at high resolution the cold-point
  ! spike is excluded from dt_outer (it cannot drive the slow climate) and is
  ! absorbed by N_sub = ceil(dt_outer·max|HR|/(cfl·dT_target)).  Applying the
  ! full dt_outer·HR in N_sub conserving sub-steps is energy-conserving exactly.
  ! (An earlier attempt keyed dt_outer on the NET tendency; that runs away —
  ! convection cancels radiation so net ΔT stays small while the column is still
  ! radiatively stiff, and the slab derails to a spurious cold attractor.)
  ! See docs/resolution_independence.md.  With f_excl_mass = 0 (default) no layer
  ! is excluded, so dt_outer = the bulk radiative CFL and N_sub = 1 — i.e. the
  ! reference single-explicit-step scheme; the scaffolding is dormant.
  real(r8), parameter :: f_excl_mass = 0.0_r8       ! column-mass fraction of stiffest layers excluded from dt_outer
  integer,  parameter :: nsub_max    = 2000         ! cap on radiation sub-steps

  ! Throwaway diagnostic: when .true., the loop dumps the stiff-layer profile
  ! and probes the radiative Jacobian (∂HR/∂T) once the stiff layer forms,
  ! then stops.  Used to design the implicit-radiation solver.  Set .false.
  ! for production.
  logical, parameter :: diag_jacobian = .false.

  ! Cold-trap toggle (default on).  apply_stratospheric_coldtrap is now
  ! water-conservative (it sources its stratospheric water from the troposphere
  ! rather than creating it), so it no longer leaks latent heat at TOA.
  logical, parameter :: use_coldtrap = .true.

  ! Non-deadlocking convective fallback (default on).  The SBM scheme gates OFF
  ! when the convecting column is net-subsaturated relative to its reference
  ! (Wvap <= 0): it cannot sustain a precipitating moist adiabat, so it makes no
  ! change.  In a 1-D column at a long timestep that gate is a trap — once the
  ! free troposphere dries (e.g. during the cold-start warming overshoot) SBM
  ! stays off, the free troposphere radiatively cools and dries further, and it
  ! cannot recover (a cold layer holds no moisture without latent warming, and
  ! cannot warm without deep convection, which needs the moisture).  The column
  ! collapses to an unphysical cold-dry state with a near-surface inversion.
  ! This collapse occurs independently of the boundary layer (it appears with
  ! the bulk and MOS surface schemes, with and without BL mixing); it is the
  ! gate, not the surface coupling.  Fix: when SBM is gated off, fall back to
  ! DRY convective adjustment, which mixes any dry-unstable layers to the dry
  ! adiabat (conserving enthalpy and water — no latent-heat or moisture
  ! creation).  This keeps the column convectively coupled and the sub-cloud
  ! layer warm, so surface evaporation can re-moisten it until the moist branch
  ! re-engages.  convadj_dry is a no-op on a moist-adiabatic or radiatively
  ! stable profile (it only fires on dry-superadiabatic layers), so the
  ! supersaturated 'bulk' reference equilibrium (Ts=288) is preserved.
  logical, parameter :: use_dry_fallback = .true.

  ! Surface-coupled mixed layer (default on).  convadj_surface roots a dry
  ! convective adjustment at the slab so the lowest layers are mixed toward the
  ! surface temperature (a fixed-pressure-depth sub-cloud mixed layer), the
  ! resolution-independent replacement for the Frierson K-profile boundary layer.
  ! Without it, a bulk/MOS surface flux at a ~400 m lowest level leaves the bottom
  ! layer radiatively decoupled (super-adiabatic gap) — a too-cold convective
  ! parcel base (cold free troposphere) and a saturating bottom layer that chokes
  ! evaporation (dry column → moist convection gates off).  It is a no-op on a
  ! stable surface layer (theta_bottom >= theta_surf, as in the 'bulk' reference),
  ! so it does not disturb that equilibrium.  Energy-conserving: the bottom-layer
  ! warming is debited from the slab.
  logical, parameter :: use_surf_couple = .true.

  ! Current adaptive timestep — set each step in run_rce_loop, read by
  ! condense() (which converts cond_heating to K/day for diagnostics).
  real(r8), save      :: dt_days = dt_max

  integer,  parameter :: nmax     = 200000          ! maximum iterations
  real(r8), parameter :: hr_tol   = 0.01_r8         ! heating-rate convergence [K/day]
  real(r8), parameter :: toa_tol  = 0.1_r8          ! TOA flux convergence [W/m²]

  ! Console diagnostics: print every ~10 model days regardless of how many
  ! steps that took (adaptive dt means step count varies).
  real(r8), parameter :: print_every_days = 10._r8

  ! Profile-stability convergence (Path B).  Snapshots are taken every
  ! stab_check_days of model time and compared to the previous snapshot.
  !
  ! prof_stab_tol = 1.0 K: 100-day window means of tmid have an irreducible
  !   noise floor of ~0.3–0.6 K driven by discrete convective bursts (τ_corr
  !   ~3–5 days, N_eff ~20–25 per window → σ ≈ 0.3 K).  The maximum observed
  !   at true equilibrium is 0.59 K.  ΔTs_win is the primary safeguard; ΔTpro
  !   is kept as a redundant guard against 1+ K profile drift.
  ! q_stab_tol = 5e-4: discrete condensation bursts create LE flicker that
  !   prevents 100-day mean h2ommr from settling below 1e-4 kg/kg.
  ! TOA condition: only |⟨TOA⟩| < toa_tol is used (no ΔTOA branch).  The
  !   ΔTOA alternative allowed false convergence during transient slow-drift
  !   phases where ΔTOA happened to be small while ⟨TOA⟩ was still 0.5–1 W/m².
  real(r8), parameter :: stab_check_days = 100._r8   ! snapshot interval [days]
  real(r8), parameter :: prof_stab_tol = 1.0_r8      ! max tmid change [K]
  real(r8), parameter :: ts_stab_tol   = 0.02_r8     ! max Ts   change [K]
  real(r8), parameter :: q_stab_tol    = 5.0e-4_r8   ! max h2ommr change [kg/kg]

  ! Convadj fixed-point iteration (within each outer step).  Each convadj
  ! pass equilibrates the lapse rate; iteration exits when no layer-midpoint
  ! T changes more than sat_T_tol between passes.  Each round honors max_sat_iter.
  integer,  parameter :: max_sat_iter = 50          ! cap on inner convadj passes
  real(r8), parameter :: sat_T_tol    = 1.0e-4_r8   ! max |ΔT| per convadj pass [K]

  ! Outer thermodynamic consistency loop: convadj and satadj are iterated
  ! together until the column simultaneously satisfies lapse-rate stability
  ! AND q ≤ qsat everywhere.  Converged when satadj removes no vapour (exact
  ! zero is achievable with Newton satadj).  Starting each radiation step from
  ! this self-consistent state reduces step-to-step T-profile variability and
  ! hence instantaneous TOA flux noise.
  integer,  parameter :: max_inner_phys = 20         ! cap on convadj+satadj cycles

  ! Slab-ocean heat capacity.  Density and cp are fixed (seawater); thickness
  ! is configurable via &exocol_nml::dz_slab.  H_slab is computed at run start.
  real(r8), parameter :: rho_w     = 1026._r8     ! seawater density [kg/m³]
  real(r8), parameter :: cp_w      = 4000._r8     ! seawater specific heat [J/kg/K]
  real(r8), save      :: H_slab                   ! ρ_w·cp_w·dz_slab [J/m²/K]

  ! Implicit slab Planck damping: the slab budget is integrated with a
  ! semi-implicit Euler step that incorporates the surface Planck feedback
  ! 4·σ·Ts³ on the right-hand side:
  !   Ts_new = Ts + dt · F_srf_total / (H_slab + dt · 4·σ·Ts³)
  ! For thick slabs (H_slab ≫ dt·λ) this reduces to the explicit Euler step;
  ! for thin slabs it suppresses overshoot when |F_srf_total| is large.  Only
  ! the Planck contribution is included — water-vapor and other feedbacks are
  ! handled implicitly through the explicit re-evaluation of LE/SH next step.

  ! Fixed-RH legacy closure (only used when moisture_scheme='fixed_rh')
  real(r8), parameter :: tau_relax = 50._r8       ! moisture relaxation [days]

contains

  subroutine run_rce_loop()
  ! Main RCE iteration.  Column state in exocol_mod is updated in-place.

    real(r8), dimension(pver)  :: LWHR, SWHR
    real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

    real(r8), dimension(pver) :: rh_init             ! fixed-RH legacy target
    real(r8) :: max_hr, toa_signed, toa_flux
    real(r8) :: F_net_srf_rad
    real(r8) :: LE, SH
    real(r8) :: precip_total                         ! column precip mass flux [kg/m²/s]
    real(r8) :: precip_iter                          ! precip from one condensation call

    ! ---- Radiation sub-cycling state ----
    real(r8) :: dt_outer_days                        ! outer (radiation) step
    real(r8) :: dt_sub_days, dt_sub_sec              ! sub-step = dt_outer / N_sub
    real(r8) :: LE_acc_outer, SH_acc_outer           ! dt-weighted flux sums over sub-steps
    real(r8) :: precip_acc_outer                     ! dt-weighted precip sum over sub-steps
    integer  :: N_sub, isub_rad

    real(r8) :: C_drag_diag                          ! MOS surface drag coefficient

    ! Window-averaged Path B accumulators.  Each step adds dt-weighted state to
    ! the sums; at end of window (model_time_days − t_window_start ≥ stab_check_days)
    ! we form time-means and compare to the previous window's means.  This
    ! suppresses step-level flicker (LE/precip/TOA spikes from discrete
    ! saturation events) in the convergence test — what matters is whether the
    ! CLIMATE (time-mean state) is steady, not whether instantaneous diagnostics
    ! are tightly held.
    real(r8), dimension(pver) :: tmid_sum, h2ommr_sum
    real(r8), dimension(pver) :: tmid_mean_prev, h2ommr_mean_prev
    real(r8) :: ts_sum, toa_signed_sum, window_dt_sum
    real(r8) :: ts_mean_prev, toa_mean_prev
    real(r8) :: ts_mean, toa_mean
    real(r8), dimension(pver) :: tmid_mean, h2ommr_mean
    real(r8) :: max_dT_mean, max_dq_mean, dTs_mean, dToa_mean
    logical  :: have_prev_window

    ! Window-mean flux accumulators for the final budget diagnostic.
    ! Parallels the ts/toa/tmid/h2ommr accumulators so the budget is evaluated
    ! from the same time-mean as the convergence criteria.
    real(r8) :: LE_sum, SH_sum, F_srf_rad_sum, precip_day_sum
    real(r8) :: LE_win, SH_win, F_srf_rad_win, precip_win

    ! Model-time tracking and last-snapshot times for adaptive-dt accounting.
    real(r8) :: model_time_days              ! cumulative model time [days]
    real(r8) :: t_window_start               ! model time at start of current window
    real(r8) :: t_last_print                 ! model time of last console print

    integer  :: it, k, inner_iter
    integer  :: sat_warn_count               ! count of outer steps where convadj hit max_sat_iter
    real(r8) :: f_zm                         ! ZM relaxation fraction for current step
    logical  :: converged, profile_stable
    logical  :: prognostic, fixed_rh

    ! Compute slab heat capacity from configured thickness.
    H_slab = rho_w * cp_w * cfg_dz_slab

    ! Initialize state.  dt_days is a module variable (the SUB-step dt, used by
    ! condense for its heating-rate diagnostic); reset here so re-entering the
    ! loop after a previous run starts fresh.
    dt_days          = dt_max
    dt_outer_days    = dt_max
    model_time_days  = 0._r8
    t_window_start   = 0._r8
    t_last_print     = -print_every_days     ! print at it=1
    converged        = .false.
    profile_stable   = .false.
    have_prev_window = .false.
    tmid_sum         = 0._r8
    h2ommr_sum       = 0._r8
    ts_sum           = 0._r8
    toa_signed_sum   = 0._r8
    window_dt_sum    = 0._r8
    ts_mean_prev     = 0._r8
    toa_mean_prev    = 0._r8
    tmid_mean_prev   = 0._r8
    h2ommr_mean_prev = 0._r8
    LE = 0._r8;  SH = 0._r8
    C_drag_diag      = C_D
    sat_warn_count   = 0
    LE_sum           = 0._r8
    SH_sum           = 0._r8
    F_srf_rad_sum    = 0._r8
    precip_day_sum   = 0._r8
    LE_win           = 0._r8
    SH_win           = 0._r8
    F_srf_rad_win    = 0._r8
    precip_win       = 0._r8

    prognostic = (trim(adjustl(moisture_scheme)) == 'prognostic')
    fixed_rh   = (trim(adjustl(moisture_scheme)) == 'fixed_rh')

    write(*,'(/,a)')      '========================================'
    write(*,'(a)')        ' ExoColumn RCE loop starting'
    write(*,'(a,f5.2,a,f7.4,a,f5.2,a,f4.2)') &
      '   adaptive dt:  dt_max =', dt_max, ' d   dt_min =', dt_min, &
      ' d   dT_target =', dT_target, ' K   safety =', cfl_safety
    write(*,'(a,i0)')     '   nmax       = ', nmax
    write(*,'(a,f7.4)')   '   hr_tol  [K/day] = ', hr_tol
    write(*,'(a,f6.3)')   '   toa_tol [W/m2]  = ', toa_tol

    write(*,'(a,f6.2,a,es9.2,a)') &
      '   slab: dz =', cfg_dz_slab, ' m   H_slab =', H_slab, ' J/m²/K'
    write(*,'(a)') '   condensation     : saturation adjustment (Newton)'
    write(*,'(a,/)')      '========================================'

    ! Capture initial RH only if the legacy fixed-RH closure is selected.
    if (fixed_rh) then
      call capture_rh_init(rh_init)
    end if

    do it = 1, nmax

      ! ---- 1. Radiation tendency (ONCE per outer step; frozen over sub-steps) ----
      call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)
      max_hr     = maxval(abs(LWHR(:) + SWHR(:)))
      toa_signed = SWDN(1) - SWUP(1) + LWDN(1) - LWUP(1)
      toa_flux   = abs(toa_signed)
      F_net_srf_rad = (SWDN(pverp) - SWUP(pverp)) + (LWDN(pverp) - LWUP(pverp))

      ! ---- DIAGNOSTIC (throwaway): radiative stiffness + Jacobian probe ----
      ! Fires once the stiff cold-point layer is well-formed, dumps the profile
      ! and probes J=∂HR/∂T, then stops.  Remove after the resolution study.
      if (diag_jacobian .and. it >= 100 .and. max_hr > 8.0_r8) then
        call diagnose_radiative_stiffness(LWHR, SWHR, it)
        stop 'diagnostic complete'
      end if

      ! ---- 2. Outer (radiation) timestep — keyed to the BULK radiative CFL ----
      ! Size dt_outer so the per-step ΔT of the climatically-relevant (mass-
      ! bearing) layers is bounded by dT_target, EXCLUDING the stiffest layers
      ! that together hold less than f_excl_mass of the column mass.  A thin
      ! cold-point / cold-trap layer concentrates a sharp flux divergence into
      ! negligible mass and has a large local |HR| that cannot drive the slow
      ! climate; excluding it keeps dt_outer at the bulk rate (≈ the reference
      ! scheme at PVER=70) instead of throttling to that one layer's CFL.  The
      ! excluded layer is still integrated stably below via the sub-cycle.
      block
        real(r8) :: hrabs(pver), excl_budget, acc_mass, hr_for_dt
        logical  :: avail(pver)
        integer  :: kk, kmax
        hrabs       = abs(LWHR + SWHR)
        excl_budget = f_excl_mass * sum(pdel)
        avail       = .true.
        acc_mass    = 0._r8
        hr_for_dt   = maxval(hrabs)
        do
          kmax = 0
          do kk = 1, pver
            if (avail(kk)) then
              if (kmax == 0) then
                kmax = kk
              else if (hrabs(kk) > hrabs(kmax)) then
                kmax = kk
              end if
            end if
          end do
          if (kmax == 0) then
            hr_for_dt = 0._r8;  exit              ! all layers excluded (unreachable)
          end if
          ! Stop at the stiffest layer whose mass would overrun the exclusion
          ! budget — that layer is KEPT and sets the bulk rate.
          if (acc_mass + pdel(kmax) > excl_budget) then
            hr_for_dt = hrabs(kmax);  exit
          end if
          acc_mass    = acc_mass + pdel(kmax)     ! exclude this thin stiff layer
          avail(kmax) = .false.
        end do

        if (hr_for_dt > 0._r8) then
          dt_outer_days = cfl_safety * dT_target / hr_for_dt
        else
          dt_outer_days = dt_max
        end if
      end block
      dt_outer_days   = max(dt_min, min(dt_max, dt_outer_days))
      model_time_days = model_time_days + dt_outer_days

      ! ---- 3. Sub-step count from the TRUE max radiative CFL ----
      ! Each sub-step applies frozen HR explicitly: require dt_sub·max|HR| ≤
      ! dT_target so even the excluded stiff layer integrates stably.  N_sub = 1
      ! whenever there is no stiff outlier (dt_outer already meets the CFL).
      if (max_hr > 0._r8) then
        N_sub = ceiling(dt_outer_days * max_hr / (cfl_safety * dT_target))
      else
        N_sub = 1
      end if
      N_sub       = max(1, min(nsub_max, N_sub))
      dt_sub_days = dt_outer_days / real(N_sub, r8)
      dt_sub_sec  = dt_sub_days * SHR_CONST_CSEC
      dt_days     = dt_sub_days     ! module var → condense heating-rate diagnostic

      ! Reset outer-step flux / precip accumulators.
      LE_acc_outer     = 0._r8
      SH_acc_outer     = 0._r8
      precip_acc_outer = 0._r8

      ! ================= Radiation sub-cycle (frozen HR) =================
      do isub_rad = 1, N_sub

      ! ---- 3a. Apply frozen radiation tendency (small, CFL-bounded step) ----
      do k = 1, pver
        tmid(k) = tmid(k) + dt_sub_days * (LWHR(k) + SWHR(k))
      end do

      ! ---- 4–6. Surface turbulent fluxes + slab ocean ----
      ! Compute LE, SH and the (MOS) drag coefficient from the current near-
      ! surface state and step the slab ocean (semi-implicit Planck damping).
      ! The mechanical (MOS) sensible and latent fluxes are deposited explicitly
      ! into the bottom grid layer; the log grid's bottom layer is thick enough
      ! that a single step is stable.  (The convective surface sensible flux that
      ! keeps the sub-cloud layer coupled to Ts is handled separately by
      ! convadj_surface in the convection block.)
      block
        real(r8) :: za_bot, p0_sfc, lambda_s, F_sfc, layer_mass_bot
        za_bot = 0.5_r8 * (zint(pver) + zint(pverp))   ! lowest-level height
        p0_sfc = pint(pverp)                           ! surface pressure

        call compute_surface_fluxes(surface_flux, ts, tmid(pver), h2ommr(pver), &
                                    pmid(pver), za_bot, p0_sfc, z0_rough, &
                                    mwdry_col, cpdry_col, wind_speed, C_D, &
                                    LE, SH, C_drag_diag)

        ! Slab budget (semi-implicit Planck)
        if (prognostic) then
          F_sfc = F_net_srf_rad - LE - SH
        else
          F_sfc = F_net_srf_rad
        end if
        lambda_s = 4._r8 * SHR_CONST_STEBOL * ts**3
        ts       = ts + dt_sub_sec * F_sfc / (H_slab + dt_sub_sec * lambda_s)

        ! Explicit deposit of the mechanical surface fluxes into the bottom layer.
        if (prognostic) then
          layer_mass_bot = pdel(pver) / gravity
          tmid(pver)   = tmid(pver) + dt_sub_sec * SH / (cpdry_col * layer_mass_bot)
          h2ommr(pver) = max(h2ommr(pver) &
                         + dt_sub_sec * (LE / Lvap_T(ts)) / layer_mass_bot, 0._r8)
        end if
      end block

      ! ---- 7. Legacy fixed-RH closure ----
      if (fixed_rh) then
        call update_h2ommr_fixed_rh(rh_init)
      end if

      ! ---- 8. Update derived (q has changed in either branch) ----
      call exocol_update_derived()

      ! ---- 9. Interface temperatures ----
      call update_tint()

      ! ---- 10. Thermodynamic consistency loop: convadj + satadj ----
      precip_total = 0._r8
      cond_heating = 0._r8

      if (trim(adjustl(conv_scheme)) == 'sbm') then
        ! Simplified Betts-Miller (Frierson 2007): relax T to the moist adiabat
        ! and q to rh_sbm·qsat over tau_conv, energy-conservingly.  A light
        ! saturation-adjustment mop-up then removes any residual (stratiform)
        ! supersaturation outside the convecting column.
        block
          real(r8) :: precip_sbm, ts_pre_couple
          real(r8), dimension(pver) :: cond_tend_sbm
          ! Surface-coupled mixed layer FIRST, so SBM lifts its moist adiabat from
          ! a bottom layer that is coupled to the surface temperature (warm,
          ! ventilated parcel base) rather than a radiatively-decoupled cold one.
          if (use_surf_couple) then
            ts_pre_couple = ts
            call convadj_surface(tmid, tint, h2ommr, ts, H_slab, &
                                 pint, pmid, pdel, cpdry_col, gravity, pver)
            ! The mixed-layer coupling carries a convective surface sensible flux
            ! (the slab debit H_slab·Δts).  Fold it into the SH diagnostic so the
            ! reported surface energy budget closes (F_srf_rad − LE − SH ≈ 0); it
            ! is a real slab→atmosphere flux, just realised by convection rather
            ! than the mechanical bulk/MOS formula.
            SH = SH + H_slab * (ts_pre_couple - ts) / dt_sub_sec
          end if
          call convadj_sbm(tmid, tint, h2ommr, zint, pint, pdel, &
                           cpdry_col, gravity, ts, dt_sub_sec, tau_conv, rh_sbm, &
                           Lvap_T(ts), prognostic, precip_sbm, cond_tend_sbm, pver)
          do k = 1, pver
            cond_heating(k) = cond_heating(k) + cond_tend_sbm(k) * SHR_CONST_CSEC
          end do
          precip_total = precip_sbm
          ! Non-deadlocking fallback: dry-adjust any dry-unstable layers SBM left
          ! behind (e.g. when its net-subsaturation gate tripped).  No-op on a
          ! moist-adiabatic / radiatively stable profile, so the reference
          ! equilibrium is untouched; prevents the cold-dry collapse otherwise.
          ! Runs before condense so any q homogenised into a colder layer is
          ! cleaned up to q <= qsat in the same step.
          if (use_dry_fallback) &
            call convadj_dry(tmid, tint, h2ommr, pint, pdel, &
                             cpdry_col, gravity, ts, pver)
          if (prognostic) then
            call condense(dt_sub_sec, precip_iter)
            precip_total = precip_total + precip_iter
          end if
        end block
      else if (trim(adjustl(conv_scheme)) == 'zm') then
        ! ZM soft scheme: one relaxed pass, condensation, then one hard cleanup.
        ! f_zm = 1 - exp(-dt/τ_conv): fraction of instability removed this step.
        ! At dt >> τ_conv (cold start) f_zm → 1 (hard); near equilibrium where
        ! dt_sub ≈ 0.06 d and τ_conv = 7200 s = 0.083 d, f_zm ≈ 0.51.
        f_zm = 1.0_r8 - exp(-dt_sub_days / (tau_conv / SHR_CONST_CSEC))
        call convadj_zm(tmid, tint, h2ommr, zint, pint, pdel, &
                        cpdry_col, gravity, ts, f_zm, cape_trigger, pver)
        if (prognostic) then
          call condense(dt_sub_sec, precip_iter)
          precip_total = precip_iter
        end if
        ! Hard cleanup (f=1, no CAPE check) removes any instability from
        ! latent heat release — this is always a fast local process.
        call convadj_zm(tmid, tint, h2ommr, zint, pint, pdel, &
                        cpdry_col, gravity, ts, 1.0_r8, 0.0_r8, pver)
      else
        ! Hard schemes: iterate convadj + satadj until q ≤ qsat everywhere.
        do inner_iter = 1, max_inner_phys
          call sat_iter_convadj(sat_warn_count)
          if (prognostic) then
            call condense(dt_sub_sec, precip_iter)
            precip_total = precip_total + precip_iter
            if (precip_iter == 0._r8) exit
          else
            exit
          end if
        end do
        call sat_iter_convadj(sat_warn_count)  ! final cleanup after last satadj
      end if

      ! ---- 10b. Stratospheric cold-point cold trap (Brewer-Dobson) ----
      ! ExoColumn has no vertical moisture transport above the convective top,
      ! so the stratosphere otherwise stays at its initialised q = 0, removing
      ! the principal stratospheric LW coolant.  Set stratospheric H2O from the
      ! cold-point freeze-drying value (see subroutine header).
      if (prognostic .and. use_coldtrap) call apply_stratospheric_coldtrap()

      ! ---- 11. Final derived update (q + T changed in sat-adjust loop) ----
      call exocol_update_derived()

      ! ---- 12. Heights ----
      call update_zint()

      ! ---- Accumulate sub-step fluxes / precip (dt-weighted) for outer means ----
      LE_acc_outer     = LE_acc_outer     + LE           * dt_sub_days
      SH_acc_outer     = SH_acc_outer     + SH           * dt_sub_days
      precip_acc_outer = precip_acc_outer + precip_total * dt_sub_days

      end do
      ! ================= end radiation sub-cycle =================

      ! Outer-step time-mean fluxes / precip for diagnostics + window accumulators.
      LE           = LE_acc_outer     / dt_outer_days
      SH           = SH_acc_outer     / dt_outer_days
      precip_total = precip_acc_outer / dt_outer_days

      ! ---- Diagnostics for output ----
      LE_diag     = LE
      SH_diag     = SH
      precip_diag = precip_total * 86400._r8   ! kg/m²/s → mm/day (ρ_water = 1000 kg/m³)

      ! ---- Window-mean accumulators (Path B convergence basis) ----
      ts_sum         = ts_sum         + ts            * dt_outer_days
      toa_signed_sum = toa_signed_sum + toa_signed    * dt_outer_days
      tmid_sum       = tmid_sum       + tmid          * dt_outer_days
      h2ommr_sum     = h2ommr_sum     + h2ommr        * dt_outer_days
      window_dt_sum  = window_dt_sum  + dt_outer_days
      LE_sum         = LE_sum         + LE            * dt_outer_days
      SH_sum         = SH_sum         + SH            * dt_outer_days
      F_srf_rad_sum  = F_srf_rad_sum  + F_net_srf_rad * dt_outer_days
      precip_day_sum = precip_day_sum + precip_diag   * dt_outer_days

      ! ---- End-of-window: form means, compare to previous, reset ----
      if (model_time_days - t_window_start >= stab_check_days) then
        ts_mean       = ts_sum         / window_dt_sum
        toa_mean      = toa_signed_sum / window_dt_sum
        tmid_mean     = tmid_sum       / window_dt_sum
        h2ommr_mean   = h2ommr_sum     / window_dt_sum
        LE_win        = LE_sum         / window_dt_sum
        SH_win        = SH_sum         / window_dt_sum
        F_srf_rad_win = F_srf_rad_sum  / window_dt_sum
        precip_win    = precip_day_sum / window_dt_sum

        if (have_prev_window) then
          dTs_mean    = abs(ts_mean  - ts_mean_prev)
          dToa_mean   = abs(toa_mean - toa_mean_prev)
          max_dT_mean = maxval(abs(tmid_mean   - tmid_mean_prev))
          max_dq_mean = maxval(abs(h2ommr_mean - h2ommr_mean_prev))
          profile_stable = ( max_dT_mean < prof_stab_tol .and. &
                             dTs_mean    < ts_stab_tol   .and. &
                             max_dq_mean < q_stab_tol    .and. &
                             abs(toa_mean) < toa_tol )
          write(*,'(a,f9.2,a,f7.3,a,f7.3,a,f7.3,a,es9.2,a,f7.4,a,es9.2)') &
            '  WINDOW@t=', model_time_days, &
            '  ⟨Ts⟩=', ts_mean, &
            '  ⟨TOA⟩=', toa_mean, &
            '  ΔTs_win=', dTs_mean, &
            '  ΔTOA_win=', dToa_mean, &
            '  ΔTpro=', max_dT_mean, &
            '  Δq=', max_dq_mean
          write(*,'(a,f8.3,a,f8.3,a,f8.3,a,f8.3,a,f8.3)') &
            '         BUDGET: ⟨TOA⟩=', toa_mean, &
            '  ⟨Fsrf_net⟩=', F_srf_rad_win - LE_win - SH_win, &
            '  (Frad=', F_srf_rad_win, ' LE=', LE_win, ' SH=', SH_win, ')'
          write(*,'(a,f8.4,a,f8.4,a,f8.4,a,f8.3)') &
            '         WATER:  ⟨evap⟩=', LE_win / Lvap_T(ts) * 86400._r8, &
            ' mm/d  ⟨precip⟩=', precip_win, &
            ' mm/d  imbal=', LE_win / Lvap_T(ts) * 86400._r8 - precip_win, &
            ' mm/d  → leak[W/m2]=', &
            (precip_win - LE_win / Lvap_T(ts) * 86400._r8) / 86400._r8 * Lvap_T(ts)
          flush(6)
        end if

        ts_mean_prev     = ts_mean
        toa_mean_prev    = toa_mean
        tmid_mean_prev   = tmid_mean
        h2ommr_mean_prev = h2ommr_mean
        have_prev_window = .true.

        ts_sum         = 0._r8
        toa_signed_sum = 0._r8
        tmid_sum       = 0._r8
        h2ommr_sum     = 0._r8
        window_dt_sum  = 0._r8
        LE_sum         = 0._r8
        SH_sum         = 0._r8
        F_srf_rad_sum  = 0._r8
        precip_day_sum = 0._r8
        t_window_start = model_time_days
      end if

      ! ---- Console print (model-time based) ----
      if (it == 1 .or. model_time_days - t_last_print >= print_every_days) then
        write(*,'(a,i7,a,f9.2,a,f7.4,a,i5,a,f7.3,a,f8.3,a,f7.2,a,f6.1,a,f6.1,a,f6.2)') &
          '  step=', it, &
          '  t[d]=', model_time_days, &
          '  dt[d]=', dt_outer_days, &
          '  Nsub=', N_sub, &
          '  max|HR|=', max_hr, &
          '  TOA=', toa_flux, &
          '  Ts=', ts, &
          '  LE=', LE, &
          '  SH=', SH, &
          '  P[mm/d]=', precip_diag
        flush(6)     ! live progress when stdout is redirected to a file
        t_last_print = model_time_days
      end if

      ! ---- Convergence checks ----
      if (max_hr < hr_tol .and. toa_flux < toa_tol) then
        converged = .true.
        exit
      end if

      if (profile_stable) then
        converged = .true.
        exit
      end if

    end do

    if (converged) then
      write(*,'(/,a,i0,a,f9.2,a)') &
        ' ExoColumn: converged after ', it, ' steps (', model_time_days, ' model days).'
    else
      write(*,'(/,a,i0,a,f9.2,a)') &
        ' ExoColumn: WARNING — not converged after ', nmax, ' steps (', &
        model_time_days, ' model days).'
    end if
    write(*,'(a,f8.3,a,f8.3)') &
      '   Final max|HR| [K/day] = ', max_hr, &
      '   Final TOA net [W/m2]  = ', toa_flux

    ! ---- Column energy & water budget closure (steady-state diagnostics) ----
    ! Uses window-mean fluxes (Path B) so step-level flicker doesn't inflate
    ! the residual.  Falls back to instantaneous values for Path A.
    if (have_prev_window) then
      write(*,'(a,f8.3,a,f8.3,a,f7.3)') &
        '   ⟨LE⟩  [W/m2]          = ', LE_win, &
        '   ⟨SH⟩ [W/m2] = ', SH_win, &
        '   ⟨P⟩ [mm/day] = ', precip_win
      block
        real(r8) :: srf_net_win, energy_residual, evap_win, water_residual
        srf_net_win    = F_srf_rad_win - LE_win - SH_win
        energy_residual = toa_mean - srf_net_win    ! both ≈ 0 at steady state
        evap_win       = LE_win / Lvap_T(ts) * 86400._r8
        water_residual = evap_win - precip_win
        write(*,'(a,f8.3,a,f8.3,a,f8.3)') &
          '   ⟨Energy budget⟩ [W/m2]: ⟨TOA⟩ = ', toa_mean, &
          '   ⟨F_srf_net⟩ = ', srf_net_win, &
          '   residual = ', energy_residual
        write(*,'(a,f7.3,a,f7.3,a,f7.3)') &
          '   ⟨Water budget⟩ [mm/day]: ⟨evap⟩ = ', evap_win, &
          '   ⟨precip⟩ = ', precip_win, &
          '   residual = ', water_residual
      end block
      if (abs(toa_mean) > toa_tol) then
        write(*,'(a,f6.3,a)') &
          '   WARNING: ⟨TOA⟩ = ', toa_mean, ' W/m² exceeds toa_tol — energy budget not fully closed.'
      end if
    else
      write(*,'(a,f8.3,a,f8.3,a,f7.3)') &
        '   Final LE [W/m2]       = ', LE_diag, &
        '   SH [W/m2] = ', SH_diag, &
        '   P [mm/day] = ', precip_diag
      block
        real(r8) :: toa_net, srf_net, energy_residual, evap_mm_day, water_residual
        toa_net        = SWDN(1) - SWUP(1) + LWDN(1) - LWUP(1)
        srf_net        = F_net_srf_rad - LE_diag - SH_diag
        energy_residual = toa_net - srf_net
        evap_mm_day    = LE_diag / Lvap_T(ts) * 86400._r8
        water_residual = evap_mm_day - precip_diag
        write(*,'(a,f8.3,a,f8.3,a,f8.3)') &
          '   Energy budget [W/m2]  : F_TOA = ', toa_net, &
          '   F_srf_net = ', srf_net, &
          '   residual = ', energy_residual
        write(*,'(a,f7.3,a,f7.3,a,f7.3)') &
          '   Water budget [mm/day] : evap = ', evap_mm_day, &
          '   precip = ', precip_diag, &
          '   residual = ', water_residual
      end block
      if (toa_flux > toa_tol) then
        write(*,'(a)') &
          '   WARNING: TOA net flux exceeds toa_tol — column may be convectively' // &
          ' active; structural imbalance possible.'
      end if
    end if
    if (sat_warn_count > 0) then
      write(*,'(a,i0,a,i0,a)') &
        '   WARNING: saturation adjustment hit max_sat_iter (', max_sat_iter, &
        ') on ', sat_warn_count, ' step(s) — consider raising max_sat_iter.'
    end if

  end subroutine run_rce_loop

  ! -----------------------------------------------------------------------
  ! Convadj fixed-point iteration (no condensation inside)
  ! -----------------------------------------------------------------------

  subroutine sat_iter_convadj(warn_count)
  ! Iterate the selected convadj scheme until the column lapse rate is stable
  ! (no layer-midpoint temperature change above sat_T_tol between passes).
  ! convadj routines mix q within adjusted pairs but do not condense — that
  ! is the caller's responsibility.  Typical convergence: 2–4 passes when
  ! called from a settled column; up to ~10 during cold-start transients.
  !
  ! ifx host-association workaround: explicit USE inside the contained
  ! subroutine (matches update_zint below).
    use exocol_mod, only: tmid, tint, h2ommr, zint, pint, pdel, cpdry_col, ts
    integer, intent(inout) :: warn_count

    real(r8), dimension(pver) :: tmid_pre
    real(r8) :: dT_sat
    integer  :: sat_iter
    logical  :: sat_converged

    sat_converged = .false.
    do sat_iter = 1, max_sat_iter
      tmid_pre = tmid

      select case (trim(adjustl(conv_scheme)))
      case ('moist')
        call convadj_moist(tmid, tint, h2ommr, zint, pint, pdel, &
                           cpdry_col, gravity, ts, pver)
      case ('manabe')
        call convadj_manabe(tmid, tint, h2ommr, zint, pint, pdel, &
                            cpdry_col, gravity, ts, pver)
      case ('zm')
        ! ZM is dispatched in run_rce_loop; this fallback applies one hard pass.
        call convadj_zm(tmid, tint, h2ommr, zint, pint, pdel, &
                        cpdry_col, gravity, ts, 1.0_r8, 0.0_r8, pver)
      case default  ! 'dry'
        call convadj_dry(tmid, tint, h2ommr, pint, pdel, &
                         cpdry_col, gravity, ts, pver)
      end select

      dT_sat = maxval(abs(tmid - tmid_pre))
      if (dT_sat < sat_T_tol) then
        sat_converged = .true.
        exit
      end if
    end do
    if (.not. sat_converged) warn_count = warn_count + 1
  end subroutine sat_iter_convadj

  ! -----------------------------------------------------------------------
  ! Condensation / precipitation with phase-aware latent-heat release
  ! -----------------------------------------------------------------------

  subroutine condense(dt_step_sec, precip_mass_flux)
  ! Saturation adjustment: where h2ommr > qsat(T,p), find the simultaneous
  ! (T_new, q_new) that satisfy energy conservation and exact saturation.
  !
  ! Newton solve for Δq (one step per layer, no outer iteration needed):
  !   T_new  = T + (L/cp) · Δq                 [enthalpy conservation]
  !   q_new  = q − Δq  =  qsat(T_new)          [exact saturation]
  !
  ! Linearising qsat(T_new) ≈ qsat(T) + (dqsat/dT)·(L/cp)·Δq and solving:
  !   Δq = (q − qsat(T)) / (1 + (L/cp) · dqsat/dT)
  ! where dqsat/dT = qsat·(1 + qsat/ε)·L/(Rv·T²)  from Clausius-Clapeyron,
  ! and ε = Mwv/Mdry.
  !
  ! Correctness of single step: the Clausius-Clapeyron curve is convex
  ! (d²qsat/dT² > 0), so qsat(T_new) ≥ the linearised value = q_new.
  ! The Newton step therefore always leaves q_new ≤ qsat(T_new) — the layer
  ! is guaranteed subsaturated or exactly saturated; no iteration is needed.
  !
  ! L is Lvap_T(ts) throughout (matching the surface evap ledger): every kg
  ! of vapour was placed in the column by surface evaporation that debited
  ! L(ts) from the slab; releasing L(tmid) at condensation would inject a
  ! free L_fusion per kg cycled through the ice phase.  Phase-aware qsat
  ! (esat_cc uses L_sub below freezing) controls WHEN condensation triggers,
  ! independently of how much heat is released per kg.
    real(r8), intent(in)  :: dt_step_sec
    real(r8), intent(out) :: precip_mass_flux

    integer  :: k
    real(r8) :: eps_wv, es_k, qsat_k, dqsat_dT, q_excess, L_release

    eps_wv    = SHR_CONST_MWWV / mwdry_col
    L_release = Lvap_T(ts)
    precip_mass_flux = 0._r8

    do k = 1, pver
      es_k   = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
      qsat_k = eps_wv * es_k / (pmid(k) - es_k)
      if (h2ommr(k) > qsat_k) then
        dqsat_dT  = qsat_k * (1._r8 + qsat_k / eps_wv) &
                    * L_release / (SHR_CONST_RWV * tmid(k)**2)
        q_excess  = (h2ommr(k) - qsat_k) / (1._r8 + (L_release / cpdry_col) * dqsat_dT)
        h2ommr(k) = h2ommr(k) - q_excess
        tmid(k)   = tmid(k)   + L_release * q_excess / cpdry_col
        cond_heating(k)  = cond_heating(k) &
                           + L_release * q_excess / cpdry_col / dt_days
        precip_mass_flux = precip_mass_flux + q_excess * pdel(k) / gravity / dt_step_sec
      end if
    end do

  end subroutine condense

  ! -----------------------------------------------------------------------
  ! Stratospheric cold-point cold trap (Brewer-Dobson freeze-drying)
  ! -----------------------------------------------------------------------

  subroutine apply_stratospheric_coldtrap()
  ! Physics-based stratospheric water-vapour closure.
  !
  ! Physical basis: tropospheric air can only enter the stratosphere by
  ! ascending through the cold-point tropopause, where it is dehydrated by
  ! condensation down to the saturation mixing ratio at the cold-point
  ! temperature ("freeze drying"; Brewer 1949).  Above the cold point the air
  ! continues to rise and warm, so it becomes subsaturated and no further
  ! condensation occurs — the water-vapour mixing ratio is therefore conserved
  ! (vertically uniform) throughout the stratosphere at the cold-point entry
  ! value.  This is the dominant control on stratospheric humidity and is the
  ! same closure used by konrad (ColdPointCoupling).
  !
  ! Why ExoColumn needs it: the model has no vertical moisture transport above
  ! the convective top (SBM only relaxes q within the convecting column;
  ! condense only removes vapour; surface evaporation sources only the bottom
  ! layer).  The stratosphere would otherwise remain at its initialised q = 0,
  ! removing the principal stratospheric LW coolant and leaving the stratosphere
  ! unphysically warm — which in turn forces a sharp single-layer tropopause
  ! where the cold convective top abuts the too-warm dry stratosphere.
  !
  ! Self-consistency: the entry value is recomputed from the current cold point
  ! every step, so it tracks the radiative-convective solution — a colder cold
  ! point gives a drier stratosphere (less H2O cooling) and vice versa.
  !
  ! Implementation: freeze-dried air at the cold point is the driest air in the
  ! whole column — every layer below is warmer and moister (the convectively-
  ! moistened troposphere), and every layer above conserves the same mixing ratio
  ! on ascent.  ABOVE the cold point we ASSIGN q = qsat(T_cp) (the conserved entry
  ! value, vertically uniform), so the upper stratosphere tracks the cold point
  ! both up and down.  AT/BELOW the cold point we apply qsat(T_cp) as a FLOOR:
  ! in the troposphere q >> qsat(T_cp) so it is a no-op, while in the gap between
  ! the convective top and the cold point it fills q = qsat(T_cp), avoiding a dry
  ! warm gap that would otherwise push the cold point unphysically high.  Since
  ! T >= T_cp everywhere, qsat >= q_cp, so no supersaturation is introduced.
  ! (A pure column-wide max() floor was used previously, but its one-way ratchet
  ! stranded stale higher humidity in the warm upper stratosphere when the cold
  ! point cooled, producing an unphysical moist bump above the cold point.)
    use exocol_mod,    only: tmid, pmid, pdel, h2ommr, mwdry_col, gravity
    use ppgrid,        only: pver
    integer  :: k, k_cp
    real(r8) :: eps_wv, es_cp, q_cp, q_old, q_new, dW_strat, W_trop, frac

    eps_wv = SHR_CONST_MWWV / mwdry_col

    ! Cold-point tropopause = coldest model level → freeze-drying value.
    k_cp  = minloc(tmid, dim=1)
    es_cp = min(esat_cc(tmid(k_cp)), 0.99_r8 * pmid(k_cp))
    q_cp  = eps_wv * es_cp / (pmid(k_cp) - es_cp)

    ! Above the cold point (k < k_cp): the conserved entry value is exactly q_cp,
    ! so ASSIGN it (not max).  A pure max() floor is a one-way ratchet — it can
    ! only raise q — so when the cold point cools over the run the warm upper
    ! stratosphere keeps a stale, higher humidity it can never shed (it is far
    ! subsaturated, so condense never removes it), producing an unphysical
    ! moist bump above the cold point.  Assigning lets the upper stratosphere
    ! track the cold point both up AND down.  Safe: q_cp = qsat(T_cp) ≤ qsat
    ! everywhere above (the cold point is the coldest level), so no supersaturation.
    !
    ! At and below the cold point (k >= k_cp): keep the floor.  It fills the dry
    ! gap between the convective top and the cold point at q_cp and is a no-op in
    ! the moist troposphere (q >> q_cp there).
    !
    ! WATER CONSERVATION (critical): the assign/floor below would otherwise CREATE
    ! water from nothing (it raises q in dry stratospheric/gap layers without a
    ! source).  That spurious vapour is later precipitated by the convection/
    ! condensation machinery, releasing ~9 W/m² of latent heat from nothing — a
    ! structural TOA energy leak (precip > evap; see the WATER diagnostic in the
    ! RCE loop).  The leak is small in a dry-stratosphere reference column but is
    ! amplified ~100x once the boundary layer mixes more moisture through the
    ! column.  Fix: accumulate the NET water the cold trap adds (dW_strat) and
    ! source it CONSERVATIVELY from the moist troposphere below — vapour is moved,
    ! not created (no phase change, no latent heat, column water unchanged), which
    ! is what the Brewer-Dobson circulation physically does.
    dW_strat = 0._r8
    do k = 1, pver
      q_old = h2ommr(k)
      if (k < k_cp) then
        q_new = q_cp                          ! assign (conserved entry value)
      else
        q_new = max(h2ommr(k), q_cp)          ! floor the dry gap
      end if
      dW_strat  = dW_strat + (q_new - q_old) * pdel(k) / gravity
      h2ommr(k) = q_new
    end do

    ! Debit dW_strat from the moist troposphere (q > q_cp), weighted by each
    ! layer's water content so the removal is a tiny uniform fraction that cannot
    ! drive any layer negative.  dW_strat < 0 (net stratospheric drying) adds the
    ! water back, so the redistribution is exactly water-conserving either way.
    W_trop = 0._r8
    do k = 1, pver
      if (h2ommr(k) > q_cp) W_trop = W_trop + h2ommr(k) * pdel(k) / gravity
    end do
    if (W_trop > dW_strat .and. W_trop > 0._r8) then
      frac = dW_strat / W_trop
      do k = 1, pver
        if (h2ommr(k) > q_cp) h2ommr(k) = h2ommr(k) * (1._r8 - frac)
      end do
    end if
  end subroutine apply_stratospheric_coldtrap

  ! -----------------------------------------------------------------------
  ! THROWAWAY DIAGNOSTIC: radiative stiffness + Jacobian structure probe
  ! -----------------------------------------------------------------------

  subroutine diagnose_radiative_stiffness(LWHR0, SWHR0, it_now)
  ! Dump the current column profile and probe the radiative Jacobian
  ! J(i,j) = ∂HR_i/∂T_j by +1 K finite difference at representative layers.
  ! Writes iofiles/diag_profile.txt and iofiles/diag_jacobian.txt and prints a
  ! console summary.  Purpose: measure where the stiff layer is and how far a
  ! single-layer T perturbation's HR response spreads (the band half-width),
  ! to decide whether a banded backward-Euler radiation solve is feasible/cheap.
    use exocol_mod, only: tmid, pmid, h2ommr, pdel
    use ppgrid,     only: pver, pverp
    real(r8), intent(in) :: LWHR0(pver), SWHR0(pver)
    integer,  intent(in) :: it_now

    real(r8) :: HR0(pver), Jcol(pver)
    real(r8) :: LWHR_p(pver), SWHR_p(pver)
    real(r8) :: LWUP(pverp), LWDN(pverp), SWUP(pverp), SWDN(pverp)
    real(r8) :: tsave, dTp, diagv, offmax, bandfrac
    integer  :: k, j, ip, kstiff, ku, kl, u, bw
    integer, parameter :: nprobe = 9
    integer :: probes(nprobe)

    HR0    = LWHR0 + SWHR0
    kstiff = maxloc(abs(HR0), dim=1)
    dTp    = 1.0_r8

    ! --- profile dump ---
    open(newunit=u, file='iofiles/diag_profile.txt', status='replace')
    write(u,'(a)') '# k  pmid[Pa]  tmid[K]  h2ommr[kg/kg]  pdel[Pa]  LWHR  SWHR  absHR[K/day]'
    do k = 1, pver
      write(u,'(i5,7es15.6)') k, pmid(k), tmid(k), h2ommr(k), pdel(k), &
                              LWHR0(k), SWHR0(k), abs(HR0(k))
    end do
    close(u)

    ! Probe layers: cluster around the stiff layer + a few references.
    probes = (/ max(1,kstiff-6), max(1,kstiff-2), max(1,kstiff-1), kstiff, &
                min(pver,kstiff+1), min(pver,kstiff+2), min(pver,kstiff+6), &
                max(1,kstiff/2), pver-1 /)

    write(*,'(/,a)') '================ JACOBIAN DIAGNOSTIC ================'
    write(*,'(a,i0,a,i0,a)') '  step=', it_now, '  pver=', pver, ''
    write(*,'(a,i0,a,es12.4,a,f8.3)') '  stiff layer kstiff=', kstiff, &
          '  pmid=', pmid(kstiff), ' Pa  HR=', HR0(kstiff)
    write(*,'(a)') '  probe j   pmid[Pa]    diag J(j,j)[1/day]  band(|J|>5%diag)  max|offdiag|/|diag|'

    open(newunit=u, file='iofiles/diag_jacobian.txt', status='replace')
    write(u,'(a,i0)') '# kstiff=', kstiff
    write(u,'(a)') '# columns: i  pmid[i]  then J(i,j) for each probe j'
    write(u,'(a,9i14)') '# probe j list: ', (probes(ip), ip=1,nprobe)

    block
      real(r8) :: Jmat(pver, nprobe)
      do ip = 1, nprobe
        j        = probes(ip)
        tsave    = tmid(j)
        tmid(j)  = tmid(j) + dTp
        call exocol_rad_tend(LWHR_p, SWHR_p, LWUP, LWDN, SWUP, SWDN)
        tmid(j)  = tsave
        Jcol     = ((LWHR_p + SWHR_p) - HR0) / dTp
        Jmat(:,ip) = Jcol

        ! band half-width: furthest layer from j with |J| > 5% of |diag|
        diagv    = abs(Jcol(j))
        bandfrac = 0.05_r8 * max(diagv, 1.0e-30_r8)
        kl = j;  ku = j
        do k = 1, pver
          if (abs(Jcol(k)) > bandfrac) then
            if (k < kl) kl = k
            if (k > ku) ku = k
          end if
        end do
        bw = max(j-kl, ku-j)
        ! max off-diagonal magnitude relative to diagonal
        offmax = 0._r8
        do k = 1, pver
          if (k /= j) offmax = max(offmax, abs(Jcol(k)))
        end do
        write(*,'(i9,es13.4,es16.4,i14,f20.4)') j, pmid(j), Jcol(j), bw, &
              offmax / max(diagv,1.0e-30_r8)
      end do

      ! raw matrix dump
      do k = 1, pver
        write(u,'(i5,es14.5,9es14.5)') k, pmid(k), (Jmat(k,ip), ip=1,nprobe)
      end do
    end block
    close(u)

    write(*,'(a)') '  wrote iofiles/diag_profile.txt, iofiles/diag_jacobian.txt'
    write(*,'(a,/)') '===================================================='
    flush(6)
  end subroutine diagnose_radiative_stiffness

  ! -----------------------------------------------------------------------
  ! Legacy fixed-RH closure (only invoked when moisture_scheme='fixed_rh')
  ! -----------------------------------------------------------------------

  subroutine capture_rh_init(rh_init)
    real(r8), intent(out) :: rh_init(pver)
    integer  :: k
    real(r8) :: eps_wv, es_k, qsat_k
    eps_wv = SHR_CONST_MWWV / mwdry_col
    do k = 1, pver
      es_k       = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
      qsat_k     = eps_wv * es_k / (pmid(k) - es_k)
      rh_init(k) = min(h2ommr(k) / max(qsat_k, 1.0e-20_r8), 1.0_r8)
    end do
  end subroutine capture_rh_init

  subroutine update_h2ommr_fixed_rh(rh_init)
  ! Relax h2ommr toward rh_init(k)·qsat(T(k),p(k)) with timescale tau_relax.
    use exocol_mod, only: h2ommr, tmid, pmid, mwdry_col
    real(r8), intent(in) :: rh_init(pver)
    integer  :: k
    real(r8) :: eps_wv, es_k, qsat_k, alpha
    eps_wv = SHR_CONST_MWWV / mwdry_col
    alpha  = min(dt_days / tau_relax, 1.0_r8)
    do k = 1, pver
      es_k      = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
      qsat_k    = eps_wv * es_k / (pmid(k) - es_k)
      h2ommr(k) = h2ommr(k) + alpha * (rh_init(k) * qsat_k - h2ommr(k))
      h2ommr(k) = max(h2ommr(k), 0.0_r8)
    end do
  end subroutine update_h2ommr_fixed_rh

  ! -----------------------------------------------------------------------
  ! Private helpers — tint and zint recomputation
  ! -----------------------------------------------------------------------

  subroutine update_tint()
  ! ifx host-association workaround: explicit USE inside contained subroutine.
    use exocol_mod,     only: tint, tmid, pint, ts
    use ppgrid,         only: pver, pverp
    use exocol_convadj, only: compute_tint_interp
    tint(pverp) = ts
    call compute_tint_interp(tmid, pint, pver, tint)
  end subroutine update_tint

  subroutine update_zint()
  ! ifx host-association workaround: explicit USE inside the contained subroutine.
    use exocol_mod, only: zint, tmid, pint, mwdry_col, gravity
    integer  :: k
    real(r8) :: R_gas
    R_gas = SHR_CONST_RGAS / mwdry_col
    do k = pver, 1, -1
      zint(k) = zint(k+1) + (R_gas / gravity) * tmid(k) * log(pint(k+1) / pint(k))
    end do
  end subroutine update_zint

end module exocol_rce_loop
