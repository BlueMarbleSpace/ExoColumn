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
!   7. Saturation adjustment — convadj and condense iterated to fixed point.
!      convadj equilibrates the lapse rate; condense relaxes q toward qsat
!      with an implicit-Euler step of timescale τ_cond (replacing the legacy
!      instant cap), releasing Lvap_T(ts)·Δq into tmid.  The iteration
!      handles the latent-feedback into convadj.
!   8. Update derived (pdeldry, pintdry), tint, zint.
!
! CFL / stability summary (Earth-like column at dt_max = 1 d):
!   Radiative ΔT/step          → dt_max · max|HR|        bounded by dT_target.
!   Surface turbulent flux     → τ_surf ≈ 8 h            implicit damping.
!   Slab radiative equilibration → τ_slab ≈ 350 d        well below dt_max.
!   Condensation onset         → τ_cond ≈ 1 h            relaxation.
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
                             SHR_CONST_STEBOL
  use exoplanet_mod,   only: exo_g
  use ppgrid,          only: pver, pverp
  use exocol_mod
  use exocol_radiation, only: exocol_rad_tend
  use exocol_config,    only: conv_scheme, moisture_scheme, wind_speed, C_D, &
                              cfg_dz_slab => dz_slab
  use exocol_convadj,   only: convadj_dry, convadj_moist, convadj_manabe, &
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
  real(r8), parameter :: dT_target   = 2.0_r8       ! target |dt·HR| per step [K]
  real(r8), parameter :: cfl_safety  = 0.8_r8       ! safety factor [-]

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
  ! T changes more than sat_T_tol between passes.  Called twice per outer
  ! step: once before condense (initial cleanup) and once after (handle
  ! latent-release destabilization).  Each round honors max_sat_iter.
  integer,  parameter :: max_sat_iter = 50          ! cap on inner iterations
  real(r8), parameter :: sat_T_tol    = 1.0e-4_r8   ! max |ΔT| per convadj pass [K]

  ! Condensation relaxation timescale.  In condense(), the excess vapor above
  ! qsat is removed via an implicit-Euler step of timescale tau_cond:
  !   q_new = (q + (dt/τ)·qsat) / (1 + dt/τ)
  ! For dt ≫ τ this recovers the legacy instant cap (alpha → 1); for dt ≲ τ
  ! it spreads the condensation across multiple outer steps, eliminating the
  ! threshold-switch flicker that previously prevented Path B from firing at
  ! short dt.  τ ≈ 1 hour matches warm-cloud autoconversion timescales.
  real(r8), parameter :: tau_cond = 3600._r8        ! condensation timescale [s]

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
    real(r8) :: dt_sec
    real(r8) :: max_hr, toa_signed, toa_flux
    real(r8) :: F_net_srf_rad
    real(r8) :: LE, SH
    real(r8) :: precip_total                         ! column precip mass flux [kg/m²/s]
    real(r8) :: precip_iter                          ! precip from one condensation call

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

    integer  :: it, k
    integer  :: sat_warn_count               ! count of outer steps where convadj hit max_sat_iter
    logical  :: converged, profile_stable
    logical  :: prognostic, fixed_rh

    ! Compute slab heat capacity from configured thickness.
    H_slab = rho_w * cp_w * cfg_dz_slab

    ! Initialize state.  dt_days is a module variable; reset to dt_max here so
    ! re-entering the loop after a previous run starts fresh.
    dt_days          = dt_max
    dt_sec           = dt_days * SHR_CONST_CSEC
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
    write(*,'(a,f6.1)')   '   tau_cond [s]    = ', tau_cond
    write(*,'(a,f6.2,a,es9.2,a)') &
      '   slab: dz =', cfg_dz_slab, ' m   H_slab =', H_slab, ' J/m²/K'
    write(*,'(a,/)')      '========================================'

    ! Capture initial RH only if the legacy fixed-RH closure is selected.
    if (fixed_rh) then
      call capture_rh_init(rh_init)
    end if

    do it = 1, nmax

      ! ---- 1. Radiation tendency (diagnostic; applied after adaptive dt) ----
      call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)
      max_hr     = maxval(abs(LWHR(:) + SWHR(:)))
      toa_signed = SWDN(1) - SWUP(1) + LWDN(1) - LWUP(1)
      toa_flux   = abs(toa_signed)

      ! ---- 2. Adaptive timestep ----
      ! Size dt so the most active layer's per-step ΔT is bounded by dT_target,
      ! keeping explicit-Euler radiation in its linear regime.  In quiescent
      ! states max|HR| → 0 and dt saturates at dt_max.
      if (max_hr > 0._r8) then
        dt_days = cfl_safety * dT_target / max_hr
      else
        dt_days = dt_max
      end if
      dt_days = max(dt_min, min(dt_max, dt_days))
      dt_sec  = dt_days * SHR_CONST_CSEC
      model_time_days = model_time_days + dt_days

      ! ---- 3. Apply radiation tendency ----
      do k = 1, pver
        tmid(k) = tmid(k) + dt_days * (LWHR(k) + SWHR(k))
      end do

      ! ---- 4–6. Surface turbulent fluxes, slab, and bottom-layer deposit ----
      ! The coupled slab-atmosphere system is subcycled at dt_sfc ≤ τ_surf so
      ! that explicit Euler stays stable for any bottom-layer thickness.
      !
      ! τ_surf = (pdel/g) / (ρ·C_D·U) — the bottom-layer moisture relaxation
      ! time.  The old scheme used implicit-Euler damping 1/(1+dt/τ), which is
      ! unconditionally stable but collapses LE → dm·(qsat−q)/dt → 0 when
      ! τ_surf << dt (thin surface layers from the hybrid grid).  The subcycle
      ! avoids this: for thick layers τ_surf >> dt, n_sub = 1 and the cost is
      ! identical to the old scheme; for thin layers n_sub grows but each
      ! subcycle is only one cheap flux call.
      !
      ! F_net_srf_rad is held fixed at the value from the current radiation
      ! call (computed once, applied to every subcycle — valid because the
      ! outer radiation timestep is already longer than τ_slab ~ 350 d).
      F_net_srf_rad = (SWDN(pverp) - SWUP(pverp)) + (LWDN(pverp) - LWUP(pverp))

      block
        real(r8) :: Rd, rho_air_bot, layer_mass_bot, tau_surf
        real(r8) :: dt_sfc, ts_sub, T_sub, q_sub
        real(r8) :: LE_sub, SH_sub, LE_acc, SH_acc, L_sub, lambda_sub, F_sfc_sub
        integer  :: n_sub, isub

        Rd             = SHR_CONST_RGAS / mwdry_col
        rho_air_bot    = pmid(pver) / (Rd * tmid(pver))
        layer_mass_bot = pdel(pver) / exo_g
        tau_surf       = layer_mass_bot / (rho_air_bot * C_D * wind_speed)

        ! n_sub so that dt_sfc < tau_surf (CFL ≤ 1 for the bottom layer).
        ! For τ >> dt (thick layer), n_sub = 1 — single undamped explicit step.
        n_sub  = max(1, int(dt_sec / tau_surf) + 1)
        dt_sfc = dt_sec / real(n_sub, r8)

        ts_sub = ts
        T_sub  = tmid(pver)
        q_sub  = h2ommr(pver)
        LE_acc = 0._r8
        SH_acc = 0._r8

        do isub = 1, n_sub
          call compute_surface_fluxes(ts_sub, T_sub, q_sub, pmid(pver), &
                                      mwdry_col, cpdry_col, wind_speed, C_D, &
                                      LE_sub, SH_sub)

          ! Slab budget (semi-implicit Planck at sub-step scale)
          if (prognostic) then
            F_sfc_sub = F_net_srf_rad - LE_sub - SH_sub
          else
            F_sfc_sub = F_net_srf_rad
          end if
          lambda_sub = 4._r8 * SHR_CONST_STEBOL * ts_sub**3
          ts_sub     = ts_sub + dt_sfc * F_sfc_sub &
                              / (H_slab + dt_sfc * lambda_sub)

          ! Bottom-layer moisture and temperature sources
          if (prognostic) then
            L_sub = Lvap_T(ts_sub)
            T_sub = T_sub + dt_sfc * SH_sub / (cpdry_col * layer_mass_bot)
            q_sub = max(q_sub + dt_sfc * (LE_sub / L_sub) / layer_mass_bot, &
                        0._r8)
          end if

          LE_acc = LE_acc + LE_sub * dt_sfc
          SH_acc = SH_acc + SH_sub * dt_sfc
        end do

        ! Commit subcycled state to the column
        ts           = ts_sub
        tmid(pver)   = T_sub
        h2ommr(pver) = q_sub

        ! Time-averaged fluxes for diagnostics and window accumulators
        LE = LE_acc / dt_sec
        SH = SH_acc / dt_sec
      end block

      ! ---- 7. Legacy fixed-RH closure ----
      if (fixed_rh) then
        call update_h2ommr_fixed_rh(rh_init)
      end if

      ! ---- 8. Update derived (q has changed in either branch) ----
      call exocol_update_derived()

      ! ---- 9. Interface temperatures ----
      call update_tint()

      ! ---- 10. Saturation adjustment: convadj → condense → convadj ----
      ! With smooth condensation (one dt-step of implicit-Euler relaxation),
      ! condense() should run ONCE per outer step — its alpha already
      ! integrates dq/dt over the full dt.  We bracket condense with two
      ! rounds of convadj iteration:
      !   Pass A: equilibrate the lapse rate (convadj only).
      !   Pass B: one smooth condensation relaxation.
      !   Pass C: cleanup convadj after latent release re-destabilizes pairs.
      ! Each convadj round iterates internally until no T change > sat_T_tol.
      ! Both rounds use the same max_sat_iter cap.
      precip_total = 0._r8
      cond_heating = 0._r8

      call sat_iter_convadj(sat_warn_count)

      if (prognostic) then
        call condense(dt_sec, precip_iter)
        precip_total = precip_total + precip_iter
      end if

      call sat_iter_convadj(sat_warn_count)

      ! ---- 11. Final derived update (q + T changed in sat-adjust loop) ----
      call exocol_update_derived()

      ! ---- 12. Heights ----
      call update_zint()

      ! ---- Diagnostics for output ----
      LE_diag     = LE
      SH_diag     = SH
      precip_diag = precip_total * 86400._r8   ! kg/m²/s → mm/day (ρ_water = 1000 kg/m³)

      ! ---- Window-mean accumulators (Path B convergence basis) ----
      ts_sum         = ts_sum         + ts            * dt_days
      toa_signed_sum = toa_signed_sum + toa_signed    * dt_days
      tmid_sum       = tmid_sum       + tmid          * dt_days
      h2ommr_sum     = h2ommr_sum     + h2ommr        * dt_days
      window_dt_sum  = window_dt_sum  + dt_days
      LE_sum         = LE_sum         + LE            * dt_days
      SH_sum         = SH_sum         + SH            * dt_days
      F_srf_rad_sum  = F_srf_rad_sum  + F_net_srf_rad * dt_days
      precip_day_sum = precip_day_sum + precip_diag   * dt_days

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
        write(*,'(a,i7,a,f9.2,a,f7.4,a,f7.3,a,f8.3,a,f7.2,a,f6.1,a,f6.1,a,f6.2)') &
          '  step=', it, &
          '  t[d]=', model_time_days, &
          '  dt[d]=', dt_days, &
          '  max|HR|=', max_hr, &
          '  TOA=', toa_flux, &
          '  Ts=', ts, &
          '  LE=', LE, &
          '  SH=', SH, &
          '  P[mm/d]=', precip_diag
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
                           cpdry_col, exo_g, ts, pver)
      case ('manabe')
        call convadj_manabe(tmid, tint, h2ommr, zint, pint, pdel, &
                            cpdry_col, exo_g, ts, pver)
      case default  ! 'dry'
        call convadj_dry(tmid, tint, h2ommr, pint, pdel, &
                         cpdry_col, exo_g, ts, pver)
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
  ! Where h2ommr > qsat(T,p), relax h2ommr toward qsat on timescale tau_cond
  ! and release the corresponding latent heat into tmid.
  !
  ! Numerics — implicit Euler on dq/dt = -(q − qsat)/τ_cond when q > qsat:
  !   q_new = (q + (dt/τ)·qsat) / (1 + dt/τ)
  !   Δq   = q − q_new = (dt/τ)/(1 + dt/τ) · (q − qsat) = alpha · (q − qsat)
  ! For dt ≫ τ, alpha → 1 — recovers the legacy instant cap.
  ! For dt ≲ τ, alpha < 1 — spreads condensation across multiple outer steps,
  ! eliminating the threshold-switch flicker that prevents Path B from firing
  ! at small adaptive dt.
  !
  ! Note on the choice of L:
  !   The latent heat used here is Lvap_T(ts), evaluated at the SURFACE
  !   temperature — not at the layer temperature.  Rationale: every kg of
  !   vapor in the column was placed there by surface evaporation, which
  !   debited L(ts) from the slab.  Releasing L(tmid) at condensation would
  !   credit L_sub for any kg condensing at T<273.16 K while only L_vap was
  !   originally consumed at the warm surface, injecting a free L_fusion per
  !   kg cycled.  Using Lvap_T(ts) keeps the latent-heat ledger balanced
  !   (evap debit ≡ condense credit) and treats the column as if precip
  !   leaves the column in the same phase that it entered.  Phase-aware qsat
  !   (esat_cc uses L_sub below freezing → lower qsat over ice) is retained
  !   — that's a separate, correct piece of physics that affects WHEN
  !   condensation triggers, not how much heat is released per kg.
  !
  ! After the latent release, qsat(T_new) > qsat(T_old), so the layer may be
  ! subsaturated relative to its NEW temperature, and the latent heating can
  ! destabilize adjacent pairs.  The caller follows this routine with a
  ! convadj-only iteration pass to clean up.
  !
  ! Outputs (one call per outer step):
  !   precip_mass_flux  — column-integrated condensation mass flux this step
  !                       [kg/m²/s]; reset to 0 here.
  !   cond_heating(k)   — ADDED to (caller zeroes once before this call); the
  !                       inc here is the only contribution per outer step.
    real(r8), intent(in)  :: dt_step_sec      ! step length [s] for diagnostics
    real(r8), intent(out) :: precip_mass_flux

    integer  :: k
    real(r8) :: eps_wv, es_k, qsat_k, q_excess, L_release, alpha

    eps_wv    = SHR_CONST_MWWV / mwdry_col
    L_release = Lvap_T(ts)         ! match surface evap ledger
    ! Implicit-Euler relaxation factor: alpha = (dt/τ)/(1+dt/τ) = dt/(dt+τ)
    alpha     = dt_step_sec / (dt_step_sec + tau_cond)
    precip_mass_flux = 0._r8

    do k = 1, pver
      es_k   = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
      qsat_k = eps_wv * es_k / (pmid(k) - es_k)
      if (h2ommr(k) > qsat_k) then
        q_excess  = alpha * (h2ommr(k) - qsat_k)
        h2ommr(k) = h2ommr(k) - q_excess
        tmid(k)   = tmid(k) + L_release * q_excess / cpdry_col
        ! Diagnostic: latent heating [K/day], column precip mass [kg/m²/s].
        cond_heating(k)  = cond_heating(k) &
                           + L_release * q_excess / cpdry_col / dt_days
        precip_mass_flux = precip_mass_flux + q_excess * pdel(k) / exo_g / dt_step_sec
      end if
    end do

  end subroutine condense

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
    use exocol_mod, only: zint, tmid, pint, mwdry_col
    integer  :: k
    real(r8) :: R_gas
    R_gas = SHR_CONST_RGAS / mwdry_col
    do k = pver, 1, -1
      zint(k) = zint(k+1) + (R_gas / exo_g) * tmid(k) * log(pint(k+1) / pint(k))
    end do
  end subroutine update_zint

end module exocol_rce_loop
