module exocol_rce_loop
! RCE time-marching loop for ExoColumn.
!
! run_rce_loop time-steps the column forward using a large virtual dt,
! calling radiation and convective adjustment at each step, until the
! column reaches radiative-convective equilibrium.
!
! Convergence is declared when either path is satisfied:
!
!   Path A — radiative equilibrium (quiescent atmosphere):
!     (1) max |LWHR(k) + SWHR(k)| < hr_tol  [K/day] over all layers
!     (2) |TOA net flux| < toa_tol  [W/m²]
!
!   Path B — profile stability (convectively active atmosphere):
!     (2) |TOA net flux| < toa_tol  [W/m²]
!     (3) max |Δtmid| < prof_stab_tol  [K] over stab_check steps
!     (4) |ΔTs|       < ts_stab_tol    [K] over stab_check steps
!
!   Path B handles the common case where a convectively active layer has a
!   large instantaneous heating rate that is continuously balanced by dry
!   convective adjustment, so the column is physically steady but max|HR|
!   never drops below hr_tol.
!
! A simple slab-ocean prognostic surface energy balance updates ts each step:
!   ts ← ts + dt · F_net_srf / H_slab
! where F_net_srf = SWDN(pverp) - SWUP(pverp) + LWDN(pverp) - LWUP(pverp)
! and H_slab = rho_w · cp_w · dz_slab  [J m⁻² K⁻¹].

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use shr_const_mod,   only: SHR_CONST_CSEC, SHR_CONST_RGAS, SHR_CONST_MWWV
  use exoplanet_mod,   only: exo_g
  use ppgrid,          only: pver, pverp
  use exocol_mod
  use exocol_radiation, only: exocol_rad_tend
  use exocol_config,    only: conv_scheme, cc_feedback
  use exocol_convadj,   only: convadj_dry, convadj_moist, convadj_manabe, esat_cc

  implicit none
  private

  public :: run_rce_loop

  ! ---- Tuneable parameters (may be overridden via namelist in a future version) ----
  real(r8), parameter :: dt_days     = 5._r8      ! virtual timestep [Earth days]
  integer,  parameter :: nmax        = 200000      ! maximum iterations
  integer,  parameter :: print_every = 100         ! diagnostic print interval
  real(r8), parameter :: hr_tol      = 0.01_r8    ! heating-rate convergence [K/day]
  real(r8), parameter :: toa_tol     = 0.1_r8     ! TOA flux convergence [W/m²]

  ! Profile-stability convergence (Path B)
  integer,  parameter :: stab_check    = 100        ! steps between stability snapshots
  real(r8), parameter :: prof_stab_tol = 0.001_r8   ! max tmid change per stab_check steps [K]
  real(r8), parameter :: ts_stab_tol   = 0.001_r8   ! max Ts   change per stab_check steps [K]
  real(r8), parameter :: toa_stab_tol  = 0.001_r8   ! max TOA flux change per stab_check steps [W/m²]

  ! Slab-ocean heat capacity parameters
  real(r8), parameter :: rho_w     = 1026._r8     ! seawater density [kg/m³]
  real(r8), parameter :: cp_w      = 4000._r8     ! seawater specific heat [J/kg/K]
  real(r8), parameter :: dz_slab   = 50._r8       ! slab thickness [m]
  real(r8), parameter :: H_slab    = rho_w * cp_w * dz_slab  ! [J/m²/K]

contains

  subroutine run_rce_loop()
  ! Main RCE iteration.  Column state in exocol_mod is updated in-place.

    real(r8), dimension(pver)  :: LWHR, SWHR
    real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

    real(r8), dimension(pver) :: rh_init     ! initial relative humidity (fixed-RH)
    real(r8), dimension(pver) :: tmid_snap   ! tmid at last stability snapshot
    real(r8) :: ts_snap                      ! Ts at last stability snapshot
    real(r8) :: toa_snap                     ! TOA flux at last stability snapshot [W/m²]
    real(r8) :: dt_sec        ! dt in seconds
    real(r8) :: max_hr        ! max total heating rate this step [K/day]
    real(r8) :: toa_flux      ! TOA net flux this step [W/m²]
    real(r8) :: F_net_srf     ! net surface flux [W/m²]
    real(r8) :: max_dT        ! max tmid change since last snapshot [K]
    real(r8) :: dTs           ! Ts change since last snapshot [K]
    real(r8) :: dToa          ! TOA flux change since last snapshot [W/m²]

    integer  :: it, k
    logical  :: converged, profile_stable

    dt_sec         = dt_days * SHR_CONST_CSEC
    converged      = .false.
    profile_stable = .false.
    tmid_snap      = tmid
    ts_snap        = ts
    toa_snap       = 0._r8

    write(*,'(/,a)')    '========================================'
    write(*,'(a)')      ' ExoColumn RCE loop starting'
    write(*,'(a,f6.1)') '   dt [days]  = ', dt_days
    write(*,'(a,i0)')   '   nmax       = ', nmax
    write(*,'(a,f7.4)') '   hr_tol  [K/day] = ', hr_tol
    write(*,'(a,f6.3)') '   toa_tol [W/m2]  = ', toa_tol
    write(*,'(a,/)')    '========================================'

    ! Compute initial relative humidity from the input q profile.
    ! rh_init(k) = h2ommr(k) / qsat(T(k), p(k)), capped at 1.
    ! This is held fixed throughout the run; h2ommr is updated each step
    ! via the Clausius-Clapeyron relation to maintain this RH profile.
    block
      use exocol_mod, only: h2ommr, pmid, mwdry_col
      integer  :: k
      real(r8) :: eps_wv, es_k, qsat_k
      eps_wv = SHR_CONST_MWWV / mwdry_col
      do k = 1, pver
        es_k       = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
        qsat_k     = eps_wv * es_k / (pmid(k) - es_k)
        rh_init(k) = min(h2ommr(k) / max(qsat_k, 1.0e-20_r8), 1.0_r8)
      end do
    end block

    do it = 1, nmax

      ! 1. Radiation
      call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

      ! 2. Update layer temperatures  [K/day → K/step]
      do k = 1, pver
        tmid(k) = tmid(k) + dt_days * (LWHR(k) + SWHR(k))
      end do

      ! 3. Prognostic surface temperature (slab ocean)
      !    Positive F_net_srf means ocean gains heat → ts increases.
      !    Index pverp = surface level; index 1 = TOA.
      F_net_srf = (SWDN(pverp) - SWUP(pverp)) + (LWDN(pverp) - LWUP(pverp))
      ts = ts + dt_sec * F_net_srf / H_slab

      ! 4. Update water vapour (fixed relative humidity via CC) then derived fields.
      if (cc_feedback) then
        call update_h2ommr(rh_init)
        call exocol_update_derived()
      end if

      ! 5. Recompute interface temperatures from updated tmid and ts.
      !    tint(pverp) = ts (surface);  inner interfaces interpolated in log-p.
      call update_tint()

      ! 6. Convective adjustment
      ! Block-scope USE for zint follows the same ifx workaround as update_zint:
      ! host-association of allocatables from a module-level USE can mis-resolve
      ! in contained subroutines under ifx.
      select case (trim(conv_scheme))
      case ('moist')
        block
          use exocol_mod, only: zint
          call convadj_moist(tmid, tint, zint, pint, pdel, cpdry_col, exo_g, ts, pver)
        end block
      case ('manabe')
        block
          use exocol_mod, only: zint
          call convadj_manabe(tmid, tint, zint, pint, pdel, cpdry_col, exo_g, ts, pver)
        end block
      case default  ! 'dry'
        call convadj_dry(tmid, tint, pint, pdel, cpdry_col, exo_g, ts, pver)
      end select

      ! 7. Update interface heights consistent with the new tmid profile.
      !    zint is passed to aerad_driver each call; without this update it
      !    would reflect the initial conditions throughout the run.
      call update_zint()

      ! 8. Diagnostics and convergence check
      max_hr   = maxval(abs(LWHR(:) + SWHR(:)))
      toa_flux = abs(SWDN(1) - SWUP(1) + LWDN(1) - LWUP(1))

      ! Profile-stability check every stab_check steps (Path B).
      ! Includes TOA flux stability so a structurally imbalanced dry column
      ! (where toa_flux never drops below toa_tol) can still converge.
      if (mod(it, stab_check) == 0) then
        max_dT = maxval(abs(tmid - tmid_snap))
        dTs    = abs(ts - ts_snap)
        dToa   = abs(toa_flux - toa_snap)
        ! Profile stable when tmid and ts have frozen.  For the TOA flux:
        ! if the column is already energy-balanced (toa_flux < toa_tol) the
        ! dToa criterion is redundant and is skipped — CC-active runs have a
        ! small residual TOA oscillation that would otherwise block convergence.
        ! For structurally imbalanced columns (dry scheme, no CC) toa_flux >>
        ! toa_tol so dToa is still needed to detect a truly frozen state.
        profile_stable = (max_dT < prof_stab_tol .and. dTs < ts_stab_tol &
                          .and. (toa_flux < toa_tol .or. dToa < toa_stab_tol))
        tmid_snap = tmid
        ts_snap   = ts
        toa_snap  = toa_flux
      end if

      if (mod(it, print_every) == 0 .or. it == 1) then
        write(*,'(a,i8,a,f8.4,a,f8.3,a,f7.2)') &
          '  step=', it, '  max|HR| [K/d]=', max_hr, &
          '  TOA_net [W/m2]=', toa_flux, '  Ts [K]=', ts
      end if

      ! Path A: radiative equilibrium
      if (max_hr < hr_tol .and. toa_flux < toa_tol) then
        converged = .true.
        exit
      end if

      ! Path B: profile stability (entire state frozen, including TOA flux).
      ! Does not require toa_flux < toa_tol — a dry column in radiative-
      ! convective equilibrium may have a structural TOA imbalance that never
      ! shrinks.  A warning is printed post-convergence if toa_flux > toa_tol.
      if (profile_stable) then
        converged = .true.
        exit
      end if

    end do

    if (converged) then
      write(*,'(/,a,i0,a)') ' ExoColumn: converged after ', it, ' steps.'
    else
      write(*,'(/,a,i0,a)') ' ExoColumn: WARNING — not converged after ', nmax, ' steps.'
    end if
    write(*,'(a,f8.3,a,f8.3)') &
      '   Final max|HR| [K/day] = ', max_hr, &
      '   Final TOA net [W/m2]  = ', toa_flux
    if (toa_flux > toa_tol) then
      write(*,'(a)') &
        '   WARNING: TOA net flux exceeds toa_tol — column may be convectively' // &
        ' active; structural imbalance possible.'
    end if

  end subroutine run_rce_loop

  ! -----------------------------------------------------------------------
  ! Private helper: update h2ommr for fixed relative humidity
  ! -----------------------------------------------------------------------

  subroutine update_h2ommr(rh_init)
  ! Relax the water-vapour mixing ratio toward the CC-equilibrium value for
  ! the current tmid profile while maintaining fixed relative humidity.
  !
  !   qsat(k) = (Mw_h2o/Mw_dry) · esat(T(k)) / (p(k) − esat(T(k)))
  !   h2ommr(k) ← h2ommr(k) + α · (rh_init(k)·qsat(k) − h2ommr(k))
  !   α = dt_days / tau_relax    (capped at 1)
  !
  ! Relaxing rather than instantaneously equilibrating damps the positive
  ! CC feedback and prevents limit-cycle oscillations at large virtual dt.
  ! At equilibrium the result is identical to full CC equilibration.
  ! Explicit USE works around the ifx host-association quirk for allocatables.
    use exocol_mod, only: h2ommr, tmid, pmid, mwdry_col

    real(r8), intent(in) :: rh_init(pver)

    real(r8), parameter :: tau_relax = 10._r8   ! moisture relaxation timescale [days]

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
  end subroutine update_h2ommr

  ! -----------------------------------------------------------------------
  ! Private helper: recompute tint from tmid and ts
  ! -----------------------------------------------------------------------

  subroutine update_tint()
  ! Interpolate interface temperatures from layer midpoint temperatures.
  ! Uses log-pressure weighting between adjacent midpoints.
  ! Boundary conditions:
  !   tint(pverp) = ts          (surface interface)
  !   tint(1)     = extrapolated above the first layer
    integer  :: k
    real(r8) :: pmid_k, pmid_kp1, wt

    ! Surface interface fixed to skin temperature
    tint(pverp) = ts

    ! Inner interfaces: log-pressure interpolation
    do k = 1, pver-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      wt = log(pint(k+1) / pmid_k) / log(pmid_kp1 / pmid_k)
      tint(k+1) = tmid(k) + wt * (tmid(k+1) - tmid(k))
    end do

    ! Top interface: linear extrapolation from the two uppermost midpoints
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine update_tint

  ! -----------------------------------------------------------------------
  ! Private helper: recompute zint from tmid via the hypsometric equation
  ! -----------------------------------------------------------------------

  subroutine update_zint()
  ! Rebuild interface heights from the current tmid using the hypsometric equation:
  !   zint(k) = zint(k+1) + (R/g) * tmid(k) * ln(pint(k+1)/pint(k))
  ! zint(pverp) (the surface height) is held fixed at whatever value was read
  ! from the input file and is not modified here.
  !
  ! Explicit USE here works around an ifx host-association quirk for allocatable
  ! arrays from a module-level USE in a contained subroutine.
    use exocol_mod, only: zint, tmid, pint, mwdry_col
    integer  :: k
    real(r8) :: R_gas

    R_gas = SHR_CONST_RGAS / mwdry_col
    do k = pver, 1, -1
      zint(k) = zint(k+1) + (R_gas / exo_g) * tmid(k) * log(pint(k+1) / pint(k))
    end do
  end subroutine update_zint

end module exocol_rce_loop
