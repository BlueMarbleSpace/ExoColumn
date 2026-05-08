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
  use shr_const_mod,   only: SHR_CONST_CSEC
  use exoplanet_mod,   only: exo_g
  use ppgrid,          only: pver, pverp
  use exocol_mod
  use exocol_radiation, only: exocol_rad_tend
  use exocol_convadj,   only: convadj_dry

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

    real(r8), dimension(pver) :: tmid_snap   ! tmid at last stability snapshot
    real(r8) :: ts_snap                      ! Ts at last stability snapshot
    real(r8) :: dt_sec        ! dt in seconds
    real(r8) :: max_hr        ! max total heating rate this step [K/day]
    real(r8) :: toa_flux      ! TOA net flux this step [W/m²]
    real(r8) :: F_net_srf     ! net surface flux [W/m²]
    real(r8) :: max_dT        ! max tmid change since last snapshot [K]
    real(r8) :: dTs           ! Ts change since last snapshot [K]

    integer  :: it, k
    logical  :: converged, profile_stable

    dt_sec         = dt_days * SHR_CONST_CSEC
    converged      = .false.
    profile_stable = .false.
    tmid_snap      = tmid
    ts_snap        = ts

    write(*,'(/,a)')    '========================================'
    write(*,'(a)')      ' ExoColumn RCE loop starting'
    write(*,'(a,f6.1)') '   dt [days]  = ', dt_days
    write(*,'(a,i0)')   '   nmax       = ', nmax
    write(*,'(a,f7.4)') '   hr_tol  [K/day] = ', hr_tol
    write(*,'(a,f6.3)') '   toa_tol [W/m2]  = ', toa_tol
    write(*,'(a,/)')    '========================================'

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

      ! 4. Recompute interface temperatures from updated tmid and ts.
      !    tint(pverp) = ts (surface);  inner interfaces interpolated in log-p.
      call update_tint()

      ! 5. Convective adjustment (dry adiabat, Phase 1)
      !    tint(pverp) is set to ts inside update_tint before this call.
      call convadj_dry(tmid, tint, pint, pdel, cpdry_col, exo_g, pver)

      ! 6. Diagnostics and convergence check
      max_hr   = maxval(abs(LWHR(:) + SWHR(:)))
      toa_flux = abs(SWDN(1) - SWUP(1) - LWUP(1))

      ! Profile-stability check every stab_check steps (Path B)
      if (mod(it, stab_check) == 0) then
        max_dT         = maxval(abs(tmid - tmid_snap))
        dTs            = abs(ts - ts_snap)
        profile_stable = (max_dT < prof_stab_tol .and. dTs < ts_stab_tol)
        tmid_snap      = tmid
        ts_snap        = ts
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

      ! Path B: profile stability with TOA balance
      if (profile_stable .and. toa_flux < toa_tol) then
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

    ! Final radiation call to return consistent flux diagnostics
    call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

  end subroutine run_rce_loop

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

end module exocol_rce_loop
