module exocol_rce_loop
! RCE time-marching loop for ExoColumn.
!
! run_rce_loop time-steps the column forward using a large virtual dt and
! the following physics per step, in order:
!
!   1. Radiation tendency on tmid (ExoRT aerad_driver).
!   2. Bulk-aerodynamic surface fluxes LE, SH at the slab-atmosphere interface.
!   3. Slab-ocean energy balance with turbulent loss:
!        ts ← ts + dt · (F_net_srf_rad − LE − SH) / H_slab
!   4. Surface fluxes deposited into the bottom layer:
!        tmid(pver)   ← tmid(pver)   + dt · SH / (cp · pdel/g)
!        h2ommr(pver) ← h2ommr(pver) + dt · LE / (L  · pdel/g)
!   5. Convective adjustment (T + q mixed conservatively in adjusted pairs).
!   6. Condensation/precipitation: where q > qsat, cap and release latent heat
!      to the layer (phase-aware via Lvap_T).
!   7. Update derived (pdeldry, pintdry), tint, zint.
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
!     max |Δtmid| < prof_stab_tol  AND  |ΔTs| < ts_stab_tol  AND
!     ( |TOA net flux| < toa_tol  OR  |ΔTOA| < toa_stab_tol )
!
!   Path B handles cases where convective adjustment continuously balances a
!   large instantaneous radiative tendency.

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use shr_const_mod,   only: SHR_CONST_CSEC, SHR_CONST_RGAS, SHR_CONST_MWWV
  use exoplanet_mod,   only: exo_g
  use ppgrid,          only: pver, pverp
  use exocol_mod
  use exocol_radiation, only: exocol_rad_tend
  use exocol_config,    only: conv_scheme, moisture_scheme, wind_speed, C_D
  use exocol_convadj,   only: convadj_dry, convadj_moist, convadj_manabe, &
                              esat_cc, Lvap_T
  use exocol_surface,   only: compute_surface_fluxes

  implicit none
  private

  public :: run_rce_loop

  ! ---- Tuneable parameters ----
  real(r8), parameter :: dt_days     = 5._r8      ! virtual timestep [Earth days]
  integer,  parameter :: nmax        = 200000     ! maximum iterations
  integer,  parameter :: print_every = 100        ! diagnostic print interval
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

  ! Fixed-RH legacy closure (only used when moisture_scheme='fixed_rh')
  real(r8), parameter :: tau_relax = 50._r8       ! moisture relaxation [days]

contains

  subroutine run_rce_loop()
  ! Main RCE iteration.  Column state in exocol_mod is updated in-place.

    real(r8), dimension(pver)  :: LWHR, SWHR
    real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

    real(r8), dimension(pver) :: rh_init     ! fixed-RH legacy target
    real(r8), dimension(pver) :: tmid_snap   ! tmid at last stability snapshot
    real(r8) :: ts_snap                      ! Ts at last stability snapshot
    real(r8) :: toa_snap                     ! TOA flux at last stability snapshot
    real(r8) :: dt_sec
    real(r8) :: max_hr, toa_flux
    real(r8) :: F_net_srf_rad, F_srf_total
    real(r8) :: LE, SH
    real(r8) :: max_dT, dTs, dToa
    real(r8) :: precip_total                 ! column precip mass flux [kg/m²/s]

    integer  :: it, k
    logical  :: converged, profile_stable
    logical  :: prognostic, fixed_rh

    dt_sec         = dt_days * SHR_CONST_CSEC
    converged      = .false.
    profile_stable = .false.
    tmid_snap      = tmid
    ts_snap        = ts
    toa_snap       = 0._r8
    LE = 0._r8;  SH = 0._r8

    prognostic = (trim(adjustl(moisture_scheme)) == 'prognostic')
    fixed_rh   = (trim(adjustl(moisture_scheme)) == 'fixed_rh')

    write(*,'(/,a)')    '========================================'
    write(*,'(a)')      ' ExoColumn RCE loop starting'
    write(*,'(a,f6.1)') '   dt [days]  = ', dt_days
    write(*,'(a,i0)')   '   nmax       = ', nmax
    write(*,'(a,f7.4)') '   hr_tol  [K/day] = ', hr_tol
    write(*,'(a,f6.3)') '   toa_tol [W/m2]  = ', toa_tol
    write(*,'(a,/)')    '========================================'

    ! Capture initial RH only if the legacy fixed-RH closure is selected.
    if (fixed_rh) then
      call capture_rh_init(rh_init)
    end if

    do it = 1, nmax

      ! ---- 1. Radiation tendency ----
      call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

      do k = 1, pver
        tmid(k) = tmid(k) + dt_days * (LWHR(k) + SWHR(k))
      end do

      ! ---- 2. Surface turbulent fluxes (bulk aerodynamic, implicit-damped) ----
      ! The raw explicit fluxes are valid only when dt << surface-layer
      ! relaxation time τ = (pdel/g) / (ρ_air · C_D · U).  At the large virtual
      ! dt used here (5 days vs τ ~ hours), an explicit Euler step would let
      ! the bottom layer overshoot ts by orders of magnitude in one step.
      ! Damping by 1/(1 + dt/τ) is the implicit-Euler equivalent: it reduces
      ! to the explicit formula for dt << τ and to a full one-step
      ! equilibration of the bottom layer to ts / qsat(ts) for dt >> τ.
      ! Equilibrium is unchanged — only the per-step rate is rate-limited.
      call compute_surface_fluxes(ts, tmid(pver), h2ommr(pver), pmid(pver), &
                                  mwdry_col, cpdry_col, wind_speed, C_D, LE, SH)
      block
        real(r8) :: Rd, rho_air_bot, layer_mass_bot, tau_surf, damping
        Rd             = SHR_CONST_RGAS / mwdry_col
        rho_air_bot    = pmid(pver) / (Rd * tmid(pver))
        layer_mass_bot = pdel(pver) / exo_g
        tau_surf       = layer_mass_bot / (rho_air_bot * C_D * wind_speed)
        damping        = 1._r8 / (1._r8 + dt_sec / tau_surf)
        LE = LE * damping
        SH = SH * damping
      end block

      ! ---- 3. Slab-ocean energy balance ----
      F_net_srf_rad = (SWDN(pverp) - SWUP(pverp)) + (LWDN(pverp) - LWUP(pverp))
      if (prognostic) then
        F_srf_total = F_net_srf_rad - LE - SH
      else
        ! Legacy slab budget: radiation only.  LE/SH still computed for diagnostics.
        F_srf_total = F_net_srf_rad
      end if
      ts = ts + dt_sec * F_srf_total / H_slab

      ! ---- 4. Bottom-layer source from surface turbulent fluxes ----
      ! With damped LE/SH, the per-step layer change is bounded by the heat
      ! capacity / mass of the bottom layer (no overshoot regardless of dt).
      if (prognostic) then
        tmid(pver)   = tmid(pver) + dt_sec * SH / (cpdry_col * pdel(pver) / exo_g)
        h2ommr(pver) = h2ommr(pver) + &
                       dt_sec * (LE / Lvap_T(ts)) / (pdel(pver) / exo_g)
        h2ommr(pver) = max(h2ommr(pver), 0._r8)
      end if

      ! ---- 5. Legacy fixed-RH closure ----
      if (fixed_rh) then
        call update_h2ommr_fixed_rh(rh_init)
      end if

      ! ---- 6. Update derived (q has changed in either branch) ----
      call exocol_update_derived()

      ! ---- 7. Interface temperatures ----
      call update_tint()

      ! ---- 8. Convective adjustment (T + q mixed where adjusted) ----
      select case (trim(adjustl(conv_scheme)))
      case ('moist')
        block
          use exocol_mod, only: zint, h2ommr
          call convadj_moist(tmid, tint, h2ommr, zint, pint, pdel, &
                             cpdry_col, exo_g, ts, pver)
        end block
      case ('manabe')
        block
          use exocol_mod, only: zint, h2ommr
          call convadj_manabe(tmid, tint, h2ommr, zint, pint, pdel, &
                              cpdry_col, exo_g, ts, pver)
        end block
      case default  ! 'dry'
        block
          use exocol_mod, only: h2ommr
          call convadj_dry(tmid, tint, h2ommr, pint, pdel, &
                           cpdry_col, exo_g, ts, pver)
        end block
      end select

      ! ---- 9. Condensation / precipitation ----
      precip_total = 0._r8
      cond_heating = 0._r8
      if (prognostic) then
        call condense(dt_sec, precip_total)
      end if

      ! ---- 10. Final derived update (q changed in convadj + condensation) ----
      call exocol_update_derived()

      ! ---- 11. Heights ----
      call update_zint()

      ! ---- Diagnostics for output ----
      LE_diag     = LE
      SH_diag     = SH
      precip_diag = precip_total * 86400._r8   ! kg/m²/s → mm/day (ρ_water = 1000 kg/m³)

      ! ---- Convergence ----
      max_hr   = maxval(abs(LWHR(:) + SWHR(:)))
      toa_flux = abs(SWDN(1) - SWUP(1) + LWDN(1) - LWUP(1))

      if (mod(it, stab_check) == 0) then
        max_dT = maxval(abs(tmid - tmid_snap))
        dTs    = abs(ts - ts_snap)
        dToa   = abs(toa_flux - toa_snap)
        profile_stable = (max_dT < prof_stab_tol .and. dTs < ts_stab_tol &
                          .and. (toa_flux < toa_tol .or. dToa < toa_stab_tol))
        tmid_snap = tmid
        ts_snap   = ts
        toa_snap  = toa_flux
      end if

      if (mod(it, print_every) == 0 .or. it == 1) then
        write(*,'(a,i7,a,f7.3,a,f8.3,a,f7.2,a,f6.1,a,f6.1,a,f6.2)') &
          '  step=', it, &
          '  max|HR|=', max_hr, &
          '  TOA=', toa_flux, &
          '  Ts=', ts, &
          '  LE=', LE, &
          '  SH=', SH, &
          '  P[mm/d]=', precip_diag
      end if

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
      write(*,'(/,a,i0,a)') ' ExoColumn: converged after ', it, ' steps.'
    else
      write(*,'(/,a,i0,a)') ' ExoColumn: WARNING — not converged after ', nmax, ' steps.'
    end if
    write(*,'(a,f8.3,a,f8.3)') &
      '   Final max|HR| [K/day] = ', max_hr, &
      '   Final TOA net [W/m2]  = ', toa_flux
    write(*,'(a,f8.3,a,f8.3,a,f7.3)') &
      '   Final LE [W/m2]       = ', LE_diag, &
      '   SH [W/m2] = ', SH_diag, &
      '   P [mm/day] = ', precip_diag
    if (toa_flux > toa_tol) then
      write(*,'(a)') &
        '   WARNING: TOA net flux exceeds toa_tol — column may be convectively' // &
        ' active; structural imbalance possible.'
    end if

  end subroutine run_rce_loop

  ! -----------------------------------------------------------------------
  ! Condensation / precipitation with phase-aware latent-heat release
  ! -----------------------------------------------------------------------

  subroutine condense(dt_step_sec, precip_mass_flux)
  ! Where h2ommr > qsat(T,p), cap q at qsat and release latent heat L(T) per
  ! kg condensed into the layer's tmid.  L = L_v above 0 °C, L_s below.
  !
  ! After the latent release, qsat(T_new) > qsat(T_old) = q_new, so the layer
  ! is subsaturated and no iteration is needed.
  !
  ! precip_mass_flux is the column-integrated condensed mass per unit time
  ! [kg/m²/s], reported back for the precip diagnostic.
    real(r8), intent(in)  :: dt_step_sec      ! step length [s] for diagnostics
    real(r8), intent(out) :: precip_mass_flux

    integer  :: k
    real(r8) :: eps_wv, es_k, qsat_k, q_excess, L_k

    eps_wv = SHR_CONST_MWWV / mwdry_col
    precip_mass_flux = 0._r8

    do k = 1, pver
      es_k   = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
      qsat_k = eps_wv * es_k / (pmid(k) - es_k)
      if (h2ommr(k) > qsat_k) then
        q_excess = h2ommr(k) - qsat_k
        L_k      = Lvap_T(tmid(k))
        h2ommr(k) = qsat_k
        tmid(k)   = tmid(k) + L_k * q_excess / cpdry_col
        ! Diagnostic: latent heating in K/day, and column-integrated precip mass.
        cond_heating(k)  = L_k * q_excess / cpdry_col / dt_days
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
    integer  :: k
    real(r8) :: pmid_k, pmid_kp1, wt

    tint(pverp) = ts
    do k = 1, pver-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      wt = log(pint(k+1) / pmid_k) / log(pmid_kp1 / pmid_k)
      tint(k+1) = tmid(k) + wt * (tmid(k+1) - tmid(k))
    end do
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)
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
