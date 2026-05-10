module exocol_convadj
! Convective adjustment schemes for ExoColumn.
!
! Three schemes are available; select via exocol_config::conv_scheme:
!
!   convadj_dry     Dry adiabatic.  Potential-temperature (θ) criterion;
!                   adjusts unstable pairs to equal θ while conserving
!                   cp·T·Δp.  (Original Phase 1 implementation.)
!
!   convadj_moist   Moist pseudo-adiabatic.  Like Manabe-Wetherald but with a
!                   dynamically computed local Γm(T,p) instead of a fixed 6.5 K/km.
!                   Γm = (g/cp)·(1+Lv·ws/(Rd·T))/(1+Lv²·ws/(cp·Rv·T²)).
!                   Requires interface heights; moisture is not modified.
!
!   convadj_manabe  Manabe-Wetherald fixed lapse rate.  Adjusts any layer pair
!                   whose geometric temperature gradient exceeds γ_crit = 6.5 K/km,
!                   restoring exactly γ_crit while conserving cp·T·Δp.  No
!                   moisture update.  Requires interface heights (zint_if) as an
!                   argument; uses the previous-timestep heights (sufficient since
!                   zint is updated by the caller after this routine returns).

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_RGAS, SHR_CONST_LATVAP, &
                            SHR_CONST_TKFRZ, SHR_CONST_RWV
  use ppgrid,        only: pver, pverp

  implicit none
  private

  public :: convadj_dry
  public :: convadj_moist
  public :: convadj_manabe
  public :: esat_cc

  ! Manabe-Wetherald critical lapse rate [K/m]
  real(r8), parameter :: gamma_crit = 6.5e-3_r8

  ! Saturation vapour pressure reference values for esat_cc / malr
  real(r8), parameter :: es0    = 611.2_r8          ! esat at T0_sat [Pa]
  real(r8), parameter :: T0_sat = SHR_CONST_TKFRZ   ! reference temperature [K]

contains

  ! -----------------------------------------------------------------------
  ! Scheme 0 — Dry adiabatic adjustment
  ! -----------------------------------------------------------------------

  subroutine convadj_dry(tmid, tint, pint, pdel, cp, g, ts, nv)
  ! Dry adiabatic convective adjustment.
  !
  ! Arguments:
  !   tmid(nv)     IN/OUT  layer-midpoint temperatures [K]
  !   tint(nv+1)   IN/OUT  interface temperatures [K]; tint(nv+1) (surface)
  !                        must be set to ts by the caller before this call
  !                        and is not modified here.
  !   pint(nv+1)   IN      interface pressures [Pa]; index 1 = TOA, nv+1 = srf
  !   pdel(nv)     IN      layer pressure thicknesses [Pa]
  !   cp           IN      specific heat of dry air [J/kg/K]
  !   g            IN      gravitational acceleration [m/s²] (reserved; unused)
  !   ts           IN      surface skin temperature [K]
  !   nv           IN      number of layers (normally = pver)
  !
  ! Stability criterion: θ(k) = T(k)·(p_ref/p(k))^κ must be non-decreasing
  ! upward.  Unstable: θ(k) < θ(k+1) with k above k+1.
  !
  ! Adjustment: conserve cp·T·Δp of the pair and restore θ(k) = θ(k+1):
  !   T'_{k+1} = H / (Δp_k·r + Δp_{k+1})
  !   T'_k     = r · T'_{k+1}
  ! where H = Δp_k·T_k + Δp_{k+1}·T_{k+1},  r = (pmid_k/pmid_{k+1})^κ.

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g
    real(r8), intent(in)    :: ts
    integer,  intent(in)    :: nv

    integer, parameter :: max_pass = 30

    integer  :: k, ipass
    logical  :: adjusted
    real(r8) :: pmid_k, pmid_kp1
    real(r8) :: kappa
    real(r8) :: ratio_kappa
    real(r8) :: theta_k_scaled
    real(r8) :: H
    real(r8) :: Tkp1_new, Tk_new

    block
      use exocol_mod, only: mwdry_col
      kappa = (SHR_CONST_RGAS / mwdry_col) / cp
    end block

    do ipass = 1, max_pass
      adjusted = .false.

      ! Surface–bottom-layer pair
      block
        real(r8) :: pmid_bot, ratio_bot
        pmid_bot  = 0.5_r8 * (pint(nv) + pint(nv+1))
        ratio_bot = (pmid_bot / pint(nv+1))**kappa
        if (tmid(nv) / ratio_bot < ts) then
          tmid(nv) = ratio_bot * ts
          adjusted = .true.
        end if
      end block

      ! Atmospheric sweep, surface → TOA
      do k = nv-1, 1, -1
        pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
        pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))

        ratio_kappa    = (pmid_k / pmid_kp1)**kappa
        theta_k_scaled = tmid(k) / ratio_kappa

        if (theta_k_scaled < tmid(k+1)) then
          H         = pdel(k)*tmid(k) + pdel(k+1)*tmid(k+1)
          Tkp1_new  = H / (pdel(k)*ratio_kappa + pdel(k+1))
          Tk_new    = ratio_kappa * Tkp1_new
          tmid(k)   = Tk_new
          tmid(k+1) = Tkp1_new
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interface temperatures from the adjusted tmid profile.
    ! tint(nv+1) = ts is fixed by the caller; do not touch it.
    do k = 1, nv-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      tint(k+1) = tmid(k) + (tmid(k+1) - tmid(k)) * &
                  log(pint(k+1)/pmid_k) / log(pmid_kp1/pmid_k)
    end do
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine convadj_dry

  ! -----------------------------------------------------------------------
  ! Scheme 1 — Moist pseudo-adiabatic adjustment
  ! -----------------------------------------------------------------------

  subroutine convadj_moist(tmid, tint, zint_if, pint, pdel, cp, g, ts, nv)
  ! Moist pseudo-adiabatic convective adjustment.
  !
  ! Stability criterion: the geometric lapse rate between midpoints of two
  ! adjacent layers must not exceed the local moist adiabatic lapse rate Γm:
  !   (T_{k+1} − T_k) / (zmid_k − zmid_{k+1}) > Γm(T̄, p̄)  →  unstable
  !
  ! Γm is evaluated at the mean T and p of the pair via the malr() helper:
  !   Γm = (g/cp) · (1 + Lv·ws/(Rd·T)) / (1 + Lv²·ws/(cp·Rv·T²))
  ! where ws = eps·esat(T̄)/(p̄ − esat(T̄)) is the saturation mixing ratio
  ! and esat follows the Clausius-Clapeyron relation.
  !
  ! Adjustment: restore Γm exactly while conserving cp·T·Δp:
  !   T'_{k+1} = (H + Γm·Δz·Δp_k) / (Δp_k + Δp_{k+1})
  !   T'_k     = T'_{k+1} − Γm·Δz
  !
  ! Moisture is not modified.  Γm is floored at 1 K/km.

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(in)    :: zint_if(nv+1)  ! interface heights [m]
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g
    real(r8), intent(in)    :: ts
    integer,  intent(in)    :: nv

    integer, parameter :: max_pass = 30
    real(r8), parameter :: dz_min = 1._r8

    integer  :: k, ipass
    logical  :: adjusted
    real(r8) :: Rd
    real(r8) :: zmid_k, zmid_kp1, dz
    real(r8) :: T_mean, p_mean, gamma_m
    real(r8) :: gamma_actual
    real(r8) :: H, Tkp1_new, Tk_new
    real(r8) :: pmid_k, pmid_kp1

    block
      use exocol_mod, only: mwdry_col
      Rd = SHR_CONST_RGAS / mwdry_col
    end block

    do ipass = 1, max_pass
      adjusted = .false.

      ! Surface–bottom-layer pair
      block
        real(r8) :: zmid_bot, dz_surf, pmid_bot
        zmid_bot = 0.5_r8 * (zint_if(nv) + zint_if(nv+1))
        dz_surf  = zmid_bot - zint_if(nv+1)
        if (dz_surf > dz_min) then
          pmid_bot = 0.5_r8 * (pint(nv) + pint(nv+1))
          T_mean   = 0.5_r8 * (tmid(nv) + ts)
          p_mean   = 0.5_r8 * (pmid_bot + pint(nv+1))
          gamma_m  = malr(T_mean, p_mean, Rd, g, cp)
          if ((ts - tmid(nv)) / dz_surf > gamma_m) then
            tmid(nv) = ts - gamma_m * dz_surf
            adjusted = .true.
          end if
        end if
      end block

      ! Atmospheric sweep, surface → TOA
      do k = nv-1, 1, -1
        zmid_k   = 0.5_r8 * (zint_if(k)   + zint_if(k+1))
        zmid_kp1 = 0.5_r8 * (zint_if(k+1) + zint_if(k+2))
        dz       = zmid_k - zmid_kp1

        if (dz < dz_min) cycle

        gamma_actual = (tmid(k+1) - tmid(k)) / dz

        T_mean  = 0.5_r8 * (tmid(k) + tmid(k+1))
        p_mean  = pint(k+1)
        gamma_m = malr(T_mean, p_mean, Rd, g, cp)

        if (gamma_actual > gamma_m) then
          H        = pdel(k)*tmid(k) + pdel(k+1)*tmid(k+1)
          Tkp1_new = (H + gamma_m*dz*pdel(k)) / (pdel(k) + pdel(k+1))
          Tk_new   = Tkp1_new - gamma_m * dz
          tmid(k)   = Tk_new
          tmid(k+1) = Tkp1_new
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interface temperatures
    do k = 1, nv-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      tint(k+1) = tmid(k) + (tmid(k+1) - tmid(k)) * &
                  log(pint(k+1)/pmid_k) / log(pmid_kp1/pmid_k)
    end do
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine convadj_moist

  ! -----------------------------------------------------------------------
  ! Scheme 2 — Manabe-Wetherald fixed lapse rate
  ! -----------------------------------------------------------------------

  subroutine convadj_manabe(tmid, tint, zint_if, pint, pdel, cp, g, ts, nv)
  ! Manabe-Wetherald (1967) convective adjustment.
  !
  ! Stability criterion: the geometric temperature lapse rate between the
  ! midpoints of two adjacent layers must not exceed γ_crit = 6.5 K/km:
  !   (T_{k+1} − T_k) / (zmid_k − zmid_{k+1}) > γ_crit  →  unstable
  ! where zmid_k = 0.5·(zint_if(k) + zint_if(k+1)) and k is the upper layer.
  !
  ! Adjustment: restore γ_crit exactly while conserving cp·T·Δp:
  !   T'_{k+1} = (H + γ_crit·Δz·Δp_k) / (Δp_k + Δp_{k+1})
  !   T'_k     = T'_{k+1} − γ_crit·Δz
  ! where H = Δp_k·T_k + Δp_{k+1}·T_{k+1},  Δz = zmid_k − zmid_{k+1}.
  !
  ! Moisture is not updated.  zint_if reflects heights from the previous
  ! timestep (the caller updates zint after convective adjustment).

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(in)    :: zint_if(nv+1)  ! interface heights [m]
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g              ! reserved; unused
    real(r8), intent(in)    :: ts
    integer,  intent(in)    :: nv

    integer, parameter :: max_pass = 30
    real(r8), parameter :: dz_min = 1._r8   ! [m] skip pairs thinner than this

    integer  :: k, ipass
    logical  :: adjusted
    real(r8) :: zmid_k, zmid_kp1, dz
    real(r8) :: gamma_actual
    real(r8) :: H, Tkp1_new, Tk_new
    real(r8) :: pmid_k, pmid_kp1

    do ipass = 1, max_pass
      adjusted = .false.

      ! Surface–bottom-layer pair.
      ! zmid(nv) is the midpoint of the bottom layer; zint_if(nv+1) is the surface.
      block
        real(r8) :: zmid_bot, dz_surf
        zmid_bot = 0.5_r8 * (zint_if(nv) + zint_if(nv+1))
        dz_surf  = zmid_bot - zint_if(nv+1)
        if (dz_surf > dz_min) then
          if ((ts - tmid(nv)) / dz_surf > gamma_crit) then
            tmid(nv) = ts - gamma_crit * dz_surf
            adjusted = .true.
          end if
        end if
      end block

      ! Atmospheric sweep, surface → TOA
      do k = nv-1, 1, -1
        zmid_k   = 0.5_r8 * (zint_if(k)   + zint_if(k+1))
        zmid_kp1 = 0.5_r8 * (zint_if(k+1) + zint_if(k+2))
        dz       = zmid_k - zmid_kp1     ! positive: k is above k+1

        if (dz < dz_min) cycle

        gamma_actual = (tmid(k+1) - tmid(k)) / dz

        if (gamma_actual > gamma_crit) then
          H         = pdel(k)*tmid(k) + pdel(k+1)*tmid(k+1)
          Tkp1_new  = (H + gamma_crit*dz*pdel(k)) / (pdel(k) + pdel(k+1))
          Tk_new    = Tkp1_new - gamma_crit * dz
          tmid(k)   = Tk_new
          tmid(k+1) = Tkp1_new
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interface temperatures (same as dry scheme).
    do k = 1, nv-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      tint(k+1) = tmid(k) + (tmid(k+1) - tmid(k)) * &
                  log(pint(k+1)/pmid_k) / log(pmid_kp1/pmid_k)
    end do
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine convadj_manabe

  ! -----------------------------------------------------------------------
  ! Private helpers for moist adiabatic lapse rate
  ! -----------------------------------------------------------------------

  pure function esat_cc(T) result(es)
  ! Clausius-Clapeyron saturation vapour pressure [Pa] relative to liquid water.
    real(r8), intent(in) :: T   ! temperature [K]
    real(r8) :: es
    es = es0 * exp((SHR_CONST_LATVAP / SHR_CONST_RWV) * (1._r8/T0_sat - 1._r8/T))
  end function esat_cc

  pure function malr(T, p, Rd, g_planet, cp_dry) result(gamma_m)
  ! Moist adiabatic lapse rate [K/m] for given T [K], p [Pa], and gas constants.
  ! Floored at 1 K/km (= 0.001 K/m) to avoid pathological values where esat → 0.
    real(r8), intent(in) :: T, p, Rd, g_planet, cp_dry
    real(r8) :: gamma_m
    real(r8) :: eps, es, ws
    real(r8), parameter :: ws_tiny = 1.0e-10_r8
    real(r8), parameter :: gamma_floor = 1.0e-3_r8  ! 1 K/km

    eps = Rd / SHR_CONST_RWV
    es  = esat_cc(T)
    if (p > es) then
      ws = eps * es / (p - es)
    else
      ws = ws_tiny
    end if

    gamma_m = (g_planet / cp_dry) * &
              (1._r8 + SHR_CONST_LATVAP * ws / (Rd * T)) / &
              (1._r8 + SHR_CONST_LATVAP**2 * ws / (cp_dry * SHR_CONST_RWV * T**2))
    gamma_m = max(gamma_m, gamma_floor)
  end function malr

end module exocol_convadj
