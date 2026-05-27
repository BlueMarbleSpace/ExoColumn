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
                            SHR_CONST_LATSUB, &
                            SHR_CONST_TKFRZ, SHR_CONST_RWV, &
                            SHR_CONST_MWWV
  use ppgrid,        only: pver, pverp

  implicit none
  private

  public :: convadj_dry
  public :: convadj_moist
  public :: convadj_manabe
  public :: esat_cc
  public :: Lvap_T
  public :: malr
  public :: compute_tint_interp

  ! Manabe-Wetherald critical lapse rate [K/m]
  real(r8), parameter :: gamma_crit = 6.5e-3_r8

  ! Saturation vapour pressure reference values for esat_cc / malr
  real(r8), parameter :: es0    = 611.2_r8          ! esat at T0_sat [Pa]
  real(r8), parameter :: T0_sat = SHR_CONST_TKFRZ   ! reference temperature [K]

contains

  ! -----------------------------------------------------------------------
  ! Scheme 0 — Dry adiabatic adjustment
  ! -----------------------------------------------------------------------

  subroutine convadj_dry(tmid, tint, h2ommr_col, pint, pdel, cp, g, ts, nv)
  ! Dry adiabatic convective adjustment.
  !
  ! Arguments:
  !   tmid(nv)        IN/OUT  layer-midpoint temperatures [K]
  !   tint(nv+1)      IN/OUT  interface temperatures [K]; tint(nv+1) (surface)
  !                           must be set to ts by the caller before this call
  !                           and is not modified here.
  !   h2ommr_col(nv)  IN/OUT  specific humidity [kg/kg]; mixed in every pair
  !                           that is convectively adjusted (mass-weighted by
  !                           pdel, conservative).  Surface-bottom pair does
  !                           not mix q (slab has no q field).
  !   pint(nv+1)      IN      interface pressures [Pa]; index 1 = TOA, nv+1 = srf
  !   pdel(nv)        IN      layer pressure thicknesses [Pa]
  !   cp              IN      specific heat of dry air [J/kg/K]
  !   g               IN      gravitational acceleration [m/s²] (reserved; unused)
  !   ts              IN      surface skin temperature [K]
  !   nv              IN      number of layers (normally = pver)

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(inout) :: h2ommr_col(nv)
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
    real(r8) :: H, q_mixed
    real(r8) :: Tkp1_new, Tk_new

    block
      use exocol_mod, only: mwdry_col
      kappa = (SHR_CONST_RGAS / mwdry_col) / cp
    end block

    do ipass = 1, max_pass
      adjusted = .false.

      ! Surface coupling handled by bulk SH flux in rce_loop; no
      ! surface-bottom-pair adjustment here (it would be non-conservative).

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
          q_mixed   = (pdel(k)*h2ommr_col(k) + pdel(k+1)*h2ommr_col(k+1)) &
                      / (pdel(k) + pdel(k+1))
          h2ommr_col(k)   = q_mixed
          h2ommr_col(k+1) = q_mixed
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interior and TOA interface temperatures; tint(nv+1) = ts untouched.
    call compute_tint_interp(tmid, pint, nv, tint)

  end subroutine convadj_dry

  ! -----------------------------------------------------------------------
  ! Scheme 1 — Moist pseudo-adiabatic adjustment
  ! -----------------------------------------------------------------------

  subroutine convadj_moist(tmid, tint, h2ommr_col, zint_if, pint, pdel, cp, g, ts, nv)
  ! Moist convective adjustment with rh-weighted local lapse rate.
  !
  ! The relevant critical lapse rate depends on the local saturation state:
  !
  !   Γeff(k) = rh(k)·Γm(T̄,p̄) + (1 − rh(k))·Γd
  !
  ! where Γd = g/cp (dry adiabat), Γm is the saturated moist adiabat, and
  ! rh = q/qsat is computed from the actual specific humidity h2ommr.
  !
  ! Limits: rh→1 (saturated) recovers the Manabe-style moist-adiabatic
  ! adjustment; rh→0 (dry) recovers a dry-adiabatic adjustment.  Earth's
  ! observed ~6.5 K/km tropospheric lapse rate emerges from this blend
  ! without prescribing it.  No saturation threshold is needed: the
  ! interpolation is smooth and parameter-free.
  !
  ! Stability criterion (per pair):
  !   (T_{k+1} − T_k) / (zmid_k − zmid_{k+1}) > Γeff  →  unstable
  !
  ! Γm is evaluated at the mean T and p of the pair via malr():
  !   Γm = (g/cp) · (1 + Lv·ws/(Rd·T)) / (1 + Lv²·ws/(cp·Rv·T²))
  ! ws = eps·esat(T̄)/(p̄ − esat(T̄)).
  !
  ! Adjustment: restore Γeff exactly while conserving cp·T·Δp:
  !   T'_{k+1} = (H + Γeff·Δz·Δp_k) / (Δp_k + Δp_{k+1})
  !   T'_k     = T'_{k+1} − Γeff·Δz
  !
  ! In any pair that is adjusted, q is also mixed — but weighted by the same
  ! rh_pair used for Γeff:
  !   q_k   ← (1 − rh)·q_k   + rh·q_homog
  !   q_k+1 ← (1 − rh)·q_k+1 + rh·q_homog
  ! where q_homog is the mass-weighted average of the pair.  Saturated pairs
  ! (rh → 1, cumulus-like) mix q fully; subsaturated pairs adjust T but
  ! leave q nearly untouched.  Conserves mass-weighted total q in the pair.
  ! The surface-bottom pair does not mix q (slab has no q field; surface
  ! moisture flux is handled by the bulk LE term).

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(inout) :: h2ommr_col(nv) ! specific humidity [kg/kg]
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
    real(r8) :: Rd, eps_wv, gamma_d
    real(r8) :: zmid_k, zmid_kp1, dz
    real(r8) :: T_mean, p_mean, gamma_m, gamma_eff
    real(r8) :: gamma_actual
    real(r8) :: q_pair, es_pair, qsat_pair, rh_pair, q_mixed
    real(r8) :: H, Tkp1_new, Tk_new

    block
      use exocol_mod, only: mwdry_col
      Rd     = SHR_CONST_RGAS / mwdry_col
      eps_wv = SHR_CONST_MWWV / mwdry_col
    end block
    gamma_d = g / cp

    do ipass = 1, max_pass
      adjusted = .false.

      ! NOTE: no surface-bottom-layer pair adjustment.  Slab-atmosphere
      ! coupling is handled by the bulk sensible heat flux in the rce loop,
      ! which is energy-conservative.  A direct tmid(nv) ← ts pull here
      ! would inject energy into the bottom layer without a matching slab
      ! debit, breaking column energy conservation at the surface.

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

        q_pair    = 0.5_r8 * (h2ommr_col(k) + h2ommr_col(k+1))
        es_pair   = min(esat_cc(T_mean), 0.99_r8 * p_mean)
        qsat_pair = eps_wv * es_pair / (p_mean - es_pair)
        rh_pair   = min(q_pair / max(qsat_pair, 1.0e-20_r8), 1.0_r8)
        gamma_eff = rh_pair * gamma_m + (1.0_r8 - rh_pair) * gamma_d

        if (gamma_actual > gamma_eff) then
          H        = pdel(k)*tmid(k) + pdel(k+1)*tmid(k+1)
          Tkp1_new = (H + gamma_eff*dz*pdel(k)) / (pdel(k) + pdel(k+1))
          Tk_new   = Tkp1_new - gamma_eff * dz
          tmid(k)   = Tk_new
          tmid(k+1) = Tkp1_new
          ! rh-weighted q-mixing: saturated pairs homogenize, dry pairs do not.
          q_mixed   = (pdel(k)*h2ommr_col(k) + pdel(k+1)*h2ommr_col(k+1)) &
                      / (pdel(k) + pdel(k+1))
          h2ommr_col(k)   = (1._r8 - rh_pair) * h2ommr_col(k)   + rh_pair * q_mixed
          h2ommr_col(k+1) = (1._r8 - rh_pair) * h2ommr_col(k+1) + rh_pair * q_mixed
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interior and TOA interface temperatures; tint(nv+1) = ts untouched.
    call compute_tint_interp(tmid, pint, nv, tint)

  end subroutine convadj_moist

  ! -----------------------------------------------------------------------
  ! Scheme 2 — Manabe-Wetherald fixed lapse rate
  ! -----------------------------------------------------------------------

  subroutine convadj_manabe(tmid, tint, h2ommr_col, zint_if, pint, pdel, cp, g, ts, nv)
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
  ! In any pair that is adjusted, q is also mixed (mass-weighted by pdel).
  ! Surface-bottom pair does not mix q (slab has no q field; surface
  ! moisture flux handled separately by the bulk-aerodynamic LE term).
  ! zint_if reflects heights from the previous timestep (the caller updates
  ! zint after convective adjustment).

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(inout) :: h2ommr_col(nv)
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
    real(r8) :: H, Tkp1_new, Tk_new, q_mixed

    do ipass = 1, max_pass
      adjusted = .false.

      ! Surface coupling handled by bulk SH flux in rce_loop.

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
          q_mixed   = (pdel(k)*h2ommr_col(k) + pdel(k+1)*h2ommr_col(k+1)) &
                      / (pdel(k) + pdel(k+1))
          h2ommr_col(k)   = q_mixed
          h2ommr_col(k+1) = q_mixed
          adjusted  = .true.
        end if
      end do

      if (.not. adjusted) exit
    end do

    ! Recompute interior and TOA interface temperatures; tint(nv+1) = ts untouched.
    call compute_tint_interp(tmid, pint, nv, tint)

  end subroutine convadj_manabe

  ! -----------------------------------------------------------------------
  ! Shared interface-temperature interpolation
  ! -----------------------------------------------------------------------

  subroutine compute_tint_interp(tmid, pint, nv, tint)
  ! Fill tint(1:nv) by log-pressure interpolation/extrapolation from tmid(1:nv).
  ! tint(nv+1) (surface interface, pinned to ts) is left unchanged — the caller
  ! is responsible for setting it before this routine and after as needed.
  !
  ! Interior interfaces (k+1, k=1..nv-1): log-p linear interpolation between
  ! the midpoints of adjacent layers.
  ! TOA interface (k=1): log-p linear extrapolation from tmid(1) and tmid(2).
    real(r8), intent(in)    :: tmid(nv)
    real(r8), intent(in)    :: pint(nv+1)
    integer,  intent(in)    :: nv
    real(r8), intent(inout) :: tint(nv+1)   ! tint(1:nv) overwritten; tint(nv+1) unchanged

    integer  :: k
    real(r8) :: pmid_k, pmid_kp1

    do k = 1, nv-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      tint(k+1) = tmid(k) + (tmid(k+1) - tmid(k)) * &
                  log(pint(k+1) / pmid_k) / log(pmid_kp1 / pmid_k)
    end do
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine compute_tint_interp

  ! -----------------------------------------------------------------------
  ! Private helpers for moist adiabatic lapse rate
  ! -----------------------------------------------------------------------

  pure function esat_cc(T) result(es)
  ! Clausius-Clapeyron saturation vapour pressure [Pa].
  ! Phase-aware: uses Lvap (over liquid) for T >= T0_sat (273.16 K) and Lsub
  ! (over ice) below.  Both branches give es = es0 = 611.2 Pa at T0_sat,
  ! so the function is continuous at the freezing point.
  !
  ! Defensive: the exponent (L/Rv)(1/T0 - 1/T) overflows IEEE double when T is
  ! negative or close to zero.  We clamp T to a safe physical range [50, 5000] K
  ! before evaluating.  Outside this range the input was unphysical anyway; the
  ! clamp prevents a floating overflow from killing the run.
    real(r8), intent(in) :: T   ! temperature [K]
    real(r8) :: es, L_use, T_use
    real(r8), parameter :: T_min_safe = 50._r8
    real(r8), parameter :: T_max_safe = 5000._r8
    T_use = min(max(T, T_min_safe), T_max_safe)
    if (T_use >= T0_sat) then
      L_use = SHR_CONST_LATVAP
    else
      L_use = SHR_CONST_LATSUB
    end if
    es = es0 * exp((L_use / SHR_CONST_RWV) * (1._r8/T0_sat - 1._r8/T_use))
  end function esat_cc

  pure function Lvap_T(T) result(L)
  ! Phase-appropriate latent heat for evaporation/sublimation [J/kg].
  !   T >= 273.16 K  →  L_v (liquid → vapor)
  !   T <  273.16 K  →  L_s (solid  → vapor, includes fusion)
    real(r8), intent(in) :: T
    real(r8) :: L
    if (T >= T0_sat) then
      L = SHR_CONST_LATVAP
    else
      L = SHR_CONST_LATSUB
    end if
  end function Lvap_T

  pure function malr(T, p, Rd, g_planet, cp_dry) result(gamma_m)
  ! Moist (or ice-) adiabatic lapse rate [K/m] for given T [K], p [Pa], and
  ! gas constants.  The latent-heat coefficient is phase-aware via Lvap_T(T)
  ! so saturated columns below 273.16 K release the sublimation latent heat.
  ! Floored at 1 K/km (= 0.001 K/m) to avoid pathological values where esat → 0.
    real(r8), intent(in) :: T, p, Rd, g_planet, cp_dry
    real(r8) :: gamma_m
    real(r8) :: eps, es, ws, L
    real(r8), parameter :: ws_tiny = 1.0e-10_r8
    real(r8), parameter :: gamma_floor = 1.0e-3_r8  ! 1 K/km

    eps = Rd / SHR_CONST_RWV
    es  = esat_cc(T)
    L   = Lvap_T(T)
    if (p > es) then
      ws = eps * es / (p - es)
    else
      ws = ws_tiny
    end if

    gamma_m = (g_planet / cp_dry) * &
              (1._r8 + L * ws / (Rd * T)) / &
              (1._r8 + L**2 * ws / (cp_dry * SHR_CONST_RWV * T**2))
    gamma_m = max(gamma_m, gamma_floor)
  end function malr

end module exocol_convadj
