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
  public :: convadj_zm
  public :: convadj_sbm
  public :: esat_cc
  public :: Lvap_T
  public :: malr
  public :: compute_tint_interp
  public :: compute_cape

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
  ! Scheme 3 — Zhang-McFarlane-style soft moist adjustment
  ! -----------------------------------------------------------------------

  subroutine convadj_zm(tmid, tint, h2ommr_col, zint_if, pint, pdel, cp, g, ts, &
                        f_relax, cape_trig, nv)
  ! Single-sweep moist adjustment with relaxation fraction f_relax.
  !
  ! For each unstable layer pair (same criterion as convadj_moist) the full
  ! moist-adiabatic hard-adjustment increments ΔT and Δq are computed but
  ! only a fraction f_relax ∈ (0,1] is applied.  Setting
  !   f_relax = 1 - exp(-dt / τ_conv)
  ! with the ZM default τ_conv = 7200 s reproduces the CAPE-relaxation
  ! closure of Zhang & McFarlane (1995) without an explicit mass-flux
  ! framework: each timestep removes 1 - e^{-dt/τ} of the layer-pair
  ! instability.  f_relax = 1 recovers a single hard-adjustment pass
  ! (used for the post-condensation cleanup in run_rce_loop).
  !
  ! CAPE trigger: when cape_trig > 0, a surface-parcel CAPE is computed
  ! first; convection is suppressed if CAPE < cape_trig.  Pass cape_trig = 0
  ! to skip the check (always active on any unstable pair).

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(inout) :: h2ommr_col(nv)
    real(r8), intent(in)    :: zint_if(nv+1)
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g
    real(r8), intent(in)    :: ts
    real(r8), intent(in)    :: f_relax    ! relaxation fraction [0,1]
    real(r8), intent(in)    :: cape_trig  ! CAPE activation threshold [J/kg]
    integer,  intent(in)    :: nv

    real(r8), parameter :: dz_min = 1._r8

    integer  :: k
    real(r8) :: Rd, eps_wv, gamma_d
    real(r8) :: pmid_local(nv)
    real(r8) :: zmid_k, zmid_kp1, dz
    real(r8) :: T_mean, p_mean, gamma_m, gamma_eff, gamma_actual
    real(r8) :: q_pair, es_pair, qsat_pair, rh_pair
    real(r8) :: H, Tkp1_hard, Tk_hard, dTk, dTkp1
    real(r8) :: q_mixed, dq_k, dq_kp1

    block
      use exocol_mod, only: mwdry_col, pmid
      Rd         = SHR_CONST_RGAS / mwdry_col
      eps_wv     = SHR_CONST_MWWV / mwdry_col
      pmid_local = pmid(1:nv)
    end block
    gamma_d = g / cp

    if (cape_trig > 0.0_r8) then
      if (compute_cape(tmid, h2ommr_col, pmid_local, pint, zint_if, &
                       Rd, g, cp, eps_wv, nv) < cape_trig) return
    end if

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
        H         = pdel(k)*tmid(k) + pdel(k+1)*tmid(k+1)
        Tkp1_hard = (H + gamma_eff*dz*pdel(k)) / (pdel(k) + pdel(k+1))
        Tk_hard   = Tkp1_hard - gamma_eff * dz
        dTk       = f_relax * (Tk_hard   - tmid(k))
        dTkp1     = f_relax * (Tkp1_hard - tmid(k+1))
        tmid(k)   = tmid(k)   + dTk
        tmid(k+1) = tmid(k+1) + dTkp1
        q_mixed   = (pdel(k)*h2ommr_col(k) + pdel(k+1)*h2ommr_col(k+1)) &
                    / (pdel(k) + pdel(k+1))
        dq_k      = f_relax * rh_pair * (q_mixed - h2ommr_col(k))
        dq_kp1    = f_relax * rh_pair * (q_mixed - h2ommr_col(k+1))
        h2ommr_col(k)   = h2ommr_col(k)   + dq_k
        h2ommr_col(k+1) = h2ommr_col(k+1) + dq_kp1
      end if
    end do

    call compute_tint_interp(tmid, pint, nv, tint)

  end subroutine convadj_zm

  ! -----------------------------------------------------------------------
  ! Scheme 4 — Simplified Betts-Miller (Frierson 2007)
  ! -----------------------------------------------------------------------

  subroutine convadj_sbm(tmid, tint, h2ommr_col, zint_if, pint, pdel, cp, g, ts, &
                         dt_sec, tau_sbm, rh_ref, L_release, do_moisture, &
                         precip_mass_flux, cond_tend, nv)
  ! Simplified Betts-Miller convective adjustment (Frierson 2007, JAS 64:1959).
  !
  ! Relaxes the convecting column toward two reference profiles over a
  ! convective timescale tau_sbm:
  !   T  →  T_ref   (moist adiabat lifted from the lowest model level)
  !   q  →  q_ref = rh_ref · qsat(T_ref)        (only when do_moisture)
  ! with relaxation fraction alpha = min(dt/tau_sbm, 1).  alpha → 1 (hard
  ! adjustment to the moist adiabat) when dt ≥ tau_sbm, recovering konrad's
  ! HardAdjustment + MoistLapseRate configuration near equilibrium; alpha < 1
  ! during fast transients smooths convective temperature jumps (TOA noise).
  !
  ! Unlike the rh-weighted convadj_moist/convadj_zm schemes, the TEMPERATURE
  ! target is the pure moist adiabat regardless of the environmental relative
  ! humidity (convective plumes are saturated even when the mean column is
  ! subsaturated).  Relative humidity enters only through the MOISTURE target
  ! q_ref — the standard separation of T and q closures (Betts-Miller 1986).
  ! This is what keeps the troposphere moist-adiabatic and avoids the dry-aloft,
  ! near-dry-adiabatic lapse rate the rh-weighted schemes produce.
  !
  ! Energy conservation (Frierson's correction): the column-integrated enthalpy
  ! added by the temperature relaxation is forced to equal the latent heat
  ! released by the condensed (precipitated) water, by shifting the reference
  ! temperature profile by a single constant dT_shift:
  !   dT_shift = ( Σ cp(T_ref−T)Δp/g − L·Σ(q−q_ref)Δp/g ) / (cp·Σ Δp/g)
  !   T_ref ← T_ref − dT_shift
  ! After the shift  Σ cp·ΔT·Δp/g = L·(precip mass)  exactly, for any alpha.
  ! Water is conserved by construction: precip mass = −Σ Δq·Δp/g = alpha·Wvap.
  !
  ! L_release is the latent heat used for the energy balance; the caller passes
  ! Lvap_T(ts) so the release matches the L(ts) debited from the slab by surface
  ! evaporation (the column latent-heat ledger; see exocol_rce_loop::condense).
  !
  ! Cloud base is the lowest model level (parcel source = tmid(nv)); cloud top
  ! is the highest contiguous buoyant level (T_ref ≥ T_env).  Above cloud top
  ! the column is left to radiative equilibrium.  When the column is net
  ! subsaturated relative to the reference (Σ(q−q_ref) ≤ 0 — no condensation to
  ! sustain deep convection) the scheme makes no change; use conv_scheme='dry'
  ! if dry convective adjustment is needed for a very dry atmosphere.
  !
  ! Simplification: the parcel is lifted moist-adiabatically from the lowest
  ! level (no explicit sub-LCL dry-adiabatic segment).  Adequate when the
  ! lowest layer is near-saturated, as in a moist RCE.

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(inout) :: h2ommr_col(nv)
    real(r8), intent(in)    :: zint_if(nv+1)
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g
    real(r8), intent(in)    :: ts            ! reserved (parcel base = tmid(nv))
    real(r8), intent(in)    :: dt_sec
    real(r8), intent(in)    :: tau_sbm
    real(r8), intent(in)    :: rh_ref
    real(r8), intent(in)    :: L_release
    logical,  intent(in)    :: do_moisture
    real(r8), intent(out)   :: precip_mass_flux    ! [kg/m²/s]
    real(r8), intent(out)   :: cond_tend(nv)       ! convective T tendency [K/s]
    integer,  intent(in)    :: nv

    real(r8), parameter :: dz_min = 1._r8

    integer  :: k, k_top
    real(r8) :: Rd, eps_wv, alpha
    real(r8) :: Tref(nv), qref(nv)
    real(r8) :: zmid_k, zmid_kp1, dz, pmid_k, es, dT
    real(r8) :: Mcol, Qheat, Wvap, dT_shift

    precip_mass_flux = 0._r8
    cond_tend(:)     = 0._r8

    block
      use exocol_mod, only: mwdry_col
      Rd     = SHR_CONST_RGAS / mwdry_col
      eps_wv = SHR_CONST_MWWV / mwdry_col
    end block

    ! 1. Reference moist adiabat, lifted from the lowest model level.
    Tref(nv) = tmid(nv)
    do k = nv-1, 1, -1
      zmid_k   = 0.5_r8 * (zint_if(k)   + zint_if(k+1))
      zmid_kp1 = 0.5_r8 * (zint_if(k+1) + zint_if(k+2))
      dz = max(zmid_k - zmid_kp1, dz_min)
      Tref(k) = Tref(k+1) - malr(Tref(k+1), pint(k+1), Rd, g, cp) * dz
    end do

    ! 2. Cloud top = highest contiguous buoyant level above the surface.
    k_top = nv
    do k = nv-1, 1, -1
      if (Tref(k) >= tmid(k)) then
        k_top = k
      else
        exit
      end if
    end do
    if (k_top >= nv) return        ! no convecting layer

    ! 3. Reference humidity over the convecting column.
    do k = k_top, nv
      pmid_k  = 0.5_r8 * (pint(k) + pint(k+1))
      es      = min(esat_cc(Tref(k)), 0.99_r8 * pmid_k)
      qref(k) = rh_ref * eps_wv * es / (pmid_k - es)
    end do

    alpha = min(dt_sec / tau_sbm, 1.0_r8)

    if (.not. do_moisture) then
      ! Temperature-only relaxation to the moist adiabat (q held fixed).
      do k = k_top, nv
        dT           = alpha * (Tref(k) - tmid(k))
        cond_tend(k) = dT / dt_sec
        tmid(k)      = tmid(k) + dT
      end do
      call compute_tint_interp(tmid, pint, nv, tint)
      return
    end if

    ! 4. Energy-conserving reference-temperature shift.
    Mcol = 0._r8; Qheat = 0._r8; Wvap = 0._r8
    do k = k_top, nv
      Mcol  = Mcol  + pdel(k) / g
      Qheat = Qheat + cp * (Tref(k) - tmid(k)) * pdel(k) / g
      Wvap  = Wvap  + (h2ommr_col(k) - qref(k)) * pdel(k) / g
    end do

    if (Wvap <= 0._r8 .or. Mcol <= 0._r8) return   ! too dry: no deep convection

    dT_shift = (Qheat - L_release * Wvap) / (cp * Mcol)

    ! 5. Apply relaxation: T → T_ref − dT_shift, q → q_ref.
    do k = k_top, nv
      dT            = alpha * ((Tref(k) - dT_shift) - tmid(k))
      cond_tend(k)  = dT / dt_sec
      tmid(k)       = tmid(k) + dT
      h2ommr_col(k) = max(h2ommr_col(k) + alpha * (qref(k) - h2ommr_col(k)), 0._r8)
    end do
    precip_mass_flux = alpha * Wvap / dt_sec

    call compute_tint_interp(tmid, pint, nv, tint)

  end subroutine convadj_sbm

  ! -----------------------------------------------------------------------
  ! CAPE diagnostic — surface parcel ascent
  ! -----------------------------------------------------------------------

  function compute_cape(tmid_col, h2ommr_col, pmid_col, pint_col, zint_if, &
                        Rd, g_planet, cp, eps_wv, nv) result(cape_out)
  ! Compute CAPE [J/kg] by lifting the surface-layer parcel.
  !
  ! Phase 1 (dry adiabat): lift from k=nv to the LCL — the first layer k
  ! going upward where q_parcel ≥ qsat(T_dry_parcel, p).
  ! Phase 2 (moist adiabat): above the LCL, follow malr(T,p) and integrate
  ! buoyancy using virtual temperatures.
  ! CAPE = Σ g·(Tv_par − Tv_env)/Tv_env·Δz  over buoyant layers;
  ! integration stops at the first non-buoyant layer (LNB).
  !
  ! Returns 0 if no LCL is found or the parcel is never buoyant above it.
  ! Virtual temperature: Tv = T·(ε + q·(1−ε))/ε  where ε = Mwv/Mdry.

    real(r8), intent(in) :: tmid_col(nv)
    real(r8), intent(in) :: h2ommr_col(nv)
    real(r8), intent(in) :: pmid_col(nv)
    real(r8), intent(in) :: pint_col(nv+1)
    real(r8), intent(in) :: zint_if(nv+1)
    real(r8), intent(in) :: Rd, g_planet, cp, eps_wv
    integer,  intent(in) :: nv
    real(r8) :: cape_out

    integer  :: k, k_lcl
    real(r8) :: kappa, T_parcel, q_parcel
    real(r8) :: T_parcel_dry, es_dry, qsat_dry
    real(r8) :: zmid_k, zmid_prev
    real(r8) :: T_par, es_par, q_par
    real(r8) :: Tv_par, Tv_env, dz

    cape_out = 0.0_r8
    kappa    = Rd / cp
    T_parcel = tmid_col(nv)
    q_parcel = h2ommr_col(nv)

    ! Phase 1: dry adiabatic lift to LCL
    k_lcl = -1
    do k = nv, 1, -1
      T_parcel_dry = T_parcel * (pmid_col(k) / pmid_col(nv))**kappa
      es_dry       = min(esat_cc(T_parcel_dry), 0.99_r8 * pmid_col(k))
      qsat_dry     = eps_wv * es_dry / (pmid_col(k) - es_dry)
      if (q_parcel >= qsat_dry) then
        k_lcl = k
        T_par = T_parcel_dry
        exit
      end if
    end do
    if (k_lcl < 0) return

    ! Phase 2: moist adiabatic lift above LCL, accumulate CAPE
    zmid_prev = 0.5_r8 * (zint_if(k_lcl) + zint_if(k_lcl+1))

    do k = k_lcl - 1, 1, -1
      zmid_k = 0.5_r8 * (zint_if(k) + zint_if(k+1))
      dz     = zmid_k - zmid_prev

      T_par = T_par - malr(T_par, pint_col(k+1), Rd, g_planet, cp) * dz

      es_par = min(esat_cc(T_par), 0.99_r8 * pmid_col(k))
      q_par  = eps_wv * es_par / (pmid_col(k) - es_par)

      Tv_par = T_par        * (eps_wv + q_par           * (1.0_r8 - eps_wv)) / eps_wv
      Tv_env = tmid_col(k)  * (eps_wv + h2ommr_col(k)   * (1.0_r8 - eps_wv)) / eps_wv

      if (Tv_par > Tv_env) then
        cape_out = cape_out + g_planet * (Tv_par - Tv_env) / Tv_env * dz
      else
        exit
      end if

      zmid_prev = zmid_k
    end do

  end function compute_cape

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
