module exocol_convadj
! Convective adjustment for ExoColumn.
!
! Phase 1 (implemented here): dry adiabatic adjustment.
!   For each pair of adjacent layers that is super-adiabatic (potential
!   temperature decreasing upward), redistributes enthalpy to restore
!   neutral stability while conserving the total column enthalpy of the
!   adjusted pair.  Sweeps from surface upward; repeats until stable or
!   max_pass is exhausted.
!
! Phase 2 (future): moist adiabatic adjustment.  Will replace the
!   dry-adiabat stability criterion with the local moist adiabatic lapse
!   rate computed from T(z) and q(z).

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_RGAS
  use ppgrid,        only: pver, pverp

  implicit none
  private

  public :: convadj_dry

contains

  subroutine convadj_dry(tmid, tint, pint, pdel, cp, g, nv)
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
  !   g            IN      gravitational acceleration [m/s²] (reserved for
  !                        moist phase; not used in this potential-T version)
  !   nv           IN      number of layers (normally = pver)
  !
  ! Stability criterion: potential temperature θ(k) = T(k)·(p_ref/p(k))^κ
  ! must be non-decreasing upward (toward smaller k / lower pressure).
  ! Instability: θ(k) < θ(k+1) where k is above k+1.
  !
  ! Adjustment: conserve enthalpy of the pair and restore θ(k) = θ(k+1).
  ! This gives (letting r = (pmid_k/pmid_{k+1})^κ, H = pdel_k*T_k + pdel_{k+1}*T_{k+1}):
  !   T'_{k+1} = H / (pdel_k * r + pdel_{k+1})
  !   T'_k     = r * T'_{k+1}

    real(r8), intent(inout) :: tmid(nv)
    real(r8), intent(inout) :: tint(nv+1)
    real(r8), intent(in)    :: pint(nv+1)
    real(r8), intent(in)    :: pdel(nv)
    real(r8), intent(in)    :: cp
    real(r8), intent(in)    :: g
    integer,  intent(in)    :: nv

    integer, parameter :: max_pass = 30

    integer  :: k, ipass
    logical  :: adjusted
    real(r8) :: pmid_k, pmid_kp1
    real(r8) :: kappa          ! R/cp for this gas mixture
    real(r8) :: ratio_kappa    ! (pmid_k / pmid_{k+1})^kappa
    real(r8) :: theta_k_scaled ! tmid(k)  scaled to pmid_{k+1} reference
    real(r8) :: H              ! total enthalpy of the pair [K·Pa]
    real(r8) :: Tkp1_new, Tk_new

    ! Specific gas constant R = SHR_CONST_RGAS / mwdry_col  [J/kg/K]
    ! SHR_CONST_RGAS is in J/kmol/K; mwdry_col in g/mol = kg/kmol
    ! κ = R/cp = (SHR_CONST_RGAS / mwdry_col) / cp
    ! Access mwdry_col through host-associated USE in the calling module;
    ! here we call the accessor via use-association from exocol_mod.
    block
      use exocol_mod, only: mwdry_col
      kappa = (SHR_CONST_RGAS / mwdry_col) / cp
    end block

    do ipass = 1, max_pass
      adjusted = .false.

      ! Sweep from surface (k = nv-1) upward (toward k = 1).
      ! Pair (k, k+1): k is the upper layer, k+1 is the lower layer.
      do k = nv-1, 1, -1
        pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
        pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))

        ! Scale tmid(k) to the pmid_{k+1} reference frame so that
        ! comparing theta_k_scaled and tmid(k+1) is equivalent to
        ! comparing potential temperatures.
        ratio_kappa    = (pmid_k / pmid_kp1)**kappa
        theta_k_scaled = tmid(k) / ratio_kappa   ! = θ(k)·(p_{k+1}/p_ref)^κ

        if (theta_k_scaled < tmid(k+1)) then     ! unstable: θ(k) < θ(k+1)
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
    ! Interpolate linearly in log-pressure between adjacent layer midpoints.
    do k = 1, nv-1
      pmid_k   = 0.5_r8 * (pint(k)   + pint(k+1))
      pmid_kp1 = 0.5_r8 * (pint(k+1) + pint(k+2))
      tint(k+1) = tmid(k) + (tmid(k+1) - tmid(k)) * &
                  log(pint(k+1)/pmid_k) / log(pmid_kp1/pmid_k)
    end do
    ! Top interface: linear extrapolation from the two uppermost midpoints.
    pmid_k   = 0.5_r8 * (pint(1) + pint(2))
    pmid_kp1 = 0.5_r8 * (pint(2) + pint(3))
    tint(1) = tmid(1) - (tmid(2) - tmid(1)) * &
              log(pmid_k / pint(1)) / log(pmid_kp1 / pmid_k)

  end subroutine convadj_dry

end module exocol_convadj
