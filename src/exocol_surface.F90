module exocol_surface
! Bulk-aerodynamic surface fluxes for ExoColumn.
!
! Computes the turbulent latent (LE) and sensible (SH) heat fluxes between
! the slab ocean and the bottom atmospheric layer using bulk formulas:
!
!   LE = ρ_air · L(Ts) · C_D · U · (qsat(Ts) − q_bot)        [W/m²]
!   SH = ρ_air · cp    · C_D · U · (Ts − T_bot)              [W/m²]
!
! ρ_air is computed from the bottom-layer ideal-gas law p/(Rd·T).  L is the
! phase-aware latent heat from exocol_convadj::Lvap_T (L_v above 0 °C,
! L_sub below).  qsat uses the phase-aware esat_cc.
!
! Sign convention: positive LE and SH transfer energy and water vapour from
! the surface into the atmosphere.  Negative values (dew/frost; warmer air
! over cooler surface) are allowed.
!
! The exchange coefficient C_D and wind speed U are read from
! exocol_config.nml.  This is a 1-D model with no resolved momentum field;
! both must be prescribed.

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_RGAS, SHR_CONST_MWWV
  use exocol_convadj, only: esat_cc, Lvap_T

  implicit none
  private

  public :: compute_surface_fluxes

contains

  subroutine compute_surface_fluxes(ts, t_bot, q_bot, p_bot, mwdry, cpdry, &
                                    wind_speed, C_D, LE, SH)
  ! Bulk aerodynamic latent and sensible heat fluxes at the slab-atmosphere
  ! interface.  All arguments scalar; result delivered through LE/SH.
    real(r8), intent(in)  :: ts          ! surface (slab) temperature [K]
    real(r8), intent(in)  :: t_bot       ! bottom-layer midpoint T    [K]
    real(r8), intent(in)  :: q_bot       ! bottom-layer specific hum. [kg/kg]
    real(r8), intent(in)  :: p_bot       ! bottom-layer midpoint p    [Pa]
    real(r8), intent(in)  :: mwdry       ! dry-air molar mass         [g/mol]
    real(r8), intent(in)  :: cpdry       ! dry-air specific heat      [J/kg/K]
    real(r8), intent(in)  :: wind_speed  ! surface wind speed         [m/s]
    real(r8), intent(in)  :: C_D         ! bulk exchange coefficient  [-]
    real(r8), intent(out) :: LE          ! latent heat flux           [W/m²]
    real(r8), intent(out) :: SH          ! sensible heat flux         [W/m²]

    real(r8) :: Rd, eps_wv, rho_air
    real(r8) :: es_surf, qsat_surf, L_surf

    Rd      = SHR_CONST_RGAS / mwdry
    eps_wv  = SHR_CONST_MWWV / mwdry
    rho_air = p_bot / (Rd * t_bot)

    ! Saturation specific humidity at the slab surface temperature.
    es_surf   = min(esat_cc(ts), 0.99_r8 * p_bot)
    qsat_surf = eps_wv * es_surf / (p_bot - es_surf)

    ! Phase-aware latent heat for evaporation (or sublimation if Ts < 0 °C).
    L_surf = Lvap_T(ts)

    LE = rho_air * L_surf * C_D * wind_speed * (qsat_surf - q_bot)
    SH = rho_air * cpdry  * C_D * wind_speed * (ts - t_bot)

  end subroutine compute_surface_fluxes

end module exocol_surface
