module exocol_surface
! Surface turbulent fluxes for ExoColumn.
!
! Two schemes, selected by exocol_config::surface_flux:
!
!   'mos'  (default) — simplified Monin-Obukhov similarity, following Frierson,
!          Held & Zurita-Gotor (2006, J. Atmos. Sci. 63, 2548), the lineage of
!          ExoColumn's SBM convection.  Fluxes use POTENTIAL-temperature (dry
!          static energy) differences and a drag coefficient that depends on the
!          lowest-level height z_a through the surface-layer log law:
!
!            C = κ² / [ln(z_a/z0)]²                       (R_ia < 0, unstable→neutral)
!            C = κ² / [ln(z_a/z0)]² · (1 − R_ia/Ri_c)²    (0 < R_ia < Ri_c, stable)
!            C = 0                                        (R_ia > Ri_c)
!
!            SH = ρ_a·cp·C·U·(Θ_s − Θ_a)        Θ_s = Ts, Θ_a = T_a·(p0/p_a)^κ
!            LE = ρ_a·L(Ts)·C·U·(q*_s − q_a)
!
!          The z_a dependence of C is what makes the flux RESOLUTION-INDEPENDENT:
!          as the grid is refined and z_a shrinks, C grows and (Θ_s − Θ_a)
!          shrinks so their product (the constant-flux-layer value) is preserved.
!          This is only valid when z_a sits inside the surface layer (z_a ≲ 0.1·h),
!          i.e. with the hybrid fine-surface grid (n_sfc_layers > 0, z_a ~ 10 m).
!          The unstable branch is treated as neutral (Φ = 1), so no Obukhov-length
!          iteration is needed.  R_ia uses the VIRTUAL dry static energy
!          s_v = cp·T_v + g·z (T_v = T(1+0.608 q)).
!
!   'bulk' (legacy) — fixed-C_D bulk aerodynamic formulas in ACTUAL temperature,
!          SH = ρ·cp·C_D·U·(Ts − T_a), LE = ρ·L·C_D·U·(q*_s − q_a).  Retained for
!          comparison; resolution-dependent near-surface state (see
!          project_resolution_root_cause).
!
! Sign convention: positive LE, SH transfer energy / water from surface to
! atmosphere.  compute_surface_fluxes also returns the exchange coefficient C so
! the boundary-layer scheme can build a consistent eddy diffusivity.

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_RGAS, SHR_CONST_MWWV
  use exoplanet_mod, only: exo_g
  use exocol_convadj, only: esat_cc, Lvap_T

  implicit none
  private

  public :: compute_surface_fluxes
  public :: mos_drag_coef

  real(r8), parameter :: vk    = 0.4_r8        ! von Karman constant
  real(r8), parameter :: ri_c  = 1.0_r8        ! critical bulk Richardson (Frierson)

contains

  subroutine compute_surface_fluxes(mode, ts, t_bot, q_bot, p_bot, za, p0, z0, &
                                    mwdry, cpdry, wind_speed, C_D, LE, SH, C_out)
  ! Surface latent and sensible heat fluxes at the slab-atmosphere interface.
    character(len=*), intent(in)  :: mode        ! 'mos' or 'bulk'
    real(r8), intent(in)  :: ts          ! surface (slab) temperature [K]
    real(r8), intent(in)  :: t_bot       ! lowest-level midpoint T     [K]
    real(r8), intent(in)  :: q_bot       ! lowest-level specific hum.  [kg/kg]
    real(r8), intent(in)  :: p_bot       ! lowest-level midpoint p     [Pa]
    real(r8), intent(in)  :: za          ! lowest-level height         [m]
    real(r8), intent(in)  :: p0          ! surface pressure            [Pa]
    real(r8), intent(in)  :: z0          ! surface roughness length    [m]
    real(r8), intent(in)  :: mwdry       ! dry-air molar mass          [g/mol]
    real(r8), intent(in)  :: cpdry       ! dry-air specific heat       [J/kg/K]
    real(r8), intent(in)  :: wind_speed  ! surface wind speed          [m/s]
    real(r8), intent(in)  :: C_D         ! bulk exchange coefficient (legacy) [-]
    real(r8), intent(out) :: LE          ! latent heat flux            [W/m²]
    real(r8), intent(out) :: SH          ! sensible heat flux          [W/m²]
    real(r8), intent(out) :: C_out       ! exchange coefficient used   [-]

    real(r8) :: Rd, eps_wv, kap, rho_a
    real(r8) :: es_surf, qsat_surf, L_surf, theta_a, theta_s
    real(r8) :: sv_s, sv_a, Ria, C

    Rd      = SHR_CONST_RGAS / mwdry
    eps_wv  = SHR_CONST_MWWV / mwdry
    kap     = Rd / cpdry
    rho_a   = p_bot / (Rd * t_bot)
    L_surf  = Lvap_T(ts)

    if (trim(adjustl(mode)) == 'bulk') then
      ! ---- legacy fixed-C_D bulk aerodynamic (actual temperature) ----
      ! qsat referenced to the bottom-midpoint pressure p_bot (matches the
      ! calibrated Ts=288.01 reference).
      es_surf   = min(esat_cc(ts), 0.99_r8 * p_bot)
      qsat_surf = eps_wv * es_surf / (p_bot - es_surf)
      C  = C_D
      LE = rho_a * L_surf * C * wind_speed * (qsat_surf - q_bot)
      SH = rho_a * cpdry  * C * wind_speed * (ts - t_bot)
    else
      ! ---- simplified Monin-Obukhov similarity (Frierson 2006) ----
      ! qsat at the true surface pressure p0.
      es_surf   = min(esat_cc(ts), 0.99_r8 * p0)
      qsat_surf = eps_wv * es_surf / (p0 - es_surf)
      ! Potential temperatures (surface reference pressure p0): Θ_s = Ts.
      theta_s = ts
      theta_a = t_bot * (p0 / p_bot)**kap

      ! Bulk Richardson at the lowest level from virtual dry static energy.
      sv_s = cpdry * ts    * (1._r8 + 0.608_r8 * qsat_surf)              ! z = 0
      sv_a = cpdry * t_bot * (1._r8 + 0.608_r8 * q_bot) + exo_g * za
      Ria  = exo_g * za * (sv_a - sv_s) / (sv_s * max(wind_speed, 1.0_r8)**2)

      C  = mos_drag_coef(za, z0, Ria)
      LE = rho_a * L_surf * C * wind_speed * (qsat_surf - q_bot)
      SH = rho_a * cpdry  * C * wind_speed * (theta_s - theta_a)
    end if

    C_out = C

  end subroutine compute_surface_fluxes

  ! -----------------------------------------------------------------------

  pure function mos_drag_coef(za, z0, Ria) result(C)
  ! Simplified Monin-Obukhov drag coefficient (Frierson 2006, eqs 12-14).
  ! Unstable (Ria<0) is treated as neutral; stable is reduced toward zero at Ri_c.
    real(r8), intent(in) :: za, z0, Ria
    real(r8) :: C, Cn
    Cn = (vk / log(za / z0))**2
    if (Ria < 0._r8) then
      C = Cn
    else if (Ria < ri_c) then
      C = Cn * (1._r8 - Ria / ri_c)**2
    else
      C = 0._r8
    end if
  end function mos_drag_coef

end module exocol_surface
