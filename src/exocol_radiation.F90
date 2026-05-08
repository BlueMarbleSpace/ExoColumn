module exocol_radiation
! Radiation wrapper for ExoColumn.
! Mirrors the role of exo_radiation_tend in ExoRT/CESM: reads column state
! from exocol_mod, calls aerad_driver, and returns heating rates (K/day) and
! broadband fluxes (W/m²) at interface levels.
!
! Caller responsibilities before invoking exocol_rad_tend:
!   - exocol_mod state arrays must be populated (exocol_io::read_initial_conditions)
!   - physconst::mwdry / cpair must be set (exocol_mod::exocol_setgas)
!   - pdeldry and pintdry must be current (exocol_mod::exocol_update_derived)
!   - ExoRT must be initialised (initialize_kcoeff, initialize_solar, etc.)

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use shr_const_mod,     only: SHR_CONST_CSEC
  use ppgrid,            only: pver, pverp
  use radgrid,           only: ntot_wavlnrng
  use exo_radiation_mod, only: aerad_driver
  use input,             only: ext_nazm_tshadow
  use exocol_mod

  implicit none
  private

  public :: exocol_rad_tend

  ! Topography / shadow parameters: flat, no shadowing (identical to
  ! the no-topography defaults used in ExoRT main.F90 and exo_radiation_tend).
  ! ext_nazm_tshadow is taken from ExoRT's input module.
  integer,  parameter :: ext_tslas_tog     = 0
  integer,  parameter :: ext_tshadow_tog   = 1
  real(r8), parameter :: ext_rtgt          = 1._r8
  real(r8), parameter :: ext_solar_azm_ang = 0._r8
  real(r8), parameter :: ext_tazm_ang      = 0._r8
  real(r8), parameter :: ext_tslope_ang    = 0._r8

contains

  subroutine exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)
  ! Call aerad_driver for the current column state and return
  ! heating rates and broadband fluxes.
  !
  !   LWHR, SWHR  [K/day]   LW and SW heating rates at layer midpoints
  !   LWUP, LWDN  [W/m²]    LW up/down fluxes at interfaces (TOA=1, srf=pverp)
  !   SWUP, SWDN  [W/m²]    SW up/down fluxes at interfaces

    real(r8), intent(out), dimension(pver)  :: LWHR
    real(r8), intent(out), dimension(pver)  :: SWHR
    real(r8), intent(out), dimension(pverp) :: LWUP
    real(r8), intent(out), dimension(pverp) :: LWDN
    real(r8), intent(out), dimension(pverp) :: SWUP
    real(r8), intent(out), dimension(pverp) :: SWDN

    ! Raw aerad_driver outputs (heating rates in K/s before conversion)
    real(r8), dimension(pver)  :: sw_dTdt, lw_dTdt

    ! Spectral flux arrays required by aerad_driver signature but not used
    ! by the RCE loop; allocated on the stack and discarded.
    real(r8), dimension(pverp, ntot_wavlnrng) :: lw_dnflux_spec
    real(r8), dimension(pverp, ntot_wavlnrng) :: lw_upflux_spec
    real(r8), dimension(pverp, ntot_wavlnrng) :: sw_upflux_spec
    real(r8), dimension(pverp, ntot_wavlnrng) :: sw_dnflux_spec

    ! Surface SW component diagnostics (not used in RCE loop)
    real(r8) :: vis_dir, vis_dif, nir_dir, nir_dif, sol_toa

    ! Shadow arrays (single-element, all zero: no topography)
    real(r8), dimension(ext_nazm_tshadow) :: ext_cosz_horizon
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCx_obstruct
    real(r8), dimension(ext_nazm_tshadow) :: ext_TCz_obstruct

    ext_cosz_horizon(:) = 0._r8
    ext_TCx_obstruct(:) = 0._r8
    ext_TCz_obstruct(:) = 0._r8

    ! Argument order matches aerad_driver signature in exo_radiation_mod.F90:
    !   gases, clouds, pressure/temperature grid, albedos,
    !   topography parameters, heights, then output arrays.
    call aerad_driver( &
        h2ommr,   co2mmr,                                         &  ! H2O, CO2
        ch4mmr,   c2h6mmr,                                        &  ! CH4, C2H6
        h2mmr,    n2mmr,    o3mmr,   o2mmr,                       &  ! H2, N2, O3, O2
        cicewp,   cliqwp,   cfrc,                                  &  ! cloud IWP, LWP, fraction
        rei,      rel,                                             &  ! cloud particle radii
        ts,       ps,       pmid,                                  &  ! sfc T, sfc P, midpoint P
        pdel,     pdeldry,  tmid,    pint,    pintdry,            &  ! dp, dp_dry, T_mid, p_int
        coszrs,   msdist,                                          &  ! zenith angle, star distance
        asdir,    aldir,                                           &  ! vis/nir direct albedo
        asdif,    aldif,                                           &  ! vis/nir diffuse albedo
        ext_rtgt, ext_solar_azm_ang, ext_tazm_ang, ext_tslope_ang, &
        ext_tslas_tog, ext_tshadow_tog, ext_nazm_tshadow,          &
        ext_cosz_horizon, ext_TCx_obstruct, ext_TCz_obstruct,      &
        zint,                                                       &  ! interface heights [m]
        sw_dTdt,  lw_dTdt,                                         &  ! OUT: heating rates [K/s]
        LWDN,     LWUP,     SWUP,    SWDN,                         &  ! OUT: broadband fluxes [W/m²]
        lw_dnflux_spec, lw_upflux_spec, sw_upflux_spec, sw_dnflux_spec, &
        vis_dir,  vis_dif,  nir_dir, nir_dif, sol_toa              )

    ! Convert raw K/s output to K/day for the RCE time-marching loop
    SWHR(:) = sw_dTdt(:) * SHR_CONST_CSEC
    LWHR(:) = lw_dTdt(:) * SHR_CONST_CSEC

  end subroutine exocol_rad_tend

end module exocol_radiation
