module exocol_radiation
! Radiation wrapper for ExoColumn.
! Mirrors the role of exo_radiation_tend in ExoRT/CESM: reads column state
! from exocol_mod, calls aerad_driver, and returns heating rates (K/day) and
! broadband fluxes (W/m²) at interface levels.
!
! Solar zenith angle (see &exocol_nml::sw_zenith_quad in exocol_config):
!   - Default (sw_zenith_quad=.false.): one aerad_driver call at the prescribed
!     single angle exocol_mod::coszrs.  Legacy/calibrated path — BIT-IDENTICAL to
!     the historical wrapper (the call is made exactly once, outputs assigned
!     directly, no weighting arithmetic).
!   - sw_zenith_quad=.true.: the shortwave is integrated over the illuminated
!     hemisphere with an sw_nquad-point Gauss-Legendre quadrature (Kopparapu 2013
!     style), giving the diurnal/global-mean (Bond) shortwave.  aerad_driver is
!     called once per Gauss node mu_i (coszrs replaced by mu_i) and every output
!     is combined as  X = sum_i w_i X(mu_i)  with Gauss-Legendre weights on [0,1]
!     (sum w_i = 1).  Because the shortwave is linear in the incident beam and the
!     TOA insolation scales as mu_i (exo_radiation_mod.F90:1760), this reproduces
!     the global-mean insolation S0/4 (= the coszrs=0.5 single-angle value, since
!     sum w_i mu_i = INT_0^1 mu dmu = 1/2) while correcting the angular albedo.
!     The longwave is zenith-independent, so sum_i w_i LW(mu_i) = LW (weights
!     sum to 1) — the LW comes out unchanged automatically.
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
  use exocol_config,     only: sw_zenith_quad, sw_nquad
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

  integer,  parameter :: MAXQ = 8   ! largest supported Gauss-Legendre order

contains

  subroutine exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)
  ! Return heating rates and broadband fluxes for the current column state,
  ! optionally hemispherically averaged over solar zenith angle.
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

    real(r8) :: mu(MAXQ), wq(MAXQ)
    integer  :: nq, iq

    ! Per-call (single zenith angle) outputs, accumulated in the quadrature loop.
    real(r8), dimension(pver)  :: LWHR1, SWHR1
    real(r8), dimension(pverp) :: LWUP1, LWDN1, SWUP1, SWDN1
    real(r8), dimension(ntot_wavlnrng) :: blwup1, bswup1, bswdn1

    if (.not. sw_zenith_quad) then
      ! Legacy single-angle path: one call, outputs assigned directly.  Writing
      ! the module band arrays in place keeps this BIT-IDENTICAL to the original.
      call aerad_once(coszrs, LWHR, SWHR, LWUP, LWDN, SWUP, SWDN, &
                      band_lwup_toa, band_swup_toa, band_swdn_toa)
      return
    end if

    ! Hemispheric Gauss-Legendre quadrature over zenith angle.
    call zenith_quad_nodes(sw_nquad, nq, mu, wq)

    LWHR = 0._r8; SWHR = 0._r8
    LWUP = 0._r8; LWDN = 0._r8; SWUP = 0._r8; SWDN = 0._r8
    band_lwup_toa = 0._r8; band_swup_toa = 0._r8; band_swdn_toa = 0._r8

    do iq = 1, nq
      call aerad_once(mu(iq), LWHR1, SWHR1, LWUP1, LWDN1, SWUP1, SWDN1, &
                      blwup1, bswup1, bswdn1)
      LWHR = LWHR + wq(iq) * LWHR1
      SWHR = SWHR + wq(iq) * SWHR1
      LWUP = LWUP + wq(iq) * LWUP1
      LWDN = LWDN + wq(iq) * LWDN1
      SWUP = SWUP + wq(iq) * SWUP1
      SWDN = SWDN + wq(iq) * SWDN1
      band_lwup_toa = band_lwup_toa + wq(iq) * blwup1
      band_swup_toa = band_swup_toa + wq(iq) * bswup1
      band_swdn_toa = band_swdn_toa + wq(iq) * bswdn1
    end do

  end subroutine exocol_rad_tend

  ! -----------------------------------------------------------------------

  subroutine aerad_once(cz, LWHR, SWHR, LWUP, LWDN, SWUP, SWDN, &
                        blwup, bswup, bswdn)
  ! Single aerad_driver call at solar zenith cosine `cz`.  Returns heating rates
  ! already converted to K/day, broadband interface fluxes [W/m²], and the
  ! band-resolved TOA fluxes [W/m²].  All column state is read from exocol_mod.
    real(r8), intent(in)  :: cz
    real(r8), intent(out), dimension(pver)  :: LWHR, SWHR
    real(r8), intent(out), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN
    real(r8), intent(out), dimension(ntot_wavlnrng) :: blwup, bswup, bswdn

    ! Raw aerad_driver outputs (heating rates in K/s before conversion)
    real(r8), dimension(pver) :: sw_dTdt, lw_dTdt

    ! Spectral flux arrays required by aerad_driver signature; only TOA row used.
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
        cz,       msdist,                                          &  ! zenith angle, star distance
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

    ! Band-resolved TOA fluxes (interface index 1 = TOA, since camtop = 1 in the
    ! n68equiv 1-D path).  band_lwup_toa is the spectrally resolved OLR; the SW
    ! pair gives per-band planetary albedo (inner-HZ continuum/window diagnostics).
    blwup(:) = lw_upflux_spec(1, :)
    bswup(:) = sw_upflux_spec(1, :)
    bswdn(:) = sw_dnflux_spec(1, :)

  end subroutine aerad_once

  ! -----------------------------------------------------------------------

  subroutine zenith_quad_nodes(n, nout, mu, wq)
  ! Gauss-Legendre quadrature nodes/weights for INT_0^1 f(mu) dmu, used to form
  ! the hemispheric (Bond) shortwave average.  Standard [-1,1] nodes x_i and
  ! weights w_i are mapped to [0,1] by  mu = (x+1)/2,  wq = w/2  (so sum wq = 1).
  ! Only the tabulated orders 2,3,4,6,8 are supported (validated upstream by
  ! exocol_config; an unsupported value falls back to 4 here as a safety net).
    integer,  intent(in)  :: n
    integer,  intent(out) :: nout
    real(r8), intent(out) :: mu(MAXQ), wq(MAXQ)

    real(r8) :: x(MAXQ), w(MAXQ)
    integer  :: nn, i

    nn = n
    mu = 0._r8; wq = 0._r8; x = 0._r8; w = 0._r8

    select case (nn)
    case (2)
      x(1:2) = (/ -0.5773502691896257_r8,  0.5773502691896257_r8 /)
      w(1:2) = (/  1.0000000000000000_r8,  1.0000000000000000_r8 /)
    case (3)
      x(1:3) = (/ -0.7745966692414834_r8,  0.0000000000000000_r8,  0.7745966692414834_r8 /)
      w(1:3) = (/  0.5555555555555556_r8,  0.8888888888888889_r8,  0.5555555555555556_r8 /)
    case (4)
      x(1:4) = (/ -0.8611363115940526_r8, -0.3399810435848563_r8, &
                   0.3399810435848563_r8,  0.8611363115940526_r8 /)
      w(1:4) = (/  0.3478548451374538_r8,  0.6521451548625461_r8, &
                   0.6521451548625461_r8,  0.3478548451374538_r8 /)
    case (6)
      x(1:6) = (/ -0.9324695142031521_r8, -0.6612093864662645_r8, -0.2386191860831969_r8, &
                   0.2386191860831969_r8,  0.6612093864662645_r8,  0.9324695142031521_r8 /)
      w(1:6) = (/  0.1713244923791704_r8,  0.3607615730481386_r8,  0.4679139345726910_r8, &
                   0.4679139345726910_r8,  0.3607615730481386_r8,  0.1713244923791704_r8 /)
    case (8)
      x(1:8) = (/ -0.9602898564975363_r8, -0.7966664774136267_r8, -0.5255324099163290_r8, &
                  -0.1834346424956498_r8,  0.1834346424956498_r8,  0.5255324099163290_r8, &
                   0.7966664774136267_r8,  0.9602898564975363_r8 /)
      w(1:8) = (/  0.1012285362903763_r8,  0.2223810344533745_r8,  0.3137066458778873_r8, &
                   0.3626837833783620_r8,  0.3626837833783620_r8,  0.3137066458778873_r8, &
                   0.2223810344533745_r8,  0.1012285362903763_r8 /)
    case default
      ! Should not happen (exocol_config validates), but stay robust.
      nn = 4
      x(1:4) = (/ -0.8611363115940526_r8, -0.3399810435848563_r8, &
                   0.3399810435848563_r8,  0.8611363115940526_r8 /)
      w(1:4) = (/  0.3478548451374538_r8,  0.6521451548625461_r8, &
                   0.6521451548625461_r8,  0.3478548451374538_r8 /)
    end select

    do i = 1, nn
      mu(i) = 0.5_r8 * (x(i) + 1.0_r8)
      wq(i) = 0.5_r8 * w(i)
    end do
    nout = nn

  end subroutine zenith_quad_nodes

end module exocol_radiation
