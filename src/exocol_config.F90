module exocol_config
! Runtime configuration for ExoColumn.
! Read from a Fortran namelist file (default: 'exocol_config.nml' in the
! working directory).  If the file is absent all settings retain their
! defaults.
!
! Two namelist blocks are supported, both optional:
!
! &exocol_nml — runtime physics options
!   conv_scheme      CHARACTER  'sbm'        simplified Betts-Miller (default,
!                                            recommended; Frierson 2007)
!                               'dry'        dry adiabatic
!                               'moist'      moist rh-weighted
!                               'manabe'     fixed 6.5 K/km
!                               'zm'         Zhang-McFarlane soft
!   moisture_scheme  CHARACTER  'prognostic' surface evap + condensation (default)
!                               'fixed_rh'   legacy RH-relaxation closure
!                               'off'        no moisture update
!   o3_profile       CHARACTER  'uniform'    well-mixed scalar from o3_vmr (default)
!                               'earth'      mid-latitude climatological profile
!                               'rcemip'     RCEMIP analytic (Wing 2018 / konrad)
!                               'none'       zero ozone everywhere
!   wind_speed       REAL       5.0          surface wind speed [m/s]
!   C_D              REAL       1.5e-3       bulk exchange coefficient [-]
!   msdist           REAL       1.0          planet-star distance factor [AU];
!                                            TOA stellar flux scales as 1/msdist²
!   dz_slab          REAL       50.0         slab-ocean thickness [m].  Sets the
!                                            heat capacity H_slab = ρ_w·cp_w·dz_slab
!                                            and therefore the radiative time
!                                            constant.  Smaller dz_slab → faster
!                                            equilibration with same final state;
!                                            5 m is sufficient for most 1-D RCE
!                                            experiments.
!
! &exocol_init — cold-start initial conditions (used when input_file is empty)
!   input_file       CHARACTER  ''           path to an ExoRT-format input.nc.
!                                            If empty (default), build the
!                                            initial column from this block;
!                                            otherwise read the file (legacy path).
!   ts               REAL       288.0        initial surface temperature [K]
!   t_strato         REAL       200.0        stratospheric (isothermal) temp [K]
!   p_top            REAL       1.0          TOA pressure [Pa]
!   rh_init          REAL       0.7          tropospheric relative humidity used
!                                            to seed h2ommr (h2ommr = rh·qsat)
!   coszrs           REAL       0.5          cosine of solar zenith angle
!   cpdry            REAL       1004.0       dry-air cp [J/kg/K]
!   asdir,asdif      REAL       0.25         visible direct / diffuse albedo
!   aldir,aldif      REAL       0.25         near-IR direct / diffuse albedo
!
! On cold start the column is built as: a log-spaced pressure grid from p_top
! to ps, a moist adiabat integrated upward from ts (capped at t_strato above
! the tropopause), h2ommr = rh_init · qsat in the troposphere (zero above),
! and dry-gas MMRs set by &exocol_composition with Earth-like fallbacks for
! unspecified species (N2=0.78, O2=0.21, CO2=4e-4, others=0).
!
! &exocol_composition — dry-air composition and surface pressure (all optional)
!   ps               REAL       <0 sentinel  surface pressure override [Pa].
!                                            If set, pmid/pint/pdel are rescaled
!                                            proportionally so the column has the
!                                            same σ-structure at the new ps.
!   co2_vmr,
!   ch4_vmr,
!   c2h6_vmr,
!   h2_vmr,
!   n2_vmr,
!   o3_vmr,
!   o2_vmr           REAL       <0 sentinel  Dry-air volume mixing ratio (mole
!                                            fraction of dry air).  If ≥ 0, the
!                                            gas is broadcast as a well-mixed
!                                            scalar to every layer (input vertical
!                                            profile discarded).  If < 0, the
!                                            input-file profile is kept (rescaled
!                                            by mwdry_in/mwdry_new when other
!                                            gases change mwdry).
!
! H2O is intentionally not in this namelist.  Water vapor is a prognostic
! variable (surface evaporation + convective transport + condensation sink)
! when moisture_scheme='prognostic', or relaxed toward the input RH when
! 'fixed_rh', or frozen at the input value when 'off'.  Set the initial
! moisture profile via the input NetCDF file, not here.
!
! When any *_vmr override is set, mwdry is recomputed from the dry composition:
!   mwdry_new = Σ_dry VMR_i · Mw_i
! (using layer-1 input dry-air values for non-overridden gases).  cpdry is NOT
! recomputed; the input file value is kept.
!
! Example exocol_config.nml:
!   &exocol_nml
!     conv_scheme     = 'moist'
!     moisture_scheme = 'prognostic'
!     wind_speed      = 5.0
!     C_D             = 1.5e-3
!     msdist          = 1.0
!     dz_slab         = 50.0
!   /
!   &exocol_composition
!     ps      = 1.0e5
!     co2_vmr = 4.0e-4
!     n2_vmr  = 0.78
!     o2_vmr  = 0.21
!   /

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  ! ---- &exocol_nml ----
  character(len=32), public, save :: conv_scheme     = 'sbm'
  character(len=32), public, save :: moisture_scheme = 'prognostic'
  character(len=32), public, save :: o3_profile      = 'uniform'
  real(r8),          public, save :: wind_speed      = 5.0_r8
  real(r8),          public, save :: C_D             = 1.5e-3_r8
  real(r8),          public, save :: msdist          = 1.0_r8
  real(r8),          public, save :: dz_slab         = 50.0_r8
  real(r8),          public, save :: tau_conv        = 3600.0_r8  ! ZM/SBM relaxation time [s]
  real(r8),          public, save :: cape_trigger    =    0.0_r8  ! CAPE activation threshold [J/kg]
  real(r8),          public, save :: rh_sbm          =    0.7_r8  ! SBM reference relative humidity [-]

  ! Vertical grid.  Default n_sfc_layers = 0 → pure log-spacing from p_top to ps
  ! (recommended; lowest level ~400 m, giving Earth-like bulk surface fluxes with
  ! the standard C_D and a numerically well-behaved column).
  !
  ! Setting n_sfc_layers > 0 builds a HYBRID grid: that many near-surface layers
  ! are geometrically spaced (bottom thickness dp_sfc_bot [Pa], ratio sfc_stretch
  ! between consecutive layers), with the remaining pver - n_sfc_layers layers
  ! log-spaced above.  This puts the lowest level at ~8.5 m (dp_sfc_bot=200 Pa).
  ! NOTE: without an explicit boundary-layer mixing scheme the thin bottom layer
  ! equilibrates to Ts (collapsing the flux gradient → suppressed LE/SH) and goes
  ! radiatively stiff in a moist column (collapsing the adaptive dt under sbm).
  ! The hybrid grid is therefore reserved for future work that adds BL mixing;
  ! use the log grid (n_sfc_layers = 0) for production runs.
  integer,  public, save :: n_sfc_layers = 0
  real(r8), public, save :: dp_sfc_bot   = 200.0_r8   ! bottom layer dp [Pa] (~8.5 m)
  real(r8), public, save :: sfc_stretch  = 1.5_r8     ! geometric ratio [-]

  ! ---- &exocol_init (cold-start initial conditions; used when input_file = '') ----
  ! Public so exocol_coldstart can USE them with renames.  Names that collide
  ! with exocol_mod state variables (ts, coszrs, asdir, ...) must be renamed
  ! by consumers; see exocol_coldstart for the pattern.
  character(len=256), public, save :: input_file = ''
  real(r8),           public, save :: ts         = 288.0_r8
  real(r8),           public, save :: t_strato   = 200.0_r8
  real(r8),           public, save :: p_top      = 100.0_r8
  real(r8),           public, save :: rh_init    =   0.7_r8
  real(r8),           public, save :: coszrs     =   0.5_r8
  real(r8),           public, save :: cpdry      = 1004.0_r8
  real(r8),           public, save :: asdir      =   0.25_r8
  real(r8),           public, save :: asdif      =   0.25_r8
  real(r8),           public, save :: aldir      =   0.25_r8
  real(r8),           public, save :: aldif      =   0.25_r8

  ! ---- &exocol_composition (sentinel = negative → use input file) ----
  ! Public so exocol_coldstart can USE them.  Note: `ps` collides with
  ! exocol_mod::ps; consumers that USE both modules must rename one side
  ! (see exocol_coldstart for the pattern).  The other names (co2_vmr,
  ! ch4_vmr, ...) do not collide with anything in exocol_mod.
  real(r8), parameter, public :: COMPOSITION_UNSET = -1.0_r8

  real(r8), public, save :: ps       = COMPOSITION_UNSET
  real(r8), public, save :: co2_vmr  = COMPOSITION_UNSET
  real(r8), public, save :: ch4_vmr  = COMPOSITION_UNSET
  real(r8), public, save :: c2h6_vmr = COMPOSITION_UNSET
  real(r8), public, save :: h2_vmr   = COMPOSITION_UNSET
  real(r8), public, save :: n2_vmr   = COMPOSITION_UNSET
  real(r8), public, save :: o3_vmr   = COMPOSITION_UNSET
  real(r8), public, save :: o2_vmr   = COMPOSITION_UNSET

  ! ---- Dry-gas molecular weights [kg/kmol = g/mol] ----
  ! IUPAC atomic weights to 4 sig figs.  H2O is excluded from mwdry by
  ! convention (water vapor is tracked separately as specific humidity).
  real(r8), parameter, public :: MW_CO2  = 44.010_r8
  real(r8), parameter, public :: MW_CH4  = 16.043_r8
  real(r8), parameter, public :: MW_C2H6 = 30.069_r8
  real(r8), parameter, public :: MW_H2   =  2.016_r8
  real(r8), parameter, public :: MW_N2   = 28.013_r8
  real(r8), parameter, public :: MW_O3   = 47.998_r8
  real(r8), parameter, public :: MW_O2   = 31.999_r8

  ! ---- Dry-gas specific heats at ~300 K [J/kg/K] ----
  ! Used to auto-compute the mass-weighted mean cpdry_col from composition.
  ! H2 dominates mixtures at high mixing ratio (cp_H2 ≈ 14x that of N2), so
  ! getting this right matters for H2-rich or CH4-rich exoplanet atmospheres.
  real(r8), parameter, public :: CP_CO2  =  844.0_r8
  real(r8), parameter, public :: CP_CH4  = 2220.0_r8
  real(r8), parameter, public :: CP_C2H6 = 1729.0_r8
  real(r8), parameter, public :: CP_H2   = 14310.0_r8
  real(r8), parameter, public :: CP_N2   = 1039.0_r8
  real(r8), parameter, public :: CP_O3   =  820.0_r8
  real(r8), parameter, public :: CP_O2   =  919.0_r8

  public :: read_config
  public :: apply_composition_overrides

contains

  ! -----------------------------------------------------------------------

  subroutine read_config(filename)
  ! Read both namelists from the same file.  Each block is optional; missing
  ! or unparseable blocks leave their variables at their defaults (sentinels
  ! for the composition block, so the input file is used unchanged).
    character(len=*), intent(in) :: filename

    namelist /exocol_nml/         conv_scheme, moisture_scheme, o3_profile, &
                                  wind_speed, C_D, msdist, dz_slab, &
                                  n_sfc_layers, dp_sfc_bot, sfc_stretch, &
                                  tau_conv, cape_trigger, rh_sbm
    namelist /exocol_init/        input_file, ts, t_strato, p_top, rh_init, &
                                  coszrs, cpdry, asdir, asdif, aldir, aldif
    namelist /exocol_composition/ ps, co2_vmr, ch4_vmr, c2h6_vmr, &
                                  h2_vmr, n2_vmr, o3_vmr, o2_vmr

    integer :: unit, ios
    logical :: exists

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
      write(*,'(3a)') '  exocol_config: ', trim(filename), &
        ' not found — using defaults.'
      call announce()
      return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') '  exocol_config: cannot open config file — using defaults.'
      call announce()
      return
    end if

    ! Both namelist blocks are optional.  A missing block leaves its
    ! variables at their defaults (sentinels for composition, so the input
    ! file is used unchanged).  iostat is intentionally ignored.
    rewind(unit); read(unit, nml=exocol_nml,         iostat=ios)
    rewind(unit); read(unit, nml=exocol_init,        iostat=ios)
    rewind(unit); read(unit, nml=exocol_composition, iostat=ios)
    close(unit)

    ! Validate conv_scheme
    select case (trim(adjustl(conv_scheme)))
    case ('dry','moist','manabe','zm','sbm')
      ! ok
    case default
      write(*,'(3a)') &
        "  WARNING: unknown conv_scheme='", trim(conv_scheme), &
        "' — defaulting to sbm"
      conv_scheme = 'sbm'
    end select

    ! Validate rh_sbm (reference RH for the SBM scheme; physical range (0,1])
    if (rh_sbm <= 0.0_r8 .or. rh_sbm > 1.0_r8) then
      write(*,'(a,f8.4,a)') &
        '  WARNING: rh_sbm must be in (0,1] (got ', rh_sbm, ') — defaulting to 0.7'
      rh_sbm = 0.7_r8
    end if

    ! Validate moisture_scheme
    select case (trim(adjustl(moisture_scheme)))
    case ('prognostic','fixed_rh','off')
      ! ok
    case default
      write(*,'(3a)') &
        "  WARNING: unknown moisture_scheme='", trim(moisture_scheme), &
        "' — defaulting to prognostic"
      moisture_scheme = 'prognostic'
    end select

    ! Validate o3_profile
    select case (trim(adjustl(o3_profile)))
    case ('uniform','earth','rcemip','none')
      ! ok
    case default
      write(*,'(3a)') &
        "  WARNING: unknown o3_profile='", trim(o3_profile), &
        "' — defaulting to uniform"
      o3_profile = 'uniform'
    end select

    ! Validate msdist (must be positive; ExoRT divides by it)
    if (msdist <= 0.0_r8) then
      write(*,'(a,es10.3,a)') &
        '  WARNING: msdist must be > 0 (got ', msdist, ') — defaulting to 1.0'
      msdist = 1.0_r8
    end if

    ! Validate dz_slab (must be positive; sets the slab heat capacity)
    if (dz_slab <= 0.0_r8) then
      write(*,'(a,es10.3,a)') &
        '  WARNING: dz_slab must be > 0 (got ', dz_slab, ') — defaulting to 50.0'
      dz_slab = 50.0_r8
    end if

    ! Validate hybrid grid parameters
    if (n_sfc_layers < 0) then
      write(*,'(a,i0,a)') &
        '  WARNING: n_sfc_layers must be >= 0 (got ', n_sfc_layers, ') — defaulting to 0'
      n_sfc_layers = 0
    end if
    if (dp_sfc_bot <= 0.0_r8) then
      write(*,'(a,es10.3,a)') &
        '  WARNING: dp_sfc_bot must be > 0 (got ', dp_sfc_bot, ') — defaulting to 200.0'
      dp_sfc_bot = 200.0_r8
    end if
    if (n_sfc_layers > 0 .and. sfc_stretch <= 1.0_r8) then
      write(*,'(a,f8.4,a)') &
        '  WARNING: sfc_stretch must be > 1 when n_sfc_layers > 0 (got ', &
        sfc_stretch, ') — defaulting to 1.5'
      sfc_stretch = 1.5_r8
    end if

    call announce()

  end subroutine read_config

  ! -----------------------------------------------------------------------

  subroutine announce()
  ! Print the active runtime configuration.
    select case (trim(adjustl(conv_scheme)))
    case ('dry');    write(*,'(a)') '  Convection scheme : dry adiabatic'
    case ('moist');  write(*,'(a)') '  Convection scheme : moist rh-weighted'
    case ('manabe'); write(*,'(a)') '  Convection scheme : Manabe-Wetherald (6.5 K/km)'
    case ('zm')
      write(*,'(a)') '  Convection scheme : Zhang-McFarlane soft (CAPE-relaxation)'
      write(*,'(a,f7.1,a)') '    tau_conv     = ', tau_conv,     ' s'
      write(*,'(a,f7.1,a)') '    cape_trigger = ', cape_trigger, ' J/kg'
    case ('sbm')
      write(*,'(a)') '  Convection scheme : simplified Betts-Miller (Frierson 2007)'
      write(*,'(a,f7.1,a)') '    tau_conv (relax) = ', tau_conv, ' s'
      write(*,'(a,f6.3)')   '    rh_sbm (ref RH)  = ', rh_sbm
    end select

    select case (trim(adjustl(moisture_scheme)))
    case ('prognostic')
      write(*,'(a)') '  Moisture scheme   : prognostic (evap + condensation)'
      write(*,'(a,f5.2,a)') '    wind_speed = ', wind_speed, ' m/s'
      write(*,'(a,es8.2)')  '    C_D        = ', C_D
    case ('fixed_rh')
      write(*,'(a)') '  Moisture scheme   : fixed-RH (legacy)'
    case ('off')
      write(*,'(a)') '  Moisture scheme   : off (h2ommr frozen at input)'
    end select

    select case (trim(adjustl(o3_profile)))
    case ('uniform'); write(*,'(a)') '  Ozone profile     : uniform (from o3_vmr)'
    case ('earth');   write(*,'(a)') '  Ozone profile     : Earth mid-latitude climatology'
    case ('none');    write(*,'(a)') '  Ozone profile     : none (zero ozone)'
    end select
    write(*,'(a,f6.3,a)') '  Star distance     : msdist = ', msdist, ' AU'
    write(*,'(a,f6.2,a)') '  Slab thickness    : dz_slab = ', dz_slab, ' m'
    if (n_sfc_layers > 0) then
      write(*,'(a,i0,a,f6.1,a,f5.3)') &
        '  Vertical grid     : hybrid — ', n_sfc_layers, &
        ' surface layers, dp_bot = ', dp_sfc_bot, ' Pa, stretch = ', sfc_stretch
    else
      write(*,'(a)') '  Vertical grid     : log-spaced (n_sfc_layers = 0)'
    end if

    if (len_trim(input_file) == 0) then
      write(*,'(a)')         '  Initialization    : cold start (from &exocol_init)'
      write(*,'(a,f6.1,a)')  '    ts       = ', ts,       ' K'
      write(*,'(a,f6.1,a)')  '    t_strato = ', t_strato, ' K'
      write(*,'(a,es8.2,a)') '    p_top    = ', p_top,    ' Pa'
      write(*,'(a,f5.2)')    '    rh_init  = ', rh_init
      write(*,'(a,f5.2)')    '    coszrs   = ', coszrs
      write(*,'(a,f7.1,a)')  '    cpdry    = ', cpdry,    ' J/kg/K'
      write(*,'(a,4(f5.2,1x))') &
        '    albedos  (asdir asdif aldir aldif) = ', asdir, asdif, aldir, aldif
    else
      write(*,'(2a)')        '  Initialization    : read from file ', trim(input_file)
    end if

    ! Report only the active composition overrides (saves log clutter when
    ! the &exocol_composition block is omitted entirely).
    if (any_composition_override()) then
      write(*,'(a)') '  Composition override (from &exocol_composition):'
      if (ps       >= 0.0_r8) write(*,'(a,es10.3,a)') '    ps       = ', ps,       ' Pa'
      if (co2_vmr  >= 0.0_r8) write(*,'(a,es10.3)')   '    co2_vmr  = ', co2_vmr
      if (ch4_vmr  >= 0.0_r8) write(*,'(a,es10.3)')   '    ch4_vmr  = ', ch4_vmr
      if (c2h6_vmr >= 0.0_r8) write(*,'(a,es10.3)')   '    c2h6_vmr = ', c2h6_vmr
      if (h2_vmr   >= 0.0_r8) write(*,'(a,es10.3)')   '    h2_vmr   = ', h2_vmr
      if (n2_vmr   >= 0.0_r8) write(*,'(a,es10.3)')   '    n2_vmr   = ', n2_vmr
      if (o3_vmr   >= 0.0_r8) write(*,'(a,es10.3)')   '    o3_vmr   = ', o3_vmr
      if (o2_vmr   >= 0.0_r8) write(*,'(a,es10.3)')   '    o2_vmr   = ', o2_vmr
    end if
  end subroutine announce

  ! -----------------------------------------------------------------------

  pure function any_composition_override() result(flag)
  ! True if any &exocol_composition override has been set (≥ 0).
    logical :: flag
    flag = (ps       >= 0.0_r8) .or. &
           (co2_vmr  >= 0.0_r8) .or. (ch4_vmr >= 0.0_r8) .or. &
           (c2h6_vmr >= 0.0_r8) .or. (h2_vmr  >= 0.0_r8) .or. &
           (n2_vmr   >= 0.0_r8) .or. (o3_vmr  >= 0.0_r8) .or. &
           (o2_vmr   >= 0.0_r8)
  end function any_composition_override

  ! -----------------------------------------------------------------------

  subroutine apply_composition_overrides()
  ! Apply &exocol_composition overrides to the column state owned by
  ! exocol_mod.  Must be called AFTER exocol_io::read_initial_conditions
  ! (so the input-file values are loaded) and BEFORE exocol_setgas (so the
  ! updated mwdry propagates to ExoRT's physconst).
  !
  ! Pressure override:
  !   ps_new = ps (override)
  !   scale  = ps_new / ps_in
  !   pmid, pint, pdel all multiplied by scale.
  !   exocol_update_derived() must be called after this to refresh pdeldry/pintdry.
  !
  ! Dry-gas composition override:
  !   For each dry gas, target VMR is the override value (if ≥ 0) or the input
  !   layer-1 VMR (else).  New mwdry = Σ_dry target_VMR · Mw_i.  If target VMRs
  !   sum to substantially less or more than 1, a warning is emitted (we use
  !   the values as-is).  For overridden gases the new MMR scalar is broadcast
  !   to all layers.  For non-overridden gases the input layer profile is
  !   kept but rescaled by (mwdry_in / mwdry_new) to preserve mass when mwdry
  !   changes.
  !
  ! H2O (h2ommr) is intentionally not touched: it is a prognostic variable
  ! handled by the moisture scheme.
    use exocol_mod,   only: ps_col      => ps,        &
                            mwdry_col, cpdry_col, pmid, pint, pdel, &
                            co2mmr, ch4mmr, c2h6mmr,    &
                            h2mmr, n2mmr, o3mmr, o2mmr
    use exocol_ozone, only: set_earth_o3_profile, set_rcemip_o3_profile

    integer, parameter :: NGAS = 7
    integer, parameter :: ICO2=1, ICH4=2, IC2H6=3, &
                          IH2=4, IN2=5, IO3=6, IO2=7

    real(r8) :: Mw(NGAS), vmr_in(NGAS), vmr_tgt(NGAS), override_in(NGAS)
    real(r8) :: mmr1_in(NGAS), mmr_new_scalar(NGAS)
    real(r8) :: mwdry_in, mwdry_new, vmr_sum, scale
    real(r8) :: cp_gas(NGAS), cpdry_new
    integer  :: i
    logical  :: any_gas_override, gas_set(NGAS)
    character(len=4) :: gas_name(NGAS)

    gas_name(ICO2)  = 'CO2'
    gas_name(ICH4)  = 'CH4'
    gas_name(IC2H6) = 'C2H6'
    gas_name(IH2)   = 'H2'
    gas_name(IN2)   = 'N2'
    gas_name(IO3)   = 'O3'
    gas_name(IO2)   = 'O2'

    Mw(ICO2)  = MW_CO2
    Mw(ICH4)  = MW_CH4
    Mw(IC2H6) = MW_C2H6
    Mw(IH2)   = MW_H2
    Mw(IN2)   = MW_N2
    Mw(IO3)   = MW_O3
    Mw(IO2)   = MW_O2

    override_in(ICO2)  = co2_vmr
    override_in(ICH4)  = ch4_vmr
    override_in(IC2H6) = c2h6_vmr
    override_in(IH2)   = h2_vmr
    override_in(IN2)   = n2_vmr
    override_in(IO3)   = o3_vmr
    override_in(IO2)   = o2_vmr

    gas_set(:) = (override_in(:) >= 0.0_r8)
    any_gas_override = any(gas_set)

    ! ---- ps override (independent of gas overrides) ----
    if (ps >= 0.0_r8) then
      if (ps_col <= 0.0_r8) then
        write(*,'(a)') '  exocol_config: input ps <= 0; ps override skipped.'
      else
        scale = ps / ps_col
        write(*,'(a,es10.3,a,es10.3,a,f7.4)') &
          '  Composition: ps override ', ps_col, ' Pa → ', ps, ' Pa (scale = ', scale, ')'
        pmid(:) = pmid(:) * scale
        pint(:) = pint(:) * scale
        pdel(:) = pdel(:) * scale
        ps_col  = ps
      end if
    end if

    ! ---- Dry-gas composition override ----
    if (.not. any_gas_override) return

    mwdry_in = mwdry_col
    if (mwdry_in <= 0.0_r8) then
      write(*,'(a)') &
        '  exocol_config: input mwdry <= 0; gas-composition override skipped.'
      return
    end if

    ! Capture input layer-1 dry-gas MMRs and convert to VMRs (the "well-mixed"
    ! representative used to compute the new mwdry).  H2O is excluded.
    mmr1_in(ICO2)  = co2mmr(1)
    mmr1_in(ICH4)  = ch4mmr(1)
    mmr1_in(IC2H6) = c2h6mmr(1)
    mmr1_in(IH2)   = h2mmr(1)
    mmr1_in(IN2)   = n2mmr(1)
    mmr1_in(IO3)   = o3mmr(1)
    mmr1_in(IO2)   = o2mmr(1)

    do i = 1, NGAS
      vmr_in(i) = mmr1_in(i) * mwdry_in / Mw(i)
    end do

    ! Build target VMRs: override where set, input layer-1 otherwise.
    do i = 1, NGAS
      if (gas_set(i)) then
        vmr_tgt(i) = override_in(i)
      else
        vmr_tgt(i) = vmr_in(i)
      end if
    end do

    vmr_sum = sum(vmr_tgt)
    if (abs(vmr_sum - 1.0_r8) > 1.0e-2_r8) then
      write(*,'(a,f8.5,a)') &
        '  WARNING: dry-air VMRs sum to ', vmr_sum, &
        ' (expected ~1.0). Values used as-is; check your namelist.'
    end if

    ! New mean molecular weight of dry air.
    mwdry_new = sum(vmr_tgt * Mw)
    if (mwdry_new <= 0.0_r8) then
      write(*,'(a)') &
        '  exocol_config: computed mwdry_new <= 0; gas-composition override skipped.'
      return
    end if

    ! New MMR scalar for overridden gases.  For non-overridden gases the
    ! layer profile is preserved (modulo the mwdry_in/mwdry_new rescale below).
    do i = 1, NGAS
      mmr_new_scalar(i) = vmr_tgt(i) * Mw(i) / mwdry_new
    end do

    ! Apply: overridden gases get a constant column; non-overridden gases
    ! keep their vertical structure but rescale by mwdry_in / mwdry_new.
    ! h2ommr is intentionally untouched (prognostic moisture variable).
    call apply_gas(co2mmr,  ICO2,  gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    call apply_gas(ch4mmr,  ICH4,  gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    call apply_gas(c2h6mmr, IC2H6, gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    call apply_gas(h2mmr,   IH2,   gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    call apply_gas(n2mmr,   IN2,   gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    select case (trim(o3_profile))
    case ('none')
      o3mmr(:) = 0.0_r8
    case ('earth')
      call set_earth_o3_profile(pmid, mwdry_new, o3mmr)
    case ('rcemip')
      call set_rcemip_o3_profile(pmid, mwdry_new, o3mmr)
    case default  ! 'uniform'
      call apply_gas(o3mmr, IO3, gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
    end select
    call apply_gas(o2mmr,   IO2,   gas_set, mmr_new_scalar, mwdry_in, mwdry_new)

    ! Recompute mass-weighted mean cp from the updated dry composition.
    cp_gas(ICO2)  = CP_CO2;  cp_gas(ICH4)  = CP_CH4;  cp_gas(IC2H6) = CP_C2H6
    cp_gas(IH2)   = CP_H2;   cp_gas(IN2)   = CP_N2
    cp_gas(IO3)   = CP_O3;   cp_gas(IO2)   = CP_O2
    cpdry_new = sum(mmr_new_scalar * cp_gas)

    write(*,'(a,f7.3,a,f7.3,a)') &
      '  Composition: mwdry ', mwdry_in, ' → ', mwdry_new, ' g/mol'
    write(*,'(a,f7.1,a,f7.1,a)') &
      '  Composition: cpdry ', cpdry_col, ' → ', cpdry_new, ' J/kg/K'
    do i = 1, NGAS
      if (gas_set(i)) then
        write(*,'(a,a4,a,es10.3,a,es10.3)') &
          '    ', gas_name(i), '  VMR = ', vmr_tgt(i), &
          '  MMR = ', mmr_new_scalar(i)
      end if
    end do

    mwdry_col = mwdry_new
    cpdry_col = cpdry_new

  end subroutine apply_composition_overrides

  ! -----------------------------------------------------------------------

  subroutine apply_gas(mmr_arr, idx, gas_set, mmr_new_scalar, mwdry_in, mwdry_new)
  ! Helper: for one gas, either broadcast the scalar override across the
  ! column or rescale the input vertical profile.
    real(r8), intent(inout) :: mmr_arr(:)
    integer,  intent(in)    :: idx
    logical,  intent(in)    :: gas_set(:)
    real(r8), intent(in)    :: mmr_new_scalar(:)
    real(r8), intent(in)    :: mwdry_in, mwdry_new

    if (gas_set(idx)) then
      mmr_arr(:) = mmr_new_scalar(idx)
    else
      mmr_arr(:) = mmr_arr(:) * (mwdry_in / mwdry_new)
    end if
  end subroutine apply_gas

end module exocol_config
