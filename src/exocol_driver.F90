PROGRAM exocol_driver
!-----------------------------------------------------------------------
! ExoColumn top-level driver.
! Initialises ExoRT radiation, sets up the initial column, runs to RCE,
! and writes the equilibrium state.
!
! Expected invocation (from the project root):
!   ./run/exocol.exe
!
! Initialization is selected by &exocol_init::input_file:
!   empty (default) → cold start from &exocol_init + &exocol_composition
!   non-empty       → read ExoRT-format NetCDF (legacy path)
!
! Output: iofiles/exocol_out.nc  (ExoRT RTprofile_out.nc format)
!
! ExoRT initialisation sequence mirrors ExoRT/source/src.main/main.F90:
!   initialize_kcoeff        (kabs)
!   initialize_solar         (initialize_rad_mod_1D)
!   init_ref                 (exo_init_ref)
!   init_model_specific      (model_specific / exo_model_specific)
!   init_planck              (exo_radiation_mod)
!   initialize_radbuffer     (initialize_rad_mod_1D)
!-----------------------------------------------------------------------

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use exoplanet_mod                       ! exo_pver, exo_g, solar_file, etc.
  use exo_radiation_mod,     only: init_planck
  use exo_init_ref,          only: init_ref
  use initialize_rad_mod_1D, only: initialize_kcoeff, initialize_solar, &
                                   initialize_radbuffer
  use radgrid,               only: use_bps_continuum  ! MT_CKD<->BPS toggle
  use exo_model_specific,    only: init_model_specific
  use ppgrid,                only: pver, pverp
  use exocol_mod
  use exocol_config,         only: read_config, apply_composition_overrides, &
                                   cfg_input_file => input_file, &
                                   latent_heat_mode, flux_only, &
                                   esat_formula, h2o_eos, h2o_continuum, &
                                   cold_trap_phase
  use exocol_io,             only: read_initial_conditions, write_output
  use exocol_coldstart,      only: cold_start_init
  use exocol_rce_loop,       only: run_rce_loop
  use exocol_radiation,      only: exocol_rad_tend   ! final flux retrieval
  use exocol_convadj,        only: set_latent_heat_mode, set_esat_mode, &
                                   set_cold_trap_phase

  implicit none

  character(len=256) :: output_file = 'iofiles/exocol_out.nc'

  real(r8), dimension(pver)  :: LWHR, SWHR
  real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

  write(*,'(/,a)') '=== ExoColumn: Radiative-Convective Equilibrium Model ==='
  write(*,'(a,i0,a)') '  pver = ', pver, ' layers (from exoplanet_mod::exo_pver)'

  ! ---- 0. Read runtime configuration (scheme selection, etc.) ----
  call read_config('exocol_config.nml')
  call set_latent_heat_mode(trim(latent_heat_mode) == 'fixed_vap')
  ! Select the saturation-vapour-pressure formula model-wide.  read_config has
  ! already forced esat_formula='steam' when h2o_eos='nonideal'.
  call set_esat_mode(trim(esat_formula) == 'steam')
  ! Sub-freezing saturation phase (cold trap): default 'ice' is bit-identical
  ! and is also CLIMA's convention; 'liquid' is a sensitivity toggle (~2x wetter
  ! cold stratosphere at 200 K).
  call set_cold_trap_phase(trim(cold_trap_phase) == 'liquid')
  ! Select the H2O continuum model in the ExoRT radiation core.  Must be set
  ! before initialize_kcoeff (which conditionally reads the BPS continuum file).
  use_bps_continuum = (trim(h2o_continuum) == 'bps')

  ! ---- 1. Allocate column state ----
  call exocol_init()

  ! ---- 2. Build initial column (file mode or cold start) ----
  if (len_trim(cfg_input_file) == 0) then
    ! cold_start_init sets mwdry_col, cpdry_col and calls exocol_update_derived().
    call cold_start_init()
  else
    call read_initial_conditions(trim(cfg_input_file))
    ! Apply &exocol_composition overrides (no-op if block absent).
    ! Must precede exocol_setgas so updated mwdry/cpdry propagate to ExoRT.
    call apply_composition_overrides()
    ! Refresh dry-pressure arrays after any ps or composition change.
    call exocol_update_derived()
  end if

  ! ---- 3. Propagate mw/cp into ExoRT physconst ----
  call exocol_setgas()

  ! ---- 4. Initialise ExoRT radiation ----
  !    Order must match ExoRT's main.F90: kcoeff → solar → ref → model → planck → buffer
  call initialize_kcoeff
  call initialize_solar
  call init_ref
  call init_model_specific
  call init_planck
  call initialize_radbuffer

  write(*,'(a)') '  ExoRT radiation initialised.'
  write(*,'(a,f7.2,a,f5.3)') '  ts = ', ts, ' K    coszrs = ', coszrs

  ! ---- 5. Run RCE loop (skipped when flux_only = .true.) ----
  if (.not. flux_only) call run_rce_loop()

  ! ---- 6. Retrieve final fluxes for output ----
  call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

  ! ---- 7. Write equilibrium state ----
  call write_output(trim(output_file), LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

  ! ---- 8. Clean up ----
  call exocol_finalize()

  write(*,'(/,a,/)') '=== ExoColumn: done ==='

END PROGRAM exocol_driver
