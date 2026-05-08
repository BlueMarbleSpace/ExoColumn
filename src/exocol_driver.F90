PROGRAM exocol_driver
!-----------------------------------------------------------------------
! ExoColumn top-level driver.
! Initialises ExoRT radiation, reads the input column, runs to RCE,
! and writes the equilibrium state.
!
! Expected invocation (from the run/ directory):
!   ./exocol.exe
!
! Input:  iofiles/exocol_in.nc   (ExoRT RTprofile_in.nc format)
! Output: iofiles/exocol_out.nc  (ExoRT RTprofile_out.nc format)
!
! ExoRT initialisation sequence mirrors ExoRT/source/src.main/main.F90:
!   initialize_kcoeff        (kabs)
!   initialize_solar         (initialize_rad_mod_1D)
!   init_ref                 (exo_init_ref)
!   init_model_specific      (model_specific / exo_model_specific)
!   init_planck              (exo_radiation_mod)
!   initialize_radbuffer     (initialize_rad_mod_1D)
!
! ExoRT's input_profile is NOT called; exocol_io reads the initial column.
!-----------------------------------------------------------------------

  use shr_kind_mod,          only: r8 => shr_kind_r8
  use exoplanet_mod                       ! exo_pver, exo_g, solar_file, etc.
  use exo_radiation_mod,     only: init_planck
  use exo_init_ref,          only: init_ref
  use initialize_rad_mod_1D, only: initialize_kcoeff, initialize_solar, &
                                   initialize_radbuffer
  use exo_model_specific,    only: init_model_specific
  use ppgrid,                only: pver, pverp
  use exocol_mod
  use exocol_io,             only: read_initial_conditions, write_output
  use exocol_rce_loop,       only: run_rce_loop
  use exocol_radiation,      only: exocol_rad_tend   ! final flux retrieval

  implicit none

  character(len=256) :: input_file  = 'iofiles/exocol_in.nc'
  character(len=256) :: output_file = 'iofiles/exocol_out.nc'

  real(r8), dimension(pver)  :: LWHR, SWHR
  real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

  write(*,'(/,a)') '=== ExoColumn: Radiative-Convective Equilibrium Model ==='
  write(*,'(a,i0,a)') '  pver = ', pver, ' layers (from exoplanet_mod::exo_pver)'

  ! ---- 1. Allocate column state ----
  call exocol_init()

  ! ---- 2. Read initial conditions ----
  call read_initial_conditions(trim(input_file))

  ! ---- 3. Propagate mw/cp into ExoRT physconst ----
  call exocol_setgas()

  ! ---- 4. Compute derived pressure arrays ----
  call exocol_update_derived()

  ! ---- 5. Initialise ExoRT radiation ----
  !    Order must match ExoRT's main.F90: kcoeff → solar → ref → model → planck → buffer
  call initialize_kcoeff
  call initialize_solar
  call init_ref
  call init_model_specific
  call init_planck
  call initialize_radbuffer

  write(*,'(a)') '  ExoRT radiation initialised.'
  write(*,'(a,f7.2,a,f5.3)') '  ts = ', ts, ' K    coszrs = ', coszrs

  ! ---- 6. Run RCE loop ----
  call run_rce_loop()

  ! ---- 7. Retrieve final fluxes for output ----
  call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

  ! ---- 8. Write equilibrium state ----
  call write_output(trim(output_file), LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

  ! ---- 9. Clean up ----
  call exocol_finalize()

  write(*,'(/,a,/)') '=== ExoColumn: done ==='

END PROGRAM exocol_driver
