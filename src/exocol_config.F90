module exocol_config
! Runtime configuration for ExoColumn.
! Read from a Fortran namelist file (default: 'exocol_config.nml' in the
! working directory).  If the file is absent all settings retain their
! defaults.
!
! Namelist &exocol_nml variables:
!   conv_scheme      CHARACTER  'dry'        dry adiabatic
!                               'moist'      moist rh-weighted (default)
!                               'manabe'     fixed 6.5 K/km
!   moisture_scheme  CHARACTER  'prognostic' surface evap + condensation (default)
!                               'fixed_rh'   legacy RH-relaxation closure
!                               'off'        no moisture update
!   wind_speed       REAL       5.0          surface wind speed [m/s]
!   C_D              REAL       1.5e-3       bulk exchange coefficient [-]
!
! Example exocol_config.nml:
!   &exocol_nml
!     conv_scheme     = 'moist'
!     moisture_scheme = 'prognostic'
!     wind_speed      = 5.0
!     C_D             = 1.5e-3
!   /

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private

  character(len=32), public, save :: conv_scheme     = 'moist'
  character(len=32), public, save :: moisture_scheme = 'prognostic'
  real(r8),          public, save :: wind_speed      = 5.0_r8
  real(r8),          public, save :: C_D             = 1.5e-3_r8

  public :: read_config

contains

  subroutine read_config(filename)
  ! Read exocol_config.nml and set module variables.
  ! Silently uses defaults if the file does not exist or cannot be parsed.
    character(len=*), intent(in) :: filename

    namelist /exocol_nml/ conv_scheme, moisture_scheme, wind_speed, C_D

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
    read(unit, nml=exocol_nml, iostat=ios)
    close(unit)

    if (ios /= 0) then
      write(*,'(a)') &
        '  exocol_config: namelist read error — using defaults.'
      call announce()
      return
    end if

    ! Validate conv_scheme
    select case (trim(adjustl(conv_scheme)))
    case ('dry','moist','manabe')
      ! ok
    case default
      write(*,'(3a)') &
        "  WARNING: unknown conv_scheme='", trim(conv_scheme), &
        "' — defaulting to moist"
      conv_scheme = 'moist'
    end select

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

    call announce()

  end subroutine read_config

  subroutine announce()
  ! Print the active configuration so the run log records what was used.
    select case (trim(adjustl(conv_scheme)))
    case ('dry');    write(*,'(a)') '  Convection scheme : dry adiabatic'
    case ('moist');  write(*,'(a)') '  Convection scheme : moist rh-weighted'
    case ('manabe'); write(*,'(a)') '  Convection scheme : Manabe-Wetherald (6.5 K/km)'
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
  end subroutine announce

end module exocol_config
