module exocol_mod
! Column state module for ExoColumn.
! Holds all arrays and scalars describing one atmospheric column, and exposes
! exocol_init / exocol_finalize for allocation / deallocation.
! Every other ExoColumn module USEs this module.
!
! Relationship to ExoRT:
!   pver / pverp come from ExoRT's ppgrid (compile-time via exoplanet_mod::exo_pver).
!   mwdry / cpair in ExoRT's physconst module are set by physconst_setgas(); call
!   exocol_setgas() after exocol_init to propagate mwdry_col / cpdry_col there.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pverp

  implicit none
  private
  save

  ! -----------------------------------------------------------------------
  ! Public interface
  ! -----------------------------------------------------------------------
  public :: exocol_init
  public :: exocol_finalize
  public :: exocol_setgas
  public :: exocol_update_derived

  ! -----------------------------------------------------------------------
  ! Public scalars — surface / atmospheric parameters
  ! -----------------------------------------------------------------------
  real(r8), public :: ts        ! surface (skin) temperature [K]
  real(r8), public :: ps        ! surface pressure [Pa]
  real(r8), public :: mwdry_col ! molecular weight of dry air [g/mol]
  real(r8), public :: cpdry_col ! specific heat of dry air [J/kg/K]
  real(r8), public :: coszrs    ! cosine of solar zenith angle
  real(r8), public :: msdist    ! planet–star distance factor (1/eccf; 1.0 = 1 AU)

  ! Surface albedos
  real(r8), public :: asdir     ! SW visible  direct  albedo (0.2–0.7 µm)
  real(r8), public :: asdif     ! SW visible  diffuse albedo
  real(r8), public :: aldir     ! SW near-IR  direct  albedo (0.7–4.0 µm)
  real(r8), public :: aldif     ! SW near-IR  diffuse albedo

  ! -----------------------------------------------------------------------
  ! Public layer-midpoint arrays (pver)
  ! -----------------------------------------------------------------------
  real(r8), allocatable, public :: tmid(:)    ! temperature at midpoints [K]
  real(r8), allocatable, public :: pmid(:)    ! pressure at midpoints [Pa]
  real(r8), allocatable, public :: pdel(:)    ! wet pressure thickness [Pa]
  real(r8), allocatable, public :: pdeldry(:) ! dry pressure thickness [Pa] (derived)

  ! Absorbing gases — all dry mass mixing ratios [kg/kg] unless noted
  real(r8), allocatable, public :: h2ommr(:)  ! H2O specific humidity [kg(wv)/kg(air)]
  real(r8), allocatable, public :: co2mmr(:)  ! CO2 dry mmr
  real(r8), allocatable, public :: ch4mmr(:)  ! CH4 dry mmr
  real(r8), allocatable, public :: c2h6mmr(:) ! C2H6 dry mmr
  real(r8), allocatable, public :: h2mmr(:)   ! H2  dry mmr
  real(r8), allocatable, public :: n2mmr(:)   ! N2  dry mmr
  real(r8), allocatable, public :: o3mmr(:)   ! O3  dry mmr
  real(r8), allocatable, public :: o2mmr(:)   ! O2  dry mmr

  ! Cloud properties (set to zero for clear-sky runs)
  real(r8), allocatable, public :: cicewp(:)  ! in-cloud ice water path [g/m²]
  real(r8), allocatable, public :: cliqwp(:)  ! in-cloud liquid water path [g/m²]
  real(r8), allocatable, public :: cfrc(:)    ! cloud fraction [0–1]
  real(r8), allocatable, public :: rei(:)     ! ice effective radius [µm]
  real(r8), allocatable, public :: rel(:)     ! liquid effective radius [µm]

  ! -----------------------------------------------------------------------
  ! Public interface arrays (pverp = pver+1)
  ! -----------------------------------------------------------------------
  real(r8), allocatable, public :: tint(:)    ! temperature at interfaces [K]
  real(r8), allocatable, public :: pint(:)    ! wet pressure at interfaces [Pa]
  real(r8), allocatable, public :: pintdry(:) ! dry pressure at interfaces [Pa] (derived)
  real(r8), allocatable, public :: zint(:)    ! geometric height at interfaces [m]

contains

  ! -----------------------------------------------------------------------

  subroutine exocol_init()
  ! Allocate all column arrays and zero-initialise them.
  ! Must be called before any other ExoColumn routine.
  ! pver / pverp are compile-time constants from ppgrid / exoplanet_mod.
    integer :: ierr

    ! --- midpoint arrays ---
    allocate(tmid(pver),    stat=ierr); if (ierr /= 0) stop 'exocol_init: tmid'
    allocate(pmid(pver),    stat=ierr); if (ierr /= 0) stop 'exocol_init: pmid'
    allocate(pdel(pver),    stat=ierr); if (ierr /= 0) stop 'exocol_init: pdel'
    allocate(pdeldry(pver), stat=ierr); if (ierr /= 0) stop 'exocol_init: pdeldry'

    allocate(h2ommr(pver),  stat=ierr); if (ierr /= 0) stop 'exocol_init: h2ommr'
    allocate(co2mmr(pver),  stat=ierr); if (ierr /= 0) stop 'exocol_init: co2mmr'
    allocate(ch4mmr(pver),  stat=ierr); if (ierr /= 0) stop 'exocol_init: ch4mmr'
    allocate(c2h6mmr(pver), stat=ierr); if (ierr /= 0) stop 'exocol_init: c2h6mmr'
    allocate(h2mmr(pver),   stat=ierr); if (ierr /= 0) stop 'exocol_init: h2mmr'
    allocate(n2mmr(pver),   stat=ierr); if (ierr /= 0) stop 'exocol_init: n2mmr'
    allocate(o3mmr(pver),   stat=ierr); if (ierr /= 0) stop 'exocol_init: o3mmr'
    allocate(o2mmr(pver),   stat=ierr); if (ierr /= 0) stop 'exocol_init: o2mmr'

    allocate(cicewp(pver),  stat=ierr); if (ierr /= 0) stop 'exocol_init: cicewp'
    allocate(cliqwp(pver),  stat=ierr); if (ierr /= 0) stop 'exocol_init: cliqwp'
    allocate(cfrc(pver),    stat=ierr); if (ierr /= 0) stop 'exocol_init: cfrc'
    allocate(rei(pver),     stat=ierr); if (ierr /= 0) stop 'exocol_init: rei'
    allocate(rel(pver),     stat=ierr); if (ierr /= 0) stop 'exocol_init: rel'

    ! --- interface arrays ---
    allocate(tint(pverp),    stat=ierr); if (ierr /= 0) stop 'exocol_init: tint'
    allocate(pint(pverp),    stat=ierr); if (ierr /= 0) stop 'exocol_init: pint'
    allocate(pintdry(pverp), stat=ierr); if (ierr /= 0) stop 'exocol_init: pintdry'
    allocate(zint(pverp),    stat=ierr); if (ierr /= 0) stop 'exocol_init: zint'

    ! --- zero-initialise everything ---
    tmid    = 0._r8;  pmid    = 0._r8;  pdel    = 0._r8;  pdeldry  = 0._r8
    h2ommr  = 0._r8;  co2mmr  = 0._r8;  ch4mmr  = 0._r8;  c2h6mmr = 0._r8
    h2mmr   = 0._r8;  n2mmr   = 0._r8;  o3mmr   = 0._r8;  o2mmr   = 0._r8
    cicewp  = 0._r8;  cliqwp  = 0._r8;  cfrc    = 0._r8
    rei     = 0._r8;  rel     = 0._r8
    tint    = 0._r8;  pint    = 0._r8;  pintdry = 0._r8;   zint    = 0._r8

    ts        = 0._r8
    ps        = 0._r8
    mwdry_col = 0._r8
    cpdry_col = 0._r8
    coszrs    = 0._r8
    msdist    = 1._r8
    asdir     = 0._r8;  asdif = 0._r8
    aldir     = 0._r8;  aldif = 0._r8

  end subroutine exocol_init

  ! -----------------------------------------------------------------------

  subroutine exocol_finalize()
  ! Deallocate all column arrays.
    if (allocated(tmid))    deallocate(tmid)
    if (allocated(pmid))    deallocate(pmid)
    if (allocated(pdel))    deallocate(pdel)
    if (allocated(pdeldry)) deallocate(pdeldry)

    if (allocated(h2ommr))  deallocate(h2ommr)
    if (allocated(co2mmr))  deallocate(co2mmr)
    if (allocated(ch4mmr))  deallocate(ch4mmr)
    if (allocated(c2h6mmr)) deallocate(c2h6mmr)
    if (allocated(h2mmr))   deallocate(h2mmr)
    if (allocated(n2mmr))   deallocate(n2mmr)
    if (allocated(o3mmr))   deallocate(o3mmr)
    if (allocated(o2mmr))   deallocate(o2mmr)

    if (allocated(cicewp))  deallocate(cicewp)
    if (allocated(cliqwp))  deallocate(cliqwp)
    if (allocated(cfrc))    deallocate(cfrc)
    if (allocated(rei))     deallocate(rei)
    if (allocated(rel))     deallocate(rel)

    if (allocated(tint))    deallocate(tint)
    if (allocated(pint))    deallocate(pint)
    if (allocated(pintdry)) deallocate(pintdry)
    if (allocated(zint))    deallocate(zint)

  end subroutine exocol_finalize

  ! -----------------------------------------------------------------------

  subroutine exocol_setgas()
  ! Propagate mwdry_col / cpdry_col into ExoRT's physconst module so that
  ! aerad_driver picks up the correct gas properties.
  ! Must be called after mwdry_col and cpdry_col have been set (e.g. after
  ! read_initial_conditions).
    use physconst, only: physconst_setgas
    call physconst_setgas(mwdry_col, cpdry_col)
  end subroutine exocol_setgas

  ! -----------------------------------------------------------------------

  subroutine exocol_update_derived()
  ! Recompute pdeldry and pintdry from the current primary-state arrays.
  ! Must be called after any change to pdel, pint, or h2ommr.
  !
  ! Derivation follows ExoRT main.F90:
  !   pdeldry(k)  = pdel(k)  * (1 - h2ommr(k))
  !   pintdry(1)  = pint(1)  * (1 - h2ommr(1))   ! top interface uses k=1 mmr
  !   pintdry(k)  = pint(k)  * (1 - h2ommr(k-1)) ! lower interfaces use layer below
    integer :: k

    pdeldry(:) = pdel(:) * (1._r8 - h2ommr(:))

    pintdry(1) = pint(1) * (1._r8 - h2ommr(1))
    do k = 2, pverp
      pintdry(k) = pint(k) * (1._r8 - h2ommr(k-1))
    end do

  end subroutine exocol_update_derived

end module exocol_mod
