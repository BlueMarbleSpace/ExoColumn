module exocol_ozone
! Ozone profiles for ExoColumn.
!
! Provides set_earth_o3_profile(pmid, mwdry, o3mmr) which interpolates a
! tabulated mid-latitude climatological ozone profile onto the model pressure
! grid and returns mass mixing ratios.
!
! The table is derived from the Anderson et al. (1986) / SPARC mid-latitude
! reference atmosphere (15 levels, 1000–0.1 hPa).  Interpolation is linear
! in log(p); values clamp at table boundaries.
!
! This module has no dependency on exocol_config; the caller is responsible
! for branching on exocol_config::o3_profile.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver

  implicit none
  private
  public :: set_earth_o3_profile

  ! Local copy of O3 molecular weight [g/mol] — matches exocol_config::MW_O3.
  real(r8), parameter :: MW_O3_LOC = 47.998_r8

  ! Mid-latitude climatological ozone profile.
  ! TAB_P decreases (surface → TOA), matching ExoColumn's index convention.
  integer,  parameter :: NTAB = 15
  real(r8), parameter :: TAB_P(NTAB) = &    ! Pa
    (/ 1.0e5_r8, 5.0e4_r8, 3.0e4_r8, 1.5e4_r8, 1.0e4_r8, &
       7.0e3_r8, 5.0e3_r8, 3.0e3_r8, 2.0e3_r8, 1.0e3_r8, &
       5.0e2_r8, 2.0e2_r8, 1.0e2_r8, 5.0e1_r8, 1.0e1_r8 /)
  real(r8), parameter :: TAB_VMR(NTAB) = &  ! mol/mol
    (/ 2.5e-8_r8, 4.0e-8_r8, 7.0e-8_r8, 1.5e-7_r8, 3.5e-7_r8, &
       7.5e-7_r8, 2.0e-6_r8, 5.5e-6_r8, 8.0e-6_r8, 7.5e-6_r8, &
       6.0e-6_r8, 3.5e-6_r8, 2.5e-6_r8, 1.5e-6_r8, 8.0e-7_r8 /)

contains

  subroutine set_earth_o3_profile(pmid, mwdry, o3mmr)
  ! Fill o3mmr(pver) with the mid-latitude ozone climatology interpolated onto
  ! pmid (Pa), converted to mass mixing ratio using mwdry [g/mol].
    real(r8), intent(in)  :: pmid(pver)
    real(r8), intent(in)  :: mwdry
    real(r8), intent(out) :: o3mmr(pver)

    real(r8) :: log_tab_p(NTAB), log_p, vmr, w
    integer  :: k, j

    log_tab_p(:) = log(TAB_P(:))  ! decreasing order

    do k = 1, pver
      log_p = log(pmid(k))

      if (log_p >= log_tab_p(1)) then
        vmr = TAB_VMR(1)
      else if (log_p <= log_tab_p(NTAB)) then
        vmr = TAB_VMR(NTAB)
      else
        ! Find j such that log_tab_p(j) >= log_p > log_tab_p(j+1).
        j = 1
        do while (j < NTAB - 1 .and. log_tab_p(j+1) > log_p)
          j = j + 1
        end do
        w   = (log_p - log_tab_p(j)) / (log_tab_p(j+1) - log_tab_p(j))
        vmr = TAB_VMR(j) + w * (TAB_VMR(j+1) - TAB_VMR(j))
      end if

      o3mmr(k) = vmr * MW_O3_LOC / mwdry
    end do

  end subroutine set_earth_o3_profile

end module exocol_ozone
