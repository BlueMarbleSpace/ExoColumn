program test_steam_adiabat
  ! Standalone driver for the Kasting Appendix-A adiabat evaluator.
  ! Prints dlnT/dlnP for saturated (A4/A5) and unsaturated (A11/A12) states;
  ! tools/check_steam_adiabat.py validates against an independent Python
  ! reimplementation and against Kasting's ideal-gas analogs A7/A8.
  use shr_kind_mod, only: r8 => shr_kind_r8
  use exocol_steam
  implicit none
  real(r8) :: Rd, cpdry, T, P, dl, alfa
  integer  :: i
  real(r8) :: satTP(2,8), dryTP(2,3)

  Rd = 287.0_r8; cpdry = 1004.0_r8        ! N2-dominated background (mwdry 28.97)

  ! saturated branch: P chosen just above Psat(T) so a real non-condensable exists
  satTP = reshape( [ &
        300._r8, 1.013e5_r8,  350._r8, 1.00e5_r8,  400._r8, 3.00e5_r8,  &
        450._r8, 1.00e6_r8,   500._r8, 2.80e6_r8,  550._r8, 6.30e6_r8,  &
        600._r8, 1.30e7_r8,   640._r8, 2.05e7_r8 ], [2,8] )
  do i = 1, 8
    T = satTP(1,i); P = satTP(2,i)
    dl = steam_dlnTdlnP_sat(T, P, Rd, cpdry)
    write(*,'(a,4es22.13)') 'SATAD ', T, P, dl, 1._r8/dl
  end do

  ! unsaturated branch: fixed composition alpha_v = rho_v/rho_n
  alfa = 0.069_r8
  dryTP = reshape( [ 450._r8, 5.0e5_r8,  500._r8, 1.0e6_r8,  600._r8, 5.0e6_r8 ], [2,3] )
  do i = 1, 3
    T = dryTP(1,i); P = dryTP(2,i)
    dl = steam_dlnTdlnP_dry(T, P, alfa, Rd, cpdry)
    write(*,'(a,5es22.13)') 'DRYAD ', T, P, alfa, dl, 1._r8/dl
  end do
end program test_steam_adiabat
