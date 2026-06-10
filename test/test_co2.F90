program test_co2
  ! Standalone validation of exocol_co2 (CO2 saturation curve + cp_CO2(T)).
  !
  ! Build (from project root, after sourcing the Intel environment):
  !   ifx -O2 -o /tmp/test_co2 \
  !     /models/ExoRT/source/src.misc/shr_kind_mod.F90 \
  !     src/exoplanet_mod.F90 \
  !     /models/ExoRT/source/src.misc/shr_const_mod.F90 \
  !     src/exocol_co2.F90 test/test_co2.F90 -module /tmp
  !
  ! Anchors:
  !   psat: triple point (216.592 K, 0.51795 MPa), dry-ice 1-atm sublimation
  !         point (194.686 K, 101.325 kPa), critical point (304.1282 K,
  !         7.3773 MPa) — Span & Wagner (1996).
  !   tsat: round-trips tsat(psat(T)) = T.
  !   cp:   NIST-JANAF ideal-gas heat capacities (J/mol/K): 32.359 @ 200 K,
  !         37.129 @ 298.15 K, 44.627 @ 500 K (harmonic model is good to <0.5%).
  use shr_kind_mod, only: r8 => shr_kind_r8
  use exocol_co2
  implicit none
  real(r8), parameter :: MW = 44.010_r8
  real(r8) :: ref, got
  integer  :: nfail
  nfail = 0

  call check('psat triple point   ', psat_co2(216.592_r8),  0.51795e6_r8, 1e-4_r8)
  call check('psat dry-ice @1 atm ', psat_co2(194.686_r8),  101325._r8,   1e-3_r8)
  call check('psat critical point ', psat_co2(304.1282_r8), 7.3773e6_r8,  1e-6_r8)
  call check('psat 150 K (subl.)  ', psat_co2(150.0_r8),    psat_co2(150.0_r8), 1.0_r8) ! report only

  call check('tsat(psat(140))     ', tsat_co2(psat_co2(140._r8)), 140._r8, 1e-9_r8)
  call check('tsat(psat(216.59))  ', tsat_co2(psat_co2(216.59_r8)), 216.59_r8, 1e-9_r8)
  call check('tsat(psat(280))     ', tsat_co2(psat_co2(280._r8)), 280._r8, 1e-9_r8)
  call check('tsat(9 bar)         ', tsat_co2(9.0e5_r8), tsat_co2(9.0e5_r8), 1.0_r8)   ! report only

  call check('cp 200 K  (JANAF)   ', cp_co2(200.0_r8),   32.359_r8/MW*1000._r8, 6e-3_r8)
  call check('cp 298.15 K (JANAF) ', cp_co2(298.15_r8),  37.129_r8/MW*1000._r8, 6e-3_r8)
  call check('cp 500 K  (JANAF)   ', cp_co2(500.0_r8),   44.627_r8/MW*1000._r8, 6e-3_r8)

  if (nfail == 0) then
    write(*,'(a)') 'test_co2: ALL PASS'
  else
    write(*,'(a,i0,a)') 'test_co2: ', nfail, ' FAILURES'
    stop 1
  end if

contains

  subroutine check(name, val, want, rtol)
    character(len=*), intent(in) :: name
    real(r8),         intent(in) :: val, want, rtol
    real(r8) :: rerr
    rerr = abs(val - want) / max(abs(want), 1e-30_r8)
    if (rerr <= rtol) then
      write(*,'(a,a,es16.8,a,es16.8,a,es9.2)') name, ' PASS  got=', val, &
        '  want=', want, '  rerr=', rerr
    else
      write(*,'(a,a,es16.8,a,es16.8,a,es9.2)') name, ' FAIL  got=', val, &
        '  want=', want, '  rerr=', rerr
      nfail = nfail + 1
    end if
  end subroutine check

end program test_co2
