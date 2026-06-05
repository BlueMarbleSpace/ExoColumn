program test_iapws95
  ! Standalone validation driver for the native IAPWS-95 port.
  ! Prints computed properties for a set of states in a machine-parseable
  ! format; tools/check_iapws95.py recomputes the same states with the
  ! reference `iapws` Python package and reports relative errors.
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use exocol_iapws95
  implicit none

  type(iapws95_props_t) :: p, vap, liq
  real(r8) :: Psat
  logical  :: ok
  integer  :: i

  ! All states chosen unambiguously single-phase (off the saturation dome),
  ! so the reference package returns a defined cv/cp/w for each.
  real(r8) :: rhoT(2,7) = reshape( [ &
        0.3000_r8,   373.15_r8,   &   ! superheated steam (< sat-vap 0.598)
        5.0000_r8,   473.15_r8,   &   ! superheated steam (< 7.86)
       30.0000_r8,   573.15_r8,   &   ! superheated steam (< 46.2)
      998.0000_r8,   300.0_r8,    &   ! compressed liquid (> sat-liq 996.5)
        2.7265_r8,   800.0_r8,    &   ! superheated, ~1 MPa
      150.0000_r8,   700.0_r8,    &   ! supercritical
        0.00805_r8,  300.0_r8 ],  &   ! low-density vapour
      [2,7] )

  real(r8) :: satT(7) = [ 300._r8, 373.15_r8, 473.15_r8, 573.15_r8, &
                          600._r8, 640._r8, 647.0_r8 ]

  real(r8) :: PT(3,4) = reshape( [ &
        1.0e5_r8,  373.15_r8, 0._r8, &   ! 0=vap  (P < Psat(373)=1.013e5)
        1.0e6_r8,  500.0_r8,  0._r8, &   !        (P < Psat(500)=2.64e6)
        1.0e7_r8,  560.0_r8,  1._r8, &   ! 1=liq  (P > Psat(560)=7.1e6)
        1.0e5_r8,  300.0_r8,  1._r8 ], &   ! 1=liq  (P > Psat(300)=3536 Pa)
      [3,4] )

  do i = 1, 7
    p = iapws95_rhoT(rhoT(1,i), rhoT(2,i))
    write(*,'(a,12es22.13)') 'RHOT ', p%rho, p%T, p%P, p%s, p%h, p%u, &
         p%cv, p%cp, p%w, p%alfav, p%Z
  end do

  do i = 1, 7
    call iapws95_sat(satT(i), Psat, vap, liq, ok)
    write(*,'(a,l2,8es22.13)') 'SAT ', ok, satT(i), Psat, vap%rho, liq%rho, &
         vap%s, liq%s, (vap%h-liq%h)
  end do

  do i = 1, 4
    if (PT(3,i) < 0.5_r8) then
      p = iapws95_PT(PT(1,i), PT(2,i), 'vap', ok)
    else
      p = iapws95_PT(PT(1,i), PT(2,i), 'liq', ok)
    end if
    write(*,'(a,l2,4es22.13)') 'PT ', ok, PT(1,i), PT(2,i), p%rho, p%P
  end do

end program test_iapws95
