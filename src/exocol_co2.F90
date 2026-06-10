module exocol_co2
! CO2 condensation thermodynamics and temperature-dependent specific heat,
! used by the outer-HZ (maximum greenhouse, Kopparapu 2013 Sec. 3.3) cold-start
! profile builder in exocol_coldstart:
!
!   psat_co2(T)  CO2 saturation vapour pressure [Pa].  Span & Wagner (1996,
!                J. Phys. Chem. Ref. Data 25, 1509) auxiliary equations: the
!                sublimation-pressure equation (their Eq. 3.12) below the triple
!                point and the vapour-pressure equation (Eq. 3.13) from the
!                triple point to the critical point.  These are the modern,
!                reference-quality replacement for the Fanale et al. (1982)
!                fits used by Kasting (1991) / Kopparapu (2013) — the same
!                modern-EOS substitution we made for water (IAPWS-95 in place
!                of steam tables).  Anchor checks (test/test_co2.F90):
!                  T_t = 216.592 K → 0.51795 MPa   (triple point, exact)
!                  194.686 K       → 1 atm          (dry-ice point, <0.01%)
!                  T_c = 304.1282 K → 7.3773 MPa    (critical point, exact)
!
!   tsat_co2(p)  Inverse of psat_co2: the CO2 frost-point/dew-point temperature
!                at partial pressure p.  Bisection on the monotone psat curve.
!                Used to pin the cold-start adiabat to the CO2 saturation curve
!                where the ascending column would supersaturate (the Kasting
!                1991 CO2-condensation treatment of the upper troposphere).
!
!   cp_co2(T)    Ideal-gas isobaric specific heat of CO2 [J/kg/K] from
!                statistical mechanics: cp = R(7/2 + Σ g_i E(θ_i/T)) with the
!                Einstein function E(x) = x² eˣ/(eˣ−1)² over the three CO2
!                normal modes (ν2 bend doubly degenerate, ν1 symmetric and ν3
!                asymmetric stretch).  Matches NIST-JANAF to <0.5% over
!                150–650 K (anharmonicity is the residual): 735 J/kg/K at
!                200 K vs the fixed 300-K value 844 used by exocol_config —
!                a 15% difference that directly tilts the dry-adiabatic lapse
!                rate of a cold CO2-dominated column.
!
! Self-contained (no dependency on other ExoColumn modules beyond shr kinds/
! constants) so it can be USEd anywhere without circular dependencies.

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_RGAS   ! J/kmol/K

  implicit none
  private

  public :: psat_co2, tsat_co2, cp_co2
  public :: CO2_TT, CO2_PT, CO2_TC, CO2_PC

  ! CO2 phase-boundary anchor points (Span & Wagner 1996)
  real(r8), parameter :: CO2_TT = 216.592_r8     ! triple-point temperature [K]
  real(r8), parameter :: CO2_PT = 0.51795e6_r8   ! triple-point pressure [Pa]
  real(r8), parameter :: CO2_TC = 304.1282_r8    ! critical temperature [K]
  real(r8), parameter :: CO2_PC = 7.3773e6_r8    ! critical pressure [Pa]

  ! Span & Wagner (1996) Eq. 3.13 vapour-pressure coefficients (T_t..T_c):
  !   ln(p/Pc) = (Tc/T) Σ a_i (1−T/Tc)^t_i
  real(r8), parameter :: VP_A(4) = (/ -7.0602087_r8, 1.9391218_r8, &
                                      -1.6463597_r8, -3.2995634_r8 /)
  real(r8), parameter :: VP_T(4) = (/ 1.0_r8, 1.5_r8, 2.0_r8, 4.0_r8 /)

  ! Span & Wagner (1996) Eq. 3.12 sublimation-pressure coefficients (T < T_t):
  !   ln(p/Pt) = (Tt/T) Σ a_i (1−T/Tt)^t_i
  real(r8), parameter :: SUB_A(3) = (/ -14.740846_r8, 2.4327015_r8, &
                                       -5.3061778_r8 /)
  real(r8), parameter :: SUB_T(3) = (/ 1.0_r8, 1.9_r8, 2.9_r8 /)

  ! CO2 vibrational characteristic temperatures θ = hcν/k [K] and degeneracies:
  ! ν2 = 667.4 cm⁻¹ (bend, g=2), ν1 = 1333.0 cm⁻¹ (sym. stretch, g=1),
  ! ν3 = 2349.1 cm⁻¹ (asym. stretch, g=1); c2 = 1.438777 cm·K.
  real(r8), parameter :: CP_THETA(3) = (/ 960.24_r8, 1917.87_r8, 3379.86_r8 /)
  real(r8), parameter :: CP_DEGEN(3) = (/ 2.0_r8, 1.0_r8, 1.0_r8 /)
  real(r8), parameter :: MW_CO2_LOC  = 44.010_r8   ! kg/kmol

  ! Validity floor: below ~80 K psat underflows double precision anyway.
  real(r8), parameter :: T_FLOOR = 80.0_r8

contains

  ! -----------------------------------------------------------------------

  pure function psat_co2(T) result(p)
  ! CO2 saturation vapour pressure [Pa]: sublimation curve below the triple
  ! point, liquid-vapour curve above, continuous at T_t by construction (both
  ! pass through (T_t, P_t) to fit accuracy).  T is clamped to [T_FLOOR, T_c];
  ! at/above T_c the liquid-vapour boundary ends and P_c is returned.
    real(r8), intent(in) :: T
    real(r8) :: p, T_use, th, s
    integer  :: i

    T_use = min(max(T, T_FLOOR), CO2_TC)
    if (T_use < CO2_TT) then
      th = 1.0_r8 - T_use / CO2_TT
      s  = 0.0_r8
      do i = 1, 3
        s = s + SUB_A(i) * th**SUB_T(i)
      end do
      p = CO2_PT * exp(s * CO2_TT / T_use)
    else
      th = 1.0_r8 - T_use / CO2_TC
      s  = 0.0_r8
      do i = 1, 4
        s = s + VP_A(i) * th**VP_T(i)
      end do
      p = CO2_PC * exp(s * CO2_TC / T_use)
    end if
  end function psat_co2

  ! -----------------------------------------------------------------------

  pure function tsat_co2(p) result(T)
  ! CO2 frost-/dew-point temperature [K] at partial pressure p [Pa]: the
  ! inverse of psat_co2 by bisection (psat is strictly monotone in T).
  ! p ≥ P_c returns T_c (supercritical: no saturation boundary); p below
  ! psat(T_FLOOR) returns T_FLOOR.
    real(r8), intent(in) :: p
    real(r8) :: T, Tlo, Thi, Tmd
    integer  :: it

    if (p >= CO2_PC) then
      T = CO2_TC
      return
    end if
    if (p <= psat_co2(T_FLOOR)) then
      T = T_FLOOR
      return
    end if

    Tlo = T_FLOOR
    Thi = CO2_TC
    do it = 1, 60   ! 60 halvings: ΔT ~ 2e-16 K, full double precision
      Tmd = 0.5_r8 * (Tlo + Thi)
      if (psat_co2(Tmd) < p) then
        Tlo = Tmd
      else
        Thi = Tmd
      end if
    end do
    T = 0.5_r8 * (Tlo + Thi)
  end function tsat_co2

  ! -----------------------------------------------------------------------

  pure function cp_co2(T) result(cp)
  ! Ideal-gas isobaric specific heat of CO2 [J/kg/K] at temperature T:
  !   cp/R = 7/2 + Σ_i g_i x_i² e^{x_i} / (e^{x_i} − 1)²,   x_i = θ_i/T
  ! (translation + rotation of a linear molecule + harmonic vibrations).
  ! Evaluated with exp(−x) to stay overflow-safe for small T.
    real(r8), intent(in) :: T
    real(r8) :: cp, x, emx, s
    integer  :: i

    s = 3.5_r8
    do i = 1, 3
      x = CP_THETA(i) / max(T, T_FLOOR)
      if (x < 200.0_r8) then
        emx = exp(-x)
        s   = s + CP_DEGEN(i) * x * x * emx / (1.0_r8 - emx)**2
      end if
    end do
    cp = s * SHR_CONST_RGAS / MW_CO2_LOC
  end function cp_co2

end module exocol_co2
