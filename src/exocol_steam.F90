module exocol_steam
!=======================================================================
! Kasting (1988) Appendix-A general non-ideal moist pseudoadiabat.
!
! Kasting, J. F. (1988), "Runaway and Moist Greenhouse Atmospheres and the
! Evolution of Earth and Venus", Icarus 74, 472-494 — Appendix A.
!
! Provides the adiabatic lapse rate dlnT/dlnP for a parcel of non-ideal
! condensable water vapour mixed with an ideal non-condensable background
! gas, in two regimes:
!
!   * SATURATED (condensation) region — eqs. A4 (d ln alpha_v/d ln T) then
!     A5 (d ln P/d ln T).  Carries the heat capacity of BOTH components,
!     the real-gas vapour/condensate entropies, and the full Pv/P weighting,
!     so it recovers the pure-steam saturated adiabat as Pn->0 and the
!     dilute textbook moist adiabat as alpha_v->0.  This is the regime that
!     governs the moist-greenhouse inner edge of the habitable zone.
!
!   * UNSATURATED (constant-composition, "runaway base") region — eqs. A11
!     (component partial pressures) and A12 ((d ln P/d ln T)_s = Cp/[P(dV/dT)p]),
!     with the non-ideal steam term evaluated at the vapour partial state.
!
! All real-water thermodynamics (Psat, rho_v, s_v, s_c, cp_v, expansivity,
! and hence beta = 1/Z implicitly) come from the native IAPWS-95 port
! (exocol_iapws95); the saturation-curve T-derivatives are obtained by
! central finite differences of iapws95_sat(T +/- dT) — the modern stand-in
! for Kasting's 5-degree differencing of the NBS/NRC steam tables.
!
! Units are SI throughout: T [K], P,Pv,Pn [Pa], rho [kg/m3], s [J/kg/K],
! cp/cv [J/kg/K], gas constants [J/kg/K].  The caller passes the
! non-condensable specific gas constant Rd (= R_universal/m_n) and its
! isobaric specific heat cpdry; the non-condensable isochoric heat is
! Cvn = cpdry - Rd.
!=======================================================================
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use exocol_iapws95, only: iapws95_props_t, iapws95_sat, iapws95_PT, &
                            iapws95_psat_aux, IAPWS_TC, IAPWS_R
  implicit none
  private

  public :: steam_esat            ! Psat(T) [Pa]  (fast Wagner-Pruss auxiliary)
  public :: steam_sat_derivs      ! sat-curve values + dln()/dlnT derivatives
  public :: steam_dlnTdlnP_sat    ! Kasting A4/A5   saturated moist pseudoadiabat
  public :: steam_dlnTdlnP_dry    ! Kasting A11/A12 unsaturated adiabat

contains

  !---------------------------------------------------------------------
  ! Saturation pressure [Pa] from the Wagner-Pruss auxiliary equation —
  ! used for the 'steam' esat option and for capping the dry-region solve.
  !---------------------------------------------------------------------
  pure function steam_esat(T) result(es)
    real(r8), intent(in) :: T
    real(r8) :: es
    es = iapws95_psat_aux(T)
  end function steam_esat

  !---------------------------------------------------------------------
  ! Saturation-curve values at T and their logarithmic T-derivatives.
  ! dsv,dsc are d s/d ln T (NOT d ln s/d ln T); dlnPsat,dlnrhov are
  ! d ln()/d ln T.  Central differences, clamped below the critical point.
  !---------------------------------------------------------------------
  subroutine steam_sat_derivs(T, Psat, rhov, sv, sc, dlnPsat, dlnrhov, dsv, dsc)
    real(r8), intent(in)  :: T
    real(r8), intent(out) :: Psat, rhov, sv, sc          ! values at T
    real(r8), intent(out) :: dlnPsat, dlnrhov, dsv, dsc  ! d ln()/d ln T (s: d/dlnT)
    real(r8) :: dT, Tlo, Thi, denom
    real(r8) :: Ps0, Plo, Phi
    type(iapws95_props_t) :: v0, l0, vlo, llo, vhi, lhi
    logical  :: ok

    call iapws95_sat(T, Ps0, v0, l0, ok)
    Psat = Ps0; rhov = v0%rho; sv = v0%s; sc = l0%s

    dT  = max(1.e-2_r8, 1.e-3_r8*T)
    Thi = min(T + dT, IAPWS_TC - 1.e-3_r8)
    Tlo = max(T - dT, 50._r8)
    if (Thi <= Tlo) then           ! degenerate (T essentially at Tc): no slope
      dlnPsat = 0._r8; dlnrhov = 0._r8; dsv = 0._r8; dsc = 0._r8
      return
    end if

    call iapws95_sat(Tlo, Plo, vlo, llo, ok)
    call iapws95_sat(Thi, Phi, vhi, lhi, ok)

    denom   = log(Thi) - log(Tlo)
    dlnPsat = (log(Phi)      - log(Plo)     ) / denom
    dlnrhov = (log(vhi%rho)  - log(vlo%rho) ) / denom
    dsv     = (vhi%s         - vlo%s        ) / denom
    dsc     = (lhi%s         - llo%s        ) / denom
  end subroutine steam_sat_derivs

  !---------------------------------------------------------------------
  ! Saturated non-ideal moist pseudoadiabat: Kasting (1988) A4 -> A5.
  ! Returns dlnT/dlnP = 1 / (dlnP/dlnT).
  !---------------------------------------------------------------------
  function steam_dlnTdlnP_sat(T, P, Rd, cpdry) result(dlnTdlnP)
    real(r8), intent(in) :: T       ! temperature           [K]
    real(r8), intent(in) :: P       ! total pressure        [Pa]
    real(r8), intent(in) :: Rd      ! non-cond. gas const   [J/kg/K]
    real(r8), intent(in) :: cpdry   ! non-cond. cp          [J/kg/K]
    real(r8) :: dlnTdlnP
    real(r8) :: Psat, rhov, sv, sc, dlnPsat, dlnrhov, dsv, dsc
    real(r8) :: Pv, Pn, rhon, alfav, Cvn, num, den, dlnalfav, dlnPn, dlnP

    call steam_sat_derivs(T, Psat, rhov, sv, sc, dlnPsat, dlnrhov, dsv, dsc)

    Pv    = Psat
    Pn    = max(P - Pv, 1.e-6_r8*P)        ! >0; ->0 in the pure-steam limit
    rhon  = Pn / (Rd*T)                    ! ideal non-condensable density
    alfav = rhov / rhon                    ! A3: rho_v/rho_n (mass ratio)
    Cvn   = cpdry - Rd

    ! A4 (alpha_c = 0, pseudoadiabatic):
    num      = Rd*dlnrhov - Cvn - alfav*dsv
    den      = alfav*(sv - sc) + Rd
    dlnalfav = num / den

    ! A5 (rigorous form): dlnPn/dlnT = 1 + dlnrhov/dlnT - dlnalfav/dlnT
    dlnPn    = 1._r8 + dlnrhov - dlnalfav
    dlnP     = (Pv/P)*dlnPsat + (Pn/P)*dlnPn

    dlnTdlnP = 1._r8 / dlnP
  end function steam_dlnTdlnP_sat

  !---------------------------------------------------------------------
  ! Unsaturated (constant-composition) non-ideal adiabat: Kasting A11/A12.
  ! alfav_fixed = rho_v/rho_n is held constant through the dry layer (set at
  ! the surface / condensation level).  Returns dlnT/dlnP.
  !---------------------------------------------------------------------
  function steam_dlnTdlnP_dry(T, P, alfav_fixed, Rd, cpdry) result(dlnTdlnP)
    real(r8), intent(in) :: T, P, alfav_fixed, Rd, cpdry
    real(r8) :: dlnTdlnP
    real(r8) :: Pv, Pn, Plo, Phi, Pmid, gmid, Pcap
    real(r8) :: cpv, rhovv, alfav_exp, dlnPv_s, dlnPn_s, dlnP
    type(iapws95_props_t) :: st
    integer  :: it
    logical  :: ok

    ! Partial vapour pressure from the fixed composition:
    !   rho_v^EOS(Pv,T) = alfav_fixed * rho_n = alfav_fixed*(P-Pv)/(Rd*T)
    ! Monotone in Pv on the unsaturated branch -> bisection in (0, Pcap].
    Pcap = 0.999_r8*P
    if (T < IAPWS_TC) Pcap = min(Pcap, 0.999_r8*iapws95_psat_aux(T))
    Plo = 1.e-3_r8; Phi = Pcap
    do it = 1, 100
      Pmid = 0.5_r8*(Plo + Phi)
      st   = iapws95_PT(Pmid, T, 'vap', ok)
      gmid = st%rho - alfav_fixed*(P - Pmid)/(Rd*T)
      if (gmid > 0._r8) then
        Phi = Pmid
      else
        Plo = Pmid
      end if
      if (Phi - Plo < 1.e-7_r8*P) exit
    end do
    Pv = 0.5_r8*(Plo + Phi)
    Pn = P - Pv

    st        = iapws95_PT(Pv, T, 'vap', ok)
    cpv       = st%cp
    rhovv     = st%rho
    alfav_exp = st%alfav                 ! (1/V)(dV/dT)_p of steam

    ! A12 per component: (dlnP/dlnT)_s = Cp / [P (dV/dT)_p] = cp*rho/(P*alfav)
    dlnPv_s = cpv*rhovv / (Pv*alfav_exp)         ! non-ideal steam
    dlnPn_s = cpdry / Rd                          ! ideal non-condensable
    ! A11: total
    dlnP    = (Pv/P)*dlnPv_s + (Pn/P)*dlnPn_s

    dlnTdlnP = 1._r8 / dlnP
  end function steam_dlnTdlnP_dry

end module exocol_steam
