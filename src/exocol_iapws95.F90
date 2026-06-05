module exocol_iapws95
!=======================================================================
! Native Fortran port of the IAPWS-95 formulation for the thermodynamic
! properties of ordinary water substance (Wagner & Pruss, J. Phys. Chem.
! Ref. Data 31, 387-535, 2002; IAPWS R6-95 rev. 2016).
!
! The formulation expresses the specific Helmholtz free energy as
!   f(rho,T)/(R T) = phi(delta,tau) = phi0(delta,tau) + phir(delta,tau)
! with delta = rho/rhoc, tau = Tc/T.  All thermodynamic properties follow
! from analytic derivatives of phi0 and phir (release Tables 1, 2, 3, 6).
!
! This module is self-contained: it carries the full coefficient tables
! (transcribed verbatim from the canonical release) and implements
!   * phi0 and phir and their delta/tau derivatives (incl. the Gaussian
!     and non-analytic critical terms),
!   * the property assembler iapws95_rhoT(rho,T),
!   * a single-phase density solver iapws95_PT(P,T,phase),
!   * a two-phase saturation solver iapws95_sat(T) (Maxwell construction)
!     using the IAPWS Supp-sat 1992 auxiliary equations as initial guess.
!
! Validated point-by-point against the reference Python implementation
! (`iapws` package) — see tools/check_iapws95.py.
!
! ExoColumn use: supplies the non-ideal steam properties (beta = 1/Z,
! saturated-vapour/liquid entropies, real-gas cp and expansivity) that the
! Kasting (1988) Appendix-A general moist pseudoadiabat needs in the
! steam-dominated inner-HZ regime.  No coupling to ExoColumn state here.
!=======================================================================
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  private

  ! ---- public interface ----
  public :: iapws95_props_t
  public :: iapws95_rhoT          ! (rho [kg/m3], T [K])          -> props
  public :: iapws95_PT            ! (P [Pa], T [K], phase)        -> props
  public :: iapws95_sat           ! (T [K])                       -> sat props
  public :: iapws95_psat_aux      ! (T [K]) Wagner-Pruss aux Psat [Pa] (fast)

  ! ---- public fixed constants of the formulation ----
  real(r8), parameter, public :: IAPWS_TC   = 647.096_r8      ! critical T   [K]
  real(r8), parameter, public :: IAPWS_RHOC = 322.0_r8        ! critical rho [kg/m3]
  real(r8), parameter, public :: IAPWS_PC   = 22.064e6_r8     ! critical P   [Pa]
  real(r8), parameter, public :: IAPWS_TT   = 273.16_r8       ! triple-point T [K]
  real(r8), parameter, public :: IAPWS_R    = 461.51805_r8    ! specific gas const [J/kg/K]
  real(r8), parameter, public :: IAPWS_MWWV = 18.015268_r8    ! molar mass [g/mol]

  ! Derived properties at a single (rho,T) state.
  type :: iapws95_props_t
    real(r8) :: rho      = 0._r8   ! density                 [kg/m3]
    real(r8) :: T        = 0._r8   ! temperature             [K]
    real(r8) :: P        = 0._r8   ! pressure                [Pa]
    real(r8) :: s        = 0._r8   ! specific entropy        [J/kg/K]
    real(r8) :: h        = 0._r8   ! specific enthalpy       [J/kg]
    real(r8) :: u        = 0._r8   ! specific internal energy[J/kg]
    real(r8) :: cv       = 0._r8   ! isochoric heat capacity [J/kg/K]
    real(r8) :: cp       = 0._r8   ! isobaric  heat capacity [J/kg/K]
    real(r8) :: w        = 0._r8   ! speed of sound          [m/s]
    real(r8) :: alfav    = 0._r8   ! isobaric expansivity (1/V)(dV/dT)p [1/K]
    real(r8) :: Z        = 0._r8   ! compressibility P/(rho R T)        [-]
    real(r8) :: dpdrho_T = 0._r8   ! (dP/drho)_T             [Pa m3/kg]
    real(r8) :: dpdT_rho = 0._r8   ! (dP/dT)_rho             [Pa/K]
  end type iapws95_props_t

  !---------------------------------------------------------------------
  ! Ideal-gas part phi0 (release Table 1)
  !---------------------------------------------------------------------
  real(r8), parameter :: n0_1 = -8.3204464837497_r8
  real(r8), parameter :: n0_2 =  6.6832105275932_r8
  real(r8), parameter :: n0_3 =  3.00632_r8
  real(r8), parameter :: n0(5) = [ 0.012436_r8, 0.97315_r8, 1.2795_r8, &
                                   0.96956_r8, 0.24873_r8 ]
  real(r8), parameter :: g0(5) = [ 1.28728967_r8, 3.53734222_r8, 7.74073708_r8, &
                                   9.24437796_r8, 27.5075105_r8 ]

  !---------------------------------------------------------------------
  ! Residual part phir (release Table 2)
  !---------------------------------------------------------------------
  ! Group 1: 7 polynomial terms   n * delta^d * tau^t
  real(r8), parameter :: n1(7) = [ 0.012533547935523_r8, 7.8957634722828_r8, &
        -8.7803203303561_r8, 0.31802509345418_r8, -0.26145533859358_r8, &
        -0.0078199751687981_r8, 0.0088089493102134_r8 ]
  integer,  parameter :: d1(7) = [ 1, 1, 1, 2, 2, 3, 4 ]
  real(r8), parameter :: t1(7) = [ -0.5_r8, 0.875_r8, 1._r8, 0.5_r8, 0.75_r8, &
        0.375_r8, 1._r8 ]

  ! Group 2: 44 exponential terms  n * delta^d * tau^t * exp(-delta^c)
  real(r8), parameter :: n2(44) = [ &
        -0.66856572307965_r8, 0.20433810950965_r8, -6.6212605039687e-05_r8, &
        -0.19232721156002_r8, -0.25709043003438_r8, 0.16074868486251_r8, &
        -0.040092828925807_r8, 3.9343422603254e-07_r8, -7.5941377088144e-06_r8, &
        0.00056250979351888_r8, -1.5608652257135e-05_r8, 1.1537996422951e-09_r8, &
        3.6582165144204e-07_r8, -1.3251180074668e-12_r8, -6.2639586912454e-10_r8, &
        -0.10793600908932_r8, 0.017611491008752_r8, 0.22132295167546_r8, &
        -0.40247669763528_r8, 0.58083399985759_r8, 0.0049969146990806_r8, &
        -0.031358700712549_r8, -0.74315929710341_r8, 0.4780732991548_r8, &
        0.020527940895948_r8, -0.13636435110343_r8, 0.014180634400617_r8, &
        0.0083326504880713_r8, -0.029052336009585_r8, 0.038615085574206_r8, &
        -0.020393486513704_r8, -0.0016554050063734_r8, 0.0019955571979541_r8, &
        0.00015870308324157_r8, -1.638856834253e-05_r8, 0.043613615723811_r8, &
        0.034994005463765_r8, -0.076788197844621_r8, 0.022446277332006_r8, &
        -6.2689710414685e-05_r8, -5.5711118565645e-10_r8, -0.19905718354408_r8, &
        0.31777497330738_r8, -0.11841182425981_r8 ]
  integer,  parameter :: d2(44) = [ 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, &
        13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, &
        3, 4, 4, 5, 14, 3, 6, 6, 6 ]
  real(r8), parameter :: t2(44) = [ 4._r8, 6._r8, 12._r8, 1._r8, 5._r8, 4._r8, &
        2._r8, 13._r8, 9._r8, 3._r8, 4._r8, 11._r8, 4._r8, 13._r8, 1._r8, 7._r8, &
        1._r8, 9._r8, 10._r8, 10._r8, 3._r8, 7._r8, 10._r8, 10._r8, 6._r8, &
        10._r8, 10._r8, 1._r8, 2._r8, 3._r8, 4._r8, 8._r8, 6._r8, 9._r8, 8._r8, &
        16._r8, 22._r8, 23._r8, 23._r8, 10._r8, 50._r8, 44._r8, 46._r8, 50._r8 ]
  integer,  parameter :: c2(44) = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
        3, 4, 6, 6, 6, 6 ]

  ! Group 3: 3 Gaussian terms  n*delta^d*tau^t*exp(-alf*(delta-eps)^2-bet*(tau-gam)^2)
  real(r8), parameter :: n3(3)   = [ -31.306260323435_r8, 31.546140237781_r8, &
        -2521.3154341695_r8 ]
  integer,  parameter :: d3(3)   = [ 3, 3, 3 ]
  real(r8), parameter :: t3(3)   = [ 0._r8, 1._r8, 4._r8 ]
  real(r8), parameter :: alf3(3) = [ 20._r8, 20._r8, 20._r8 ]
  real(r8), parameter :: bet3(3) = [ 150._r8, 150._r8, 250._r8 ]
  real(r8), parameter :: gam3(3) = [ 1.21_r8, 1.21_r8, 1.25_r8 ]
  real(r8), parameter :: eps3(3) = [ 1._r8, 1._r8, 1._r8 ]

  ! Group 4: 2 non-analytic critical terms  n * Delta^b * delta * psi
  real(r8), parameter :: n4(2)   = [ -0.14874640856724_r8, 0.31806110878444_r8 ]
  real(r8), parameter :: a4(2)   = [ 3.5_r8, 3.5_r8 ]        ! exponent 'a' in Delta
  real(r8), parameter :: b4(2)   = [ 0.85_r8, 0.95_r8 ]      ! exponent 'b' on Delta
  real(r8), parameter :: BB4(2)  = [ 0.2_r8, 0.2_r8 ]        ! 'B'
  real(r8), parameter :: CC4(2)  = [ 28._r8, 32._r8 ]        ! 'C'
  real(r8), parameter :: DD4(2)  = [ 700._r8, 800._r8 ]      ! 'D'
  real(r8), parameter :: AA4(2)  = [ 0.32_r8, 0.32_r8 ]      ! 'A'
  real(r8), parameter :: bt4(2)  = [ 0.3_r8, 0.3_r8 ]        ! 'beta'

  !---------------------------------------------------------------------
  ! Supp-sat 1992 auxiliary equations (IAPWS, Rev. Suppl. Release on
  ! Saturation Properties) — used only for initial guesses / fast Psat.
  !---------------------------------------------------------------------
  ! Saturation pressure: ln(Ps/Pc) = (Tc/T) * sum a_i * theta^e_i, theta=1-T/Tc
  real(r8), parameter :: pv_a(6) = [ -7.85951783_r8, 1.84408259_r8, &
        -11.7866497_r8, 22.6807411_r8, -15.9618719_r8, 1.80122502_r8 ]
  real(r8), parameter :: pv_e(6) = [ 1._r8, 1.5_r8, 3._r8, 3.5_r8, 4._r8, 7.5_r8 ]
  ! Saturated-liquid density: rho'/rhoc = 1 + sum b_i * theta^(e_i/3)
  real(r8), parameter :: rl_a(6) = [ 1.99274064_r8, 1.09965342_r8, &
        -0.510839303_r8, -1.75493479_r8, -45.5170352_r8, -674694.45_r8 ]
  real(r8), parameter :: rl_e(6) = [ 1._r8, 2._r8, 5._r8, 16._r8, 43._r8, 110._r8 ]
  ! Saturated-vapour density: rho''/rhoc = exp( sum c_i * theta^(e_i/3) )
  real(r8), parameter :: rg_a(6) = [ -2.0315024_r8, -2.6830294_r8, &
        -5.38626492_r8, -17.2991605_r8, -44.7586581_r8, -63.9201063_r8 ]
  real(r8), parameter :: rg_e(6) = [ 1._r8, 2._r8, 4._r8, 9._r8, 18.5_r8, 35.5_r8 ]

contains

  !=====================================================================
  ! Ideal-gas dimensionless Helmholtz energy and its tau-derivatives.
  ! (delta-derivatives are trivial: phi0_d = 1/delta, phi0_dd = -1/delta^2,
  !  phi0_dt = 0, and are applied directly by the caller.)
  !=====================================================================
  pure subroutine phi0_derivs(tau, delta, fio, fiot, fiott)
    real(r8), intent(in)  :: tau, delta
    real(r8), intent(out) :: fio, fiot, fiott
    integer  :: i
    real(r8) :: e, om
    fio   = log(delta) + n0_1 + n0_2*tau + n0_3*log(tau)
    fiot  = n0_2 + n0_3/tau
    fiott = -n0_3/tau**2
    do i = 1, 5
      e   = exp(-g0(i)*tau)
      om  = 1._r8 - e
      fio   = fio   + n0(i)*log(om)
      fiot  = fiot  + n0(i)*g0(i)*(1._r8/om - 1._r8)
      fiott = fiott - n0(i)*g0(i)**2 * e/om**2
    end do
  end subroutine phi0_derivs

  !=====================================================================
  ! Residual dimensionless Helmholtz energy and all derivatives needed
  ! for the full property set (release Table 6).
  !=====================================================================
  pure subroutine phir_derivs(tau, delta, fir, fird, firdd, firt, firtt, firdt)
    real(r8), intent(in)  :: tau, delta
    real(r8), intent(out) :: fir, fird, firdd, firt, firtt, firdt
    integer  :: i, d, c
    real(r8) :: nn, t, dp, ex, cdc, fac
    real(r8) :: al, be, ga, ep, gd, gt, G
    real(r8) :: dm, tm, dm2, psi, the, Del, Delb
    real(r8) :: psid, psidd, psit, psitt, psidt
    real(r8) :: dDel_d, dDel_dd, dDelb_d, dDelb_dd, dDelb_t, dDelb_tt, dDelb_dt
    real(r8) :: pw1, pw2

    fir = 0._r8; fird = 0._r8; firdd = 0._r8
    firt = 0._r8; firtt = 0._r8; firdt = 0._r8

    ! ---- Group 1: polynomial ----
    do i = 1, 7
      nn = n1(i); d = d1(i); t = t1(i)
      dp = nn * delta**d * tau**t
      fir   = fir   + dp
      fird  = fird  + dp * d / delta
      firdd = firdd + dp * d*(d-1) / delta**2
      firt  = firt  + dp * t / tau
      firtt = firtt + dp * t*(t-1._r8) / tau**2
      firdt = firdt + dp * d*t / (delta*tau)
    end do

    ! ---- Group 2: exponential exp(-delta^c) ----
    do i = 1, 44
      nn = n2(i); d = d2(i); t = t2(i); c = c2(i)
      ex  = exp(-delta**c)
      dp  = nn * delta**d * tau**t * ex
      cdc = real(c,r8) * delta**c           ! c*delta^c
      fac = real(d,r8) - cdc                 ! (d - c*delta^c)
      fir   = fir   + dp
      fird  = fird  + dp * fac / delta
      firdd = firdd + dp * ( fac*(real(d,r8)-1._r8-cdc) - real(c,r8)*cdc ) / delta**2
      firt  = firt  + dp * t / tau
      firtt = firtt + dp * t*(t-1._r8) / tau**2
      firdt = firdt + dp * fac * t / (delta*tau)
    end do

    ! ---- Group 3: Gaussian bell ----
    do i = 1, 3
      nn = n3(i); d = d3(i); t = t3(i)
      al = alf3(i); be = bet3(i); ga = gam3(i); ep = eps3(i)
      G  = exp( -al*(delta-ep)**2 - be*(tau-ga)**2 )
      dp = nn * delta**d * tau**t * G
      gd = real(d,r8)/delta - 2._r8*al*(delta-ep)   ! (d/delta - 2 al (delta-eps))
      gt = t/tau - 2._r8*be*(tau-ga)                ! (t/tau - 2 be (tau-gam))
      fir   = fir   + dp
      fird  = fird  + dp * gd
      firdd = firdd + dp * ( gd*gd - real(d,r8)/delta**2 - 2._r8*al )
      firt  = firt  + dp * gt
      firtt = firtt + dp * ( gt*gt - t/tau**2 - 2._r8*be )
      firdt = firdt + dp * gd * gt
    end do

    ! ---- Group 4: non-analytic critical terms ----
    ! Guard the exact critical point (delta=tau=1) where fractional powers of
    ! (delta-1)^2 diverge; physical states essentially never land there.
    dm = delta - 1._r8
    tm = tau   - 1._r8
    if (abs(dm) < 1.e-12_r8) dm = sign(1.e-12_r8, dm + 1.e-300_r8)
    dm2 = dm*dm
    do i = 1, 2
      nn = n4(i)
      al = CC4(i); be = DD4(i)                ! C, D in psi
      psi   = exp( -al*dm2 - be*tm*tm )
      the   = tm * (-1._r8) + AA4(i) * dm2**(0.5_r8/bt4(i))    ! theta = (1-tau)+A*(dm2)^(1/(2beta))
      ! NB sign: theta = (1 - tau) + A*(...) ; (1-tau) = -tm
      the   = (1._r8 - tau) + AA4(i) * dm2**(0.5_r8/bt4(i))
      Del   = the*the + BB4(i) * dm2**a4(i)                    ! Delta
      if (Del <= 0._r8) Del = 1.e-300_r8
      Delb  = Del**b4(i)                                       ! Delta^b

      ! psi derivatives
      psid  = -2._r8*al*dm * psi
      psidd = ( 2._r8*al*dm2 - 1._r8) * 2._r8*al * psi
      psit  = -2._r8*be*tm * psi
      psitt = ( 2._r8*be*tm*tm - 1._r8) * 2._r8*be * psi
      psidt =  4._r8*al*be*dm*tm * psi

      ! Delta derivatives
      pw1 = dm2**(0.5_r8/bt4(i) - 1._r8)          ! (dm2)^(1/(2b)-1)
      pw2 = dm2**(a4(i) - 1._r8)                  ! (dm2)^(a-1)
      dDel_d  = dm * ( AA4(i)*the*(2._r8/bt4(i))*pw1 + 2._r8*BB4(i)*a4(i)*pw2 )
      dDel_dd = dDel_d/dm + dm2 * ( 4._r8*BB4(i)*a4(i)*(a4(i)-1._r8)*dm2**(a4(i)-2._r8) &
                + 2._r8*AA4(i)**2*(1._r8/bt4(i))**2 * pw1*pw1 &
                + AA4(i)*the*(4._r8/bt4(i))*(0.5_r8/bt4(i)-1._r8)*dm2**(0.5_r8/bt4(i)-2._r8) )

      ! Delta^b derivatives
      dDelb_d  = b4(i)*Del**(b4(i)-1._r8)*dDel_d
      dDelb_dd = b4(i)*( Del**(b4(i)-1._r8)*dDel_dd &
                 + (b4(i)-1._r8)*Del**(b4(i)-2._r8)*dDel_d*dDel_d )
      dDelb_t  = -2._r8*the*b4(i)*Del**(b4(i)-1._r8)
      dDelb_tt =  2._r8*b4(i)*Del**(b4(i)-1._r8) &
                 + 4._r8*the*the*b4(i)*(b4(i)-1._r8)*Del**(b4(i)-2._r8)
      dDelb_dt = -AA4(i)*b4(i)*(2._r8/bt4(i))*Del**(b4(i)-1._r8)*dm*pw1 &
                 - 2._r8*the*b4(i)*(b4(i)-1._r8)*Del**(b4(i)-2._r8)*dDel_d

      fir   = fir   + nn * Delb * delta * psi
      fird  = fird  + nn * ( Delb*(psi + delta*psid) + dDelb_d*delta*psi )
      firdd = firdd + nn * ( Delb*(2._r8*psid + delta*psidd) &
                + 2._r8*dDelb_d*(psi + delta*psid) + dDelb_dd*delta*psi )
      firt  = firt  + nn * delta * ( dDelb_t*psi + Delb*psit )
      firtt = firtt + nn * delta * ( dDelb_tt*psi + 2._r8*dDelb_t*psit + Delb*psitt )
      firdt = firdt + nn * ( Delb*(psit + delta*psidt) + delta*dDelb_d*psit &
                + dDelb_t*(psi + delta*psid) + dDelb_dt*delta*psi )
    end do
  end subroutine phir_derivs

  !=====================================================================
  ! Full property set from a single-phase (rho,T) state.
  !=====================================================================
  pure function iapws95_rhoT(rho, T) result(p)
    real(r8), intent(in) :: rho, T
    type(iapws95_props_t) :: p
    real(r8) :: delta, tau
    real(r8) :: fio, fiot, fiott
    real(r8) :: fir, fird, firdd, firt, firtt, firdt
    real(r8) :: num, den, csum

    delta = rho / IAPWS_RHOC
    tau   = IAPWS_TC / T
    call phi0_derivs(tau, delta, fio, fiot, fiott)
    call phir_derivs(tau, delta, fir, fird, firdd, firt, firtt, firdt)

    p%rho = rho
    p%T   = T
    p%Z   = 1._r8 + delta*fird
    p%P   = p%Z * rho * IAPWS_R * T
    p%s   = IAPWS_R * ( tau*(fiot+firt) - fio - fir )
    p%h   = IAPWS_R * T * ( 1._r8 + tau*(fiot+firt) + delta*fird )
    p%u   = p%h - p%P/rho
    p%cv  = -IAPWS_R * tau**2 * (fiott + firtt)

    num   = (1._r8 + delta*fird - delta*tau*firdt)**2
    den   =  1._r8 + 2._r8*delta*fird + delta**2*firdd
    p%cp  = p%cv + IAPWS_R * num/den

    csum  = tau**2 * (fiott + firtt)          ! < 0
    ! Speed of sound; argument can go negative for mechanically-unstable or
    ! out-of-range (e.g. sub-triple-point metastable) states the saturation
    ! solver may transiently probe.  Guard it so the EOS never raises a floating
    ! invalid under -fpe0 (w is a diagnostic; the pseudoadiabat does not use it).
    p%w   = sqrt( max( IAPWS_R * T * ( den - num/csum ), 0._r8 ) )

    p%alfav    = (1._r8 + delta*fird - delta*tau*firdt) / (T*den)
    p%dpdrho_T = IAPWS_R * T * den
    p%dpdT_rho = rho * IAPWS_R * (1._r8 + delta*fird - delta*tau*firdt)
  end function iapws95_rhoT

  !=====================================================================
  ! Fast Wagner-Pruss auxiliary saturation pressure [Pa].  Good to ~0.01%;
  ! used where only Psat(T) is needed (e.g. the esat 'steam' option).
  !=====================================================================
  pure function iapws95_psat_aux(T) result(Ps)
    real(r8), intent(in) :: T
    real(r8) :: Ps, th, s, Tuse
    integer  :: i
    Tuse = min(max(T, 50._r8), IAPWS_TC)
    th = 1._r8 - Tuse/IAPWS_TC
    s  = 0._r8
    do i = 1, 6
      s = s + pv_a(i) * th**pv_e(i)
    end do
    Ps = IAPWS_PC * exp( (IAPWS_TC/Tuse) * s )
  end function iapws95_psat_aux

  ! Supp-sat auxiliary saturated-liquid density [kg/m3] (initial guess).
  pure function rho_liq_aux(T) result(rho)
    real(r8), intent(in) :: T
    real(r8) :: rho, th, s, Tuse
    integer  :: i
    Tuse = min(max(T, IAPWS_TT), IAPWS_TC)
    th = (1._r8 - Tuse/IAPWS_TC)**(1._r8/3._r8)
    s  = 1._r8
    do i = 1, 6
      s = s + rl_a(i) * th**rl_e(i)
    end do
    rho = s * IAPWS_RHOC
  end function rho_liq_aux

  ! Supp-sat auxiliary saturated-vapour density [kg/m3] (initial guess).
  pure function rho_vap_aux(T) result(rho)
    real(r8), intent(in) :: T
    real(r8) :: rho, th, s, Tuse
    integer  :: i
    Tuse = min(max(T, IAPWS_TT), IAPWS_TC)
    th = (1._r8 - Tuse/IAPWS_TC)**(1._r8/3._r8)
    s  = 0._r8
    do i = 1, 6
      s = s + rg_a(i) * th**rg_e(i)
    end do
    rho = exp(s) * IAPWS_RHOC
  end function rho_vap_aux

  !=====================================================================
  ! Single-phase density solve for (P,T).  phase = 'vap' | 'liq' | 'auto'.
  !
  ! The sub-critical liquid branch is razor-stiff (water is nearly
  ! incompressible: P swings by >1 bar over Δrho ~ 0.05 kg/m3 at the
  ! saturation density) and the single-phase isotherm has a non-monotone
  ! van der Waals loop inside the two-phase dome.  Plain Newton can fall
  ! off the physical branch into the loop, so we set up a monotone bracket
  ! per phase and use a safeguarded Newton/bisection root find (rtsafe).
  !=====================================================================
  function iapws95_PT(P, T, phase, ok) result(pr)
    real(r8),         intent(in)  :: P, T
    character(len=*), intent(in)  :: phase
    logical,          intent(out) :: ok
    type(iapws95_props_t) :: pr
    real(r8) :: rho, x1, x2, rsat
    character(len=4) :: ph

    ph = trim(phase)
    if (T >= IAPWS_TC) ph = 'auto'   ! supercritical: single monotone branch

    if (ph == 'auto') then
      ! pick phase from the auxiliary saturation pressure
      if (T < IAPWS_TC .and. P > iapws95_psat_aux(T)) then
        ph = 'liq'
      else
        ph = 'vap'
      end if
    end if

    if (T >= IAPWS_TC) then
      ! monotone supercritical isotherm: bracket [tiny, very dense]
      x1 = 1.e-8_r8
      x2 = 1300._r8
    else if (ph == 'liq') then
      rsat = rho_liq_aux(T)
      x1   = rsat                 ! at/below true sat-liquid density: P <= Psat < target
      x2   = 1300._r8             ! compressed liquid: P huge
    else  ! vapour
      rsat = rho_vap_aux(T)
      x1   = 1.e-10_r8            ! rho->0: P->0 < target
      x2   = rsat                 ! sat vapour: P = Psat >= target
    end if

    rho = solve_rho_bracket(T, P, x1, x2, ok)
    pr  = iapws95_rhoT(rho, T)
  end function iapws95_PT

  !---------------------------------------------------------------------
  ! Safeguarded Newton/bisection (rtsafe) for P(rho,T)=Ptar on a bracket
  ! [x1,x2] assumed to contain the root with the function monotone there.
  ! Falls back gracefully if the endpoints do not bracket.
  !---------------------------------------------------------------------
  function solve_rho_bracket(T, Ptar, x1, x2, ok) result(rho)
    real(r8), intent(in)  :: T, Ptar, x1, x2
    logical,  intent(out) :: ok
    real(r8) :: rho
    real(r8) :: a, b, fa, fb, fr, dfr, dx, dxold, rtemp
    type(iapws95_props_t) :: pp
    integer  :: it
    integer, parameter :: maxit = 120
    real(r8), parameter :: xtol = 1.e-9_r8

    a = x1; b = x2
    pp = iapws95_rhoT(a, T); fa = pp%P - Ptar
    pp = iapws95_rhoT(b, T); fb = pp%P - Ptar
    ok = .false.

    if (fa == 0._r8) then; rho = a; ok = .true.; return; end if
    if (fb == 0._r8) then; rho = b; ok = .true.; return; end if

    if (fa*fb > 0._r8) then
      ! not bracketed (e.g. target outside the branch range): damped Newton
      ! from the endpoint with the smaller residual.
      rho = merge(a, b, abs(fa) < abs(fb))
      do it = 1, maxit
        pp  = iapws95_rhoT(rho, T)
        dfr = pp%dpdrho_T
        if (dfr <= 0._r8) dfr = IAPWS_R*T
        dx  = (pp%P - Ptar)/dfr
        if (dx >  0.5_r8*rho) dx =  0.5_r8*rho
        if (dx < -1.0_r8*rho) dx = -1.0_r8*rho
        rho = max(rho - dx, 1.e-12_r8)
        if (abs(dx) <= xtol*rho) then; ok = .true.; exit; end if
      end do
      return
    end if

    ! orient so that fa < 0
    if (fa > 0._r8) then
      rtemp = a; a = b; b = rtemp
      rtemp = fa; fa = fb; fb = rtemp
    end if

    rho   = 0.5_r8*(x1 + x2)
    dxold = abs(x2 - x1)
    dx    = dxold
    pp  = iapws95_rhoT(rho, T); fr = pp%P - Ptar; dfr = pp%dpdrho_T

    do it = 1, maxit
      ! bisect if Newton out of range or not decreasing fast enough
      if ( ((rho-b)*dfr-fr)*((rho-a)*dfr-fr) > 0._r8 .or. &
           abs(2._r8*fr) > abs(dxold*dfr) ) then
        dxold = dx
        dx    = 0.5_r8*(b - a)
        rho   = a + dx
      else
        dxold = dx
        dx    = fr/dfr
        rho   = rho - dx
      end if
      if (abs(dx) <= xtol*rho) then; ok = .true.; exit; end if
      pp = iapws95_rhoT(rho, T); fr = pp%P - Ptar; dfr = pp%dpdrho_T
      if (fr < 0._r8) then
        a = rho
      else
        b = rho
      end if
    end do
  end function solve_rho_bracket

  !=====================================================================
  ! Two-phase saturation at temperature T by the Maxwell construction.
  ! Solves, for the reduced vapour/liquid densities (dv,dl):
  !   F1 = dv(1+dv*firdv) - dl(1+dl*firdl)                = 0   (equal P)
  !   F2 = dv*firdv - dl*firdl + ln(dv/dl) + (firv-firl)  = 0   (equal Gibbs)
  ! with the Supp-sat 1992 auxiliary densities as the initial guess.
  ! Returns Psat and the saturated vapour/liquid property structs.
  !=====================================================================
  subroutine iapws95_sat(T, Psat, vap, liq, ok)
    real(r8),              intent(in)  :: T
    real(r8),              intent(out) :: Psat
    type(iapws95_props_t), intent(out) :: vap, liq
    logical,               intent(out) :: ok
    real(r8) :: tau, dv, dl
    real(r8) :: firv, firdv, firddv, firtv, firttv, firdtv
    real(r8) :: firl, firdl, firddl, firtl, firttl, firdtl
    real(r8) :: dum1, dum2, dum3
    real(r8) :: F1, F2, j11, j12, j21, j22, det, ddv, ddl
    integer  :: it
    integer, parameter :: maxit = 100
    real(r8), parameter :: ftol = 1.e-11_r8

    ok = .false.
    if (T >= IAPWS_TC .or. T < 50._r8) then
      Psat = iapws95_psat_aux(T)
      vap  = iapws95_rhoT(IAPWS_RHOC, T)
      liq  = vap
      return
    end if

    tau = IAPWS_TC / T
    dv  = rho_vap_aux(T) / IAPWS_RHOC
    dl  = rho_liq_aux(T) / IAPWS_RHOC

    do it = 1, maxit
      call phir_derivs(tau, dv, firv, firdv, firddv, firtv, firttv, firdtv)
      call phir_derivs(tau, dl, firl, firdl, firddl, firtl, firttl, firdtl)

      F1 = dv*(1._r8 + dv*firdv) - dl*(1._r8 + dl*firdl)
      F2 = dv*firdv - dl*firdl + log(dv/dl) + (firv - firl)

      ! Jacobian
      j11 =  1._r8 + 2._r8*dv*firdv + dv*dv*firddv     ! dF1/ddv = den(dv)
      j12 = -(1._r8 + 2._r8*dl*firdl + dl*dl*firddl)   ! dF1/ddl
      j21 =  2._r8*firdv + dv*firddv + 1._r8/dv        ! dF2/ddv
      j22 = -(2._r8*firdl + dl*firddl + 1._r8/dl)      ! dF2/ddl

      det = j11*j22 - j12*j21
      if (abs(det) < 1.e-30_r8) exit
      ddv = ( F1*j22 - F2*j12) / det
      ddl = (-F1*j21 + F2*j11) / det
      ! bounded, sign-preserving update
      if (ddv >  0.4_r8*dv) ddv =  0.4_r8*dv
      if (ddv < -0.4_r8*dv) ddv = -0.4_r8*dv
      if (ddl >  0.4_r8*dl) ddl =  0.4_r8*dl
      if (ddl < -0.4_r8*dl) ddl = -0.4_r8*dl
      dv = dv - ddv
      dl = dl - ddl
      if (max(abs(F1),abs(F2)) < ftol) then
        ok = .true.
        exit
      end if
    end do

    vap  = iapws95_rhoT(dv*IAPWS_RHOC, T)
    liq  = iapws95_rhoT(dl*IAPWS_RHOC, T)
    Psat = 0.5_r8*(vap%P + liq%P)
    ! keep unused derivative temporaries referenced for clarity
    dum1 = firttv; dum2 = firttl; dum3 = firdtv + firdtl
    if (dum1 + dum2 + dum3 /= dum1 + dum2 + dum3) ok = .false.   ! NaN guard
  end subroutine iapws95_sat

end module exocol_iapws95
