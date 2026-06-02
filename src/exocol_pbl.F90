module exocol_pbl
! Boundary-layer vertical mixing for ExoColumn — Frierson, Held & Zurita-Gotor
! (2006) simplified diffusive boundary layer (the lineage of ExoColumn's SBM
! convection).
!
! WHY THIS EXISTS
! ----------------
! Without boundary-layer mixing the near-surface T and q — which set the bulk
! surface fluxes — are determined by the arbitrary thickness of the bottom grid
! layer, so the equilibrium climate is RESOLUTION-DEPENDENT (a ~4x swing in SH
! and a ~9 K swing in Ts between pver=70 and 140; see
! project_resolution_root_cause).  This scheme mixes dry static energy s = cp*T
! + g*z and specific humidity q through the diagnosed boundary layer with the
! Monin-Obukhov-consistent eddy diffusivity of Frierson et al. (2006).  Together
! with the MOS surface fluxes (exocol_surface, 'mos'), the near-surface state
! becomes resolution-independent.
!
! FORMULATION (Frierson et al. 2006, eqs 16-20)
! ---------------------------------------------
! BL depth h = lowest height where the bulk Richardson number relative to the
! lowest model level reaches Ri_c = 1:
!     Ri(z) = g·z·(s_v(z) − s_v(z_a)) / (s_v(z_a)·U²),   s_v = cp·T_v + g·z.
! Eddy diffusivity (mechanical, wind-driven; matched to the MOS surface layer):
!     z < fb·h :  K = κ·U·√C·z
!     fb·h<z<h :  K = κ·U·√C·(fb·h)·(z/(fb·h))·[1 − (z−fb·h)/((1−fb)h)]²
! with κ=0.4, fb=0.1, C the MOS drag coefficient from exocol_surface, U the
! (prescribed) surface wind.  Momentum is not diffused (no resolved wind).
!
! NUMERICS
! --------
! The diffusion of s and q is solved with an implicit (backward-Euler)
! conservative tridiagonal sweep (no-flux top/bottom): it conserves column dry
! static energy Sum(s·dp/g) and water Sum(q·dp/g) to machine precision and is
! unconditionally stable for any K and dt (a thin, vigorously-mixed bottom layer
! never constrains the timestep).

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: SHR_CONST_RGAS
  use exoplanet_mod,  only: exo_g

  implicit none
  private

  public :: pbl_diffuse

  ! ---- Closure constants (Frierson et al. 2006, Table 1) ----
  real(r8), parameter :: vk    = 0.4_r8       ! von Karman constant
  real(r8), parameter :: ri_c  = 1.0_r8       ! critical bulk Richardson number
  real(r8), parameter :: fb    = 0.1_r8       ! surface-layer fraction of BL depth
  real(r8), parameter :: eps_v = 0.608_r8     ! virtual-temperature coefficient
  real(r8), parameter :: u_min = 1.0_r8       ! floor on wind speed in Ri [m/s]

contains

  subroutine pbl_diffuse(tmid, h2ommr, pmid, pint, pdel, zint, &
                         mwdry, cpdry, wind, C_drag, surf_sh, surf_e, dt_sec, &
                         mix_q, pver, h_diag, kmax_diag)
  ! Mix dry static energy and (optionally) q through the boundary layer.
  ! tmid and h2ommr are updated in place.  Conserves column DSE and water.
  !
  ! Index convention: k=1 TOA, k=pver surface layer; zint(pver+1)=0 at surface.
    integer,  intent(in)    :: pver
    real(r8), intent(inout) :: tmid(pver)     ! temperature [K]
    real(r8), intent(inout) :: h2ommr(pver)   ! specific humidity [kg/kg]
    real(r8), intent(in)    :: pmid(pver)     ! midpoint pressure [Pa]
    real(r8), intent(in)    :: pint(pver+1)   ! interface pressure [Pa]
    real(r8), intent(in)    :: pdel(pver)     ! layer pressure thickness [Pa]
    real(r8), intent(in)    :: zint(pver+1)   ! interface height [m] (zint(pver+1)=0)
    real(r8), intent(in)    :: mwdry          ! dry-air molar mass [g/mol]
    real(r8), intent(in)    :: cpdry          ! dry-air specific heat [J/kg/K]
    real(r8), intent(in)    :: wind           ! surface wind speed [m/s]
    real(r8), intent(in)    :: C_drag         ! MOS drag coefficient [-]
    real(r8), intent(in)    :: surf_sh        ! surface sensible heat flux [W/m2] (BL bottom source)
    real(r8), intent(in)    :: surf_e         ! surface evaporation [kg/m2/s]  (BL bottom source)
    real(r8), intent(in)    :: dt_sec         ! timestep [s]
    logical,  intent(in)    :: mix_q          ! also mix moisture
    real(r8), intent(out)   :: h_diag         ! diagnosed BL height [m]
    integer,  intent(out)   :: kmax_diag      ! highest (smallest-k) mixed layer

    real(r8) :: Rd, U, sqrtC
    real(r8) :: zmid(pver), sv(pver)
    real(r8) :: gco(pver), Mlay(pver), s(pver), qn(pver)
    real(r8) :: z_base, sv_base, hbl, hsl, Kb_hsl
    real(r8) :: rib, z_if, rho_if, t_if, dz_mid, Kc, zz
    integer  :: k, ktop

    Rd = SHR_CONST_RGAS / mwdry
    U  = max(wind, u_min)

    h_diag    = 0._r8
    kmax_diag = pver

    ! No exchange if the surface layer carries no drag (e.g. statically stable
    ! cutoff, C = 0) or there is no wind.
    if (C_drag <= 0._r8 .or. wind <= 0._r8) return
    sqrtC = sqrt(C_drag)

    ! --- midpoint heights and virtual dry static energy ---
    do k = 1, pver
      zmid(k) = 0.5_r8 * (zint(k) + zint(k+1))
      sv(k)   = cpdry * tmid(k) * (1._r8 + eps_v * h2ommr(k)) + exo_g * zmid(k)
    end do
    z_base  = zmid(pver)
    sv_base = sv(pver)

    ! --- BL height: bulk Richardson (base = lowest level) reaches Ri_c ---
    ktop = pver
    do k = pver-1, 1, -1
      rib = exo_g * zmid(k) * (sv(k) - sv_base) / (sv_base * U*U)
      if (rib >= ri_c) exit
      ktop = k
    end do

    ! BL top from the bulk-Richardson criterion alone.  (An earlier LCL cap, added
    ! to tame a depth runaway, is unnecessary once the local physics is sub-stepped
    ! — SBM then restratifies between mixings so the bulk-Ri depth is stable — and
    ! the cap itself caused a q_bot-saturation collapse.  Kept out deliberately.)
    hbl       = zmid(ktop)
    h_diag    = hbl
    kmax_diag = ktop
    ! NB: no early return when ktop == pver — the implicit solve below still
    ! applies the surface flux source (it degenerates to a bottom-layer deposit
    ! when there are no interior conductances to spread it).

    hsl    = fb * hbl                          ! surface-layer top
    Kb_hsl = vk * U * sqrtC * hsl              ! K at the surface-layer top

    ! --- interface eddy diffusivity and conductance ---
    ! gco(k) couples layer k (above) and k+1 (below) across interface zint(k+1).
    gco = 0._r8
    do k = 1, pver-1
      z_if = zint(k+1)
      if (z_if >= hbl) then
        Kc = 0._r8
      else if (z_if <= hsl) then
        Kc = vk * U * sqrtC * z_if                                  ! surface layer
      else
        zz = (z_if - hsl) / ((1._r8 - fb) * hbl)
        Kc = Kb_hsl * (z_if / hsl) * (1._r8 - zz)**2                ! outer (parabolic)
      end if
      if (Kc > 0._r8) then
        t_if   = 0.5_r8 * (tmid(k) + tmid(k+1))
        rho_if = pint(k+1) / (Rd * t_if)
        dz_mid = zmid(k) - zmid(k+1)                                ! > 0
        gco(k) = rho_if * Kc / dz_mid
      end if
    end do

    ! --- layer mass and conserved scalars (dry static energy, q) ---
    do k = 1, pver
      Mlay(k) = pdel(k) / exo_g
      s(k)    = cpdry * tmid(k) + exo_g * zmid(k)
      qn(k)   = h2ommr(k)
    end do

    ! --- implicit tridiagonal diffusion with surface flux as bottom source ---
    ! The surface sensible flux (W/m²) and evaporation (kg/m²/s) enter as a
    ! source into the bottom layer and are spread through the BL within the same
    ! implicit solve, so the thin bottom layer is never spiked.
    call tridiag_diffuse(Mlay, gco, dt_sec, surf_sh, s, pver)
    if (mix_q) call tridiag_diffuse(Mlay, gco, dt_sec, surf_e, qn, pver)

    ! --- recover temperature; commit ---
    do k = 1, pver
      tmid(k) = (s(k) - exo_g * zmid(k)) / cpdry
      if (mix_q) h2ommr(k) = max(qn(k), 0._r8)
    end do

  end subroutine pbl_diffuse

  ! -----------------------------------------------------------------------

  subroutine tridiag_diffuse(M, gco, dt, src_bot, phi, pver)
  ! Backward-Euler conservative vertical diffusion of phi in place, with a
  ! surface source src_bot ([phi]·kg/m²/s) injected into the bottom layer:
  !   M_k (phi_k^{n+1}-phi_k^n)/dt = gco_{k-1}(phi_{k-1}-phi_k)
  !                                 - gco_k(phi_k-phi_{k+1}) + δ_{k,pver}·src_bot
  ! with gco(0)=gco(pver)=0 (no-flux top, surface-source bottom).  Conserves
  ! Sum(M_k phi_k) up to the injected dt·src_bot.
    integer,  intent(in)    :: pver
    real(r8), intent(in)    :: M(pver)        ! layer mass [kg/m2]
    real(r8), intent(in)    :: gco(pver)      ! interface conductance; gco(k)=k|k+1, gco(pver)=0
    real(r8), intent(in)    :: dt
    real(r8), intent(in)    :: src_bot        ! surface source into layer pver [phi·kg/m2/s]
    real(r8), intent(inout) :: phi(pver)

    real(r8) :: a(pver), b(pver), c(pver), d(pver)
    real(r8) :: gl, gu, w
    integer  :: k

    do k = 1, pver
      gu = 0._r8;  if (k > 1)    gu = gco(k-1)
      gl = 0._r8;  if (k < pver) gl = gco(k)
      a(k) = -dt * gu
      c(k) = -dt * gl
      b(k) =  M(k) + dt * (gu + gl)
      d(k) =  M(k) * phi(k)
    end do
    d(pver) = d(pver) + dt * src_bot          ! surface flux into the bottom layer

    c(1) = c(1) / b(1)
    d(1) = d(1) / b(1)
    do k = 2, pver
      w    = b(k) - a(k) * c(k-1)
      c(k) = c(k) / w
      d(k) = (d(k) - a(k) * d(k-1)) / w
    end do
    phi(pver) = d(pver)
    do k = pver-1, 1, -1
      phi(k) = d(k) - c(k) * phi(k+1)
    end do

  end subroutine tridiag_diffuse

end module exocol_pbl
