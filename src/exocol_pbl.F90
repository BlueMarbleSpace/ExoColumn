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
  public :: pbl_reset

  ! ---- Closure constants (Frierson et al. 2006, Table 1) ----
  real(r8), parameter :: vk    = 0.4_r8       ! von Karman constant
  real(r8), parameter :: ri_c  = 1.0_r8       ! critical bulk Richardson number
  real(r8), parameter :: fb    = 0.1_r8       ! surface-layer fraction of BL depth
  real(r8), parameter :: eps_v = 0.608_r8     ! virtual-temperature coefficient
  real(r8), parameter :: u_min = 1.0_r8       ! floor on wind speed in Ri [m/s]

  ! ---- Prognostic, fixed-depth-anchored BL depth (resolution-independence fix) ----
  ! The bulk-Richardson BL depth has no stable fixed point at the physical
  ! (~500 m) depth in single-column RCE.  Re-diagnosed fresh every step it is
  ! caught between two runaways with no self-limiting mechanism (a GCM's resolved
  ! variability + horizontal averaging supply that; a single column has neither):
  !   * DEEP runaway (raw bulk-Ri / no cap): once the sub-cloud layer is mixed the
  !     diagnosed top sits just above it and creeps to the tropopause; whole-column
  !     mixing brings Θ_BL→Ts and the sensible flux collapses (SH→0).
  !   * COLLAPSE (LCL / surface-parcel cap): a shallow BL over-moistens the bottom
  !     layer → q_bot saturates → the LCL drops → shallower still → LE→0.
  ! Relaxing a prognostic depth toward either target does not cure it — the
  ! *target* itself runs away, so the prognostic depth just tracks it.
  !
  ! Fix (two parts):
  !  1. ANCHOR the nominal mixing depth to a FIXED pressure depth dp_mix (~500 m,
  !     a typical marine boundary layer), and let the bulk-Richardson height only
  !     REDUCE it when the surface layer is statically stable:
  !         target = max( min(bulk-Ri height, h_fixed), h_floor ).
  !     h_fixed caps the deep runaway; h_floor prevents collapse onto the
  !     radiatively-stiff bottom grid layer.  All three are fixed pressure
  !     intervals → the same physical height at every vertical resolution, which
  !     is precisely what makes the near-surface state (and hence LE/SH and the
  !     climate) resolution-independent.
  !  2. Carry the depth as PROGNOSTIC state and RELAX it toward that target over
  !     tau_h rather than snapping to the jumpy instantaneous value (the lag
  !     smooths step-to-step jitter in the bulk-Ri diagnosis).  Deepening is
  !     additionally rate-limited by an entrainment velocity, and h is capped at a
  !     pressure-depth guardrail.  (Frierson et al. call the GCM depth "prognostic"
  !     = diagnosed each step; here it is genuinely prognostic — integrated state.)
  real(r8), parameter :: tau_h   = 3._r8 * 3600._r8  ! BL-depth relaxation timescale [s]
  real(r8), parameter :: we_max  = 0.1_r8            ! max entrainment (deepening) velocity [m/s]
  real(r8), parameter :: dp_mix   = 6.0e3_r8         ! nominal BL pressure depth below surface [Pa] (~500 m)
  real(r8), parameter :: dp_floor = 2.0e3_r8         ! floor: min BL pressure depth below surface [Pa] (~170 m)
  real(r8), parameter :: dp_cap   = 3.0e4_r8         ! guardrail: max BL pressure depth below surface [Pa]

  ! Prognostic boundary-layer depth [m].  < 0 = uninitialised (set to the first
  ! diagnosed target).  Reset by pbl_reset() when the RCE loop (re-)starts so a
  ! fresh run does not inherit the previous run's boundary-layer state.
  real(r8), save :: h_bl_prog = -1._r8

contains

  subroutine pbl_reset()
  ! Clear the prognostic BL depth so the next RCE run initialises it afresh.
    h_bl_prog = -1._r8
  end subroutine pbl_reset

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
    real(r8) :: h_bulk, h_fixed, h_target, h_floor, h_cap
    real(r8) :: rib_below, frac, dh
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

    ! --- bulk-Richardson height (base = lowest level), interpolated ---
    ! Ri = 0 at the base and increases upward; h_bulk = height where Ri first
    ! reaches Ri_c.  This is the textbook BL top, but in single-column RCE it has
    ! no stable fixed point (it drifts up into the weakly dry-stable moist-adiabat
    ! region), so it is bounded below by the LCL (next block).
    h_bulk    = zint(1)            ! no crossing → whole column (capped below)
    rib_below = 0._r8             ! Ri at the base (z = z_base)
    do k = pver-1, 1, -1
      rib = exo_g * zmid(k) * (sv(k) - sv_base) / (sv_base * U*U)
      if (rib >= ri_c) then
        if (rib > rib_below) then
          frac = (ri_c - rib_below) / (rib - rib_below)
        else
          frac = 0._r8
        end if
        h_bulk = zmid(k+1) + frac * (zmid(k) - zmid(k+1))
        exit
      end if
      rib_below = rib
    end do

    ! --- nominal mixing depth + floor / guardrail (fixed pressure depths) ---
    ! Anchor the BL to a fixed physical depth dp_mix (~500 m); bulk-Ri only makes
    ! it shallower under stable conditions.  Floor and guardrail bound it against
    ! the collapse / deep-runaway pathologies.  All are fixed pressure intervals,
    ! hence identical physical heights at every vertical resolution.
    h_fixed = bl_height_at(pint(pver+1) - dp_mix,   pint, zint, pver)
    h_floor = bl_height_at(pint(pver+1) - dp_floor, pint, zint, pver)
    h_cap   = bl_height_at(pint(pver+1) - dp_cap,   pint, zint, pver)

    h_target = min(h_bulk, h_fixed)      ! bulk-Ri can only make it shallower
    h_target = max(h_target, h_floor)    ! never collapse onto the bottom layer

    ! --- relax the prognostic depth toward the target (the stability fix) ---
    if (h_bl_prog < 0._r8) then
      h_bl_prog = h_target                              ! initialise on first call
    else
      dh = (h_target - h_bl_prog) * min(dt_sec / tau_h, 1._r8)
      if (dh > we_max * dt_sec) dh = we_max * dt_sec    ! entrainment-limited deepening
      h_bl_prog = h_bl_prog + dh
    end if
    h_bl_prog = max(h_bl_prog, h_floor)                 ! pressure-depth floor
    h_bl_prog = min(h_bl_prog, h_cap)                   ! pressure-depth guardrail
    hbl       = h_bl_prog
    h_diag    = hbl

    ! ktop = highest (smallest-k) layer whose midpoint lies within the BL.
    ktop = pver
    do k = pver-1, 1, -1
      if (zmid(k) <= hbl) then
        ktop = k
      else
        exit
      end if
    end do
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

  pure function bl_height_at(p_target, pint, zint, nv) result(z_out)
  ! Height [m] of the p_target pressure level, linearly interpolated in pressure
  ! from the interface (pint, zint) arrays.  pint increases downward (index 1 =
  ! TOA, small p); zint decreases downward (zint(nv+1) = 0 at the surface).
  ! Clamps to the column top when p_target lies above p_top.  Used to convert the
  ! fixed BL-depth floor/guardrail pressure depths into heights (so they are the
  ! same physical depth at every vertical resolution).
    integer,  intent(in) :: nv
    real(r8), intent(in) :: p_target
    real(r8), intent(in) :: pint(nv+1), zint(nv+1)
    real(r8) :: z_out, f
    integer  :: k
    z_out = zint(1)               ! default: target at/above column top
    do k = nv, 1, -1
      if (pint(k) <= p_target) then
        f = (pint(k+1) - p_target) / (pint(k+1) - pint(k))
        z_out = zint(k+1) + f * (zint(k) - zint(k+1))
        return
      end if
    end do
  end function bl_height_at

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
