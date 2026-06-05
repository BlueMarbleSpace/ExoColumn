module exocol_coldstart
! Self-contained cold-start initialization for ExoColumn.
!
! Called when &exocol_init::input_file is empty.  Builds a complete initial
! column state from the namelist:
!
!   1. Dry-air composition: built-in Earth-like defaults (N2 = 0.78,
!      O2 = 0.21, CO2 = 4e-4, others = 0), overridden by any &exocol_composition
!      VMR settings.  mwdry recomputed from the resulting composition.
!   2. Pressure grid: log-spaced interfaces from p_top to ps over pverp levels.
!      ps comes from &exocol_composition::ps if set, else defaults to 1 bar.
!   3. Temperature: moist adiabat integrated upward from ts using the local
!      Γm(T,p) from exocol_convadj::malr.  Capped at t_strato above the
!      tropopause (isothermal stratosphere).  Numerical scheme is a sub-stepped
!      Euler in log(p):  dT/d(ln p) = Γm · R_d · T_v / g,  where T_v is the
!      virtual temperature T_v = T·(1 + r_s/ε)/(1 + r_s), r_s = saturation
!      mixing ratio at (T,p), ε = R_d/R_v.  Using T_v (rather than T) gives the
!      correct moist scale height when H2O is non-negligible (inner HZ / runaway).
!   4. Moisture: Kasting (1988, Icarus) cold-trap convention.
!      Troposphere (T > t_strato): h2ommr = rh_init · qsat(T, p).
!      Stratosphere (T ≤ t_strato): h2ommr = rh_init · qsat(t_strato, p_tropo),
!      uniform throughout the stratosphere at the value air had when it last
!      passed through the cold point.  p_tropo is the upper interface of the
!      topmost tropospheric layer (= pint(k_top_conv), the tropopause pressure).
!      This is physically correct: above the cold point there is no mechanism
!      to add water, so the mixing ratio is frozen at the tropopause value.
!   5. Heights: hydrostatic integration  zint(k) = zint(k+1) + R_d·T_v/g · ln(p(k+1)/p(k))
!      using the virtual temperature T_v = T·(1 + h2ommr/ε)/(1 + h2ommr).
!   6. Surface scalars (ts, coszrs, albedos), msdist (from &exocol_nml),
!      and clouds (zeroed) are written from the namelist values.
!
! After this routine returns, the column is ready for exocol_setgas() and
! the RCE loop.  exocol_update_derived() is called internally so pdeldry /
! pintdry are consistent on exit.

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_const_mod,  only: SHR_CONST_RGAS, SHR_CONST_MWWV
  use exoplanet_mod,  only: exo_g
  use ppgrid,         only: pver, pverp
  use exocol_mod
  use exocol_config,  only: cs_ts       => ts,      &
                            cs_coszrs   => coszrs,  &
                            cs_asdir    => asdir,   &
                            cs_asdif    => asdif,   &
                            cs_aldir    => aldir,   &
                            cs_aldif    => aldif,   &
                            cs_ps       => ps,      &  ! &exocol_composition
                            cfg_msdist  => msdist,  &  ! &exocol_nml
                            t_strato, p_top, rh_init,             &
                            co2_vmr, ch4_vmr, c2h6_vmr,           &
                            h2_vmr,  n2_vmr,  o3_vmr, o2_vmr,     &
                            ar_vmr,  o3_profile, variable_ps, ihz_profile, &
                            h2o_inventory_bar, h2o_eos,            &
                            MW_CO2, MW_CH4, MW_C2H6, MW_H2,       &
                            MW_N2,  MW_O3,  MW_O2,  MW_AR,         &
                            CP_CO2, CP_CH4, CP_C2H6, CP_H2,       &
                            CP_N2,  CP_O3,  CP_O2,  CP_AR
  use exocol_convadj, only: esat, malr
  use exocol_steam,   only: steam_dlnTdlnP_sat, steam_dlnTdlnP_dry
  use exocol_iapws95, only: IAPWS_TT   ! water triple-point T (L-V saturation floor)
  use exocol_ozone,   only: set_earth_o3_profile, set_rcemip_o3_profile

  implicit none
  private

  public :: cold_start_init

  ! Earth-like dry-air defaults used when a gas VMR is not specified in
  ! &exocol_composition.  N2/O2/Ar/CO2 are the standard dry-air mole fractions and
  ! sum to ~1.0.  Argon (DEF_AR) is radiatively inert so it carries no ExoRT gas
  ! tracer, but it IS included in mwdry and cpdry (mwdry → ~28.97 g/mol, cpdry →
  ! ~1004 J/kg/K), matching real dry air and konrad's constants.
  real(r8), parameter :: DEF_CO2  = 4.0e-4_r8
  real(r8), parameter :: DEF_CH4  = 0.0_r8
  real(r8), parameter :: DEF_C2H6 = 0.0_r8
  real(r8), parameter :: DEF_H2   = 0.0_r8
  real(r8), parameter :: DEF_N2   = 0.78084_r8
  real(r8), parameter :: DEF_O3   = 0.0_r8
  real(r8), parameter :: DEF_O2   = 0.20946_r8
  real(r8), parameter :: DEF_AR   = 0.00934_r8   ! argon (inert; mw/cp only)

  ! Cold-start ps default if &exocol_composition::ps is unset.
  real(r8), parameter :: DEFAULT_PS = 1.0e5_r8   ! 1 bar

  ! Substeps per layer for moist-adiabat integration.
  integer,  parameter :: NSUB = 20

contains

  ! -----------------------------------------------------------------------

  subroutine cold_start_init()
    integer  :: k, isub
    real(r8) :: ps_use, vmr_sum, mwdry_new, cpdry_new, Rd, eps_wv
    real(r8) :: vmr_dry(7), Mw_dry(7), mmr_dry(7), cp_dry(7)
    real(r8) :: ar_vmr_use, ar_mmr
    real(r8) :: T_at_int(pverp)
    real(r8) :: T_lev, p_lev, dlogp_layer, dlogp_sub, dT_sub, Gm
    real(r8) :: es_sub, r_sub, T_v_lev
    real(r8) :: es_k, qsat_k, e_sfc
    integer  :: k_top_conv
    real(r8) :: es_tropo, q_cold_trap, p_coldpoint, p_tropo_use
    real(r8) :: T_prev_sub, p_prev_sub
    character(len=4) :: dry_name(7)
    ! Kasting IHZ dry-adiabat variables (initialised to safe no-dry-layer defaults)
    real(r8) :: f_vmr, Mw_mix, w_h2o, Rv, R_mix, cp_h2o_steam, cp_mix, kappa_mix
    real(r8) :: e_sat_sfc, p_h2o_inv, f_vmr_surf
    real(r8) :: h2ommr_dry = 0._r8  ! constant H2O MMR in dry layer
    integer  :: k_cond    = 0       ! condensation level index (0 = no dry layer)
    logical  :: in_dry_layer
    logical  :: surface_saturated   ! .true. if esat(Ts) <= water inventory (ocean present)
    ! Non-ideal (Kasting 1988 Appendix-A) moist pseudoadiabat selector + the
    ! fixed mass mixing ratio alpha_v = rho_v/rho_n carried through the dry layer.
    logical  :: use_nonideal_adiabat
    real(r8) :: alfav_dry = 0._r8

    write(*,'(/,a)') '  cold_start_init: building initial column from &exocol_init'

    ! ---- 1. Composition (defaults + &exocol_composition overrides) ----
    vmr_dry  = (/ DEF_CO2, DEF_CH4, DEF_C2H6, DEF_H2, DEF_N2, DEF_O3, DEF_O2 /)
    Mw_dry   = (/ MW_CO2,  MW_CH4,  MW_C2H6,  MW_H2,  MW_N2,  MW_O3,  MW_O2  /)
    cp_dry   = (/ CP_CO2,  CP_CH4,  CP_C2H6,  CP_H2,  CP_N2,  CP_O3,  CP_O2  /)
    dry_name = (/ 'CO2 ',  'CH4 ',  'C2H6',   'H2  ', 'N2  ', 'O3  ', 'O2  ' /)

    if (co2_vmr  >= 0._r8) vmr_dry(1) = co2_vmr
    if (ch4_vmr  >= 0._r8) vmr_dry(2) = ch4_vmr
    if (c2h6_vmr >= 0._r8) vmr_dry(3) = c2h6_vmr
    if (h2_vmr   >= 0._r8) vmr_dry(4) = h2_vmr
    if (n2_vmr   >= 0._r8) vmr_dry(5) = n2_vmr
    if (o3_vmr   >= 0._r8) vmr_dry(6) = o3_vmr
    if (o2_vmr   >= 0._r8) vmr_dry(7) = o2_vmr

    ! Argon (inert): Earth default unless overridden.  Not carried as a tracer;
    ! enters mwdry and cpdry only.
    ar_vmr_use = DEF_AR
    if (ar_vmr >= 0._r8) ar_vmr_use = ar_vmr

    vmr_sum = sum(vmr_dry) + ar_vmr_use
    if (abs(vmr_sum - 1.0_r8) > 1.0e-2_r8) then
      write(*,'(a,f8.5,a)') &
        '  WARNING: cold-start dry-air VMRs sum to ', vmr_sum, &
        ' (expected ~1.0). Values used as-is.'
    end if

    mwdry_new = sum(vmr_dry * Mw_dry) + ar_vmr_use * MW_AR
    if (mwdry_new <= 0._r8) then
      write(*,'(a)') '  cold_start_init: computed mwdry <= 0. Aborting.'
      stop 'cold_start_init: bad composition'
    end if

    do k = 1, 7
      mmr_dry(k) = vmr_dry(k) * Mw_dry(k) / mwdry_new
    end do
    ar_mmr = ar_vmr_use * MW_AR / mwdry_new

    ! Mass-weighted mean cp from the actual composition (argon included).  This
    ! replaces the namelist cpdry so non-Earth compositions (especially H2-rich)
    ! get the correct specific heat automatically.
    cpdry_new = sum(mmr_dry * cp_dry) + ar_mmr * CP_AR

    write(*,'(a,f7.3,a)') '    mwdry  = ', mwdry_new, ' g/mol'
    write(*,'(a,f7.1,a)') '    cpdry  = ', cpdry_new, ' J/kg/K (auto-computed from composition)'
    if (ar_vmr_use > 0._r8) write(*,'(a,es10.3,a,es10.3)') &
      '    Ar    VMR = ', ar_vmr_use, '  (inert)  MMR = ', ar_mmr
    do k = 1, 7
      if (vmr_dry(k) > 0._r8) then
        write(*,'(a,a4,a,es10.3,a,es10.3)') &
          '    ', dry_name(k), '  VMR = ', vmr_dry(k), '  MMR = ', mmr_dry(k)
      end if
    end do

    ! ---- 2. Pressure grid ----
    if (cs_ps >= 0._r8) then
      ps_use = cs_ps
    else
      ps_use = DEFAULT_PS
    end if

    ! With variable_ps, ps_use is the dry background; total ps = p_dry + pH2O,
    ! where the surface water partial pressure is CAPPED at the finite ocean
    ! inventory (Kasting 1988 Eq. 2): pH2O = min(esat(Ts), inventory).  The
    ! surface is saturated only while esat(Ts) <= inventory (ocean present);
    ! above that the ocean is fully evaporated and the surface is subsaturated.
    e_sat_sfc        = esat(cs_ts)
    p_h2o_inv        = h2o_inventory_bar * 1.0e5_r8   ! inventory pressure [Pa]
    surface_saturated = .true.
    if (variable_ps) then
      e_sfc             = min(e_sat_sfc, p_h2o_inv)
      surface_saturated = (e_sat_sfc <= p_h2o_inv)
      write(*,'(a,f9.2,a,f9.2,a,f9.2,a)') &
        '    variable_ps: p_dry = ', ps_use/100._r8, ' hPa  + pH2O = ', &
        e_sfc/100._r8, ' hPa  →  ps = ', (ps_use + e_sfc)/100._r8, ' hPa'
      if (.not. surface_saturated) &
        write(*,'(a,f9.2,a,f9.2,a)') &
          '      surface SUBSATURATED: esat(Ts) = ', e_sat_sfc/1.0e5_r8, &
          ' bar > inventory ', h2o_inventory_bar, ' bar (ocean evaporated → dry base)'
      ps_use = ps_use + e_sfc
    else
      ! Fixed-ps mode: surface water cannot exceed the column; saturated if it fits.
      e_sfc             = min(e_sat_sfc, 0.99_r8 * ps_use)
      surface_saturated = (e_sat_sfc <= 0.99_r8 * ps_use)
    end if
    f_vmr_surf = e_sfc / ps_use                       ! surface H2O mole fraction

    if (p_top >= ps_use) then
      write(*,'(a)') '  cold_start_init: p_top >= ps. Aborting.'
      stop 'cold_start_init: bad pressure range'
    end if

    call build_pressure_grid(ps_use)

    write(*,'(a,es10.3,a,es10.3,a)') &
      '    p_top → ps : ', p_top, ' Pa → ', ps_use, ' Pa'
    write(*,'(a,f8.2,a,f8.2,a)') &
      '    pmid(pver) : ', pmid(pver)/100._r8, ' hPa  (bottom midpoint ≈', &
      (ps_use - pmid(pver)) / (1.2_r8 * exo_g), ' m altitude)'

    ! ---- 3. Temperature: Kasting IHZ (dry + moist) or standard moist adiabat ----
    Rd     = SHR_CONST_RGAS / mwdry_new
    Rv     = SHR_CONST_RGAS / SHR_CONST_MWWV   ! H2O gas constant [J/kg/K]
    eps_wv = SHR_CONST_MWWV / mwdry_new

    ! Kasting (1988) Appendix-A general non-ideal moist pseudoadiabat (IAPWS-95
    ! EOS) when h2o_eos='nonideal'; otherwise the dilute ideal-gas malr.
    use_nonideal_adiabat = (trim(h2o_eos) == 'nonideal')

    T_at_int(pverp) = cs_ts
    p_coldpoint     = -1._r8   ! true cold-point pressure; interpolated below

    if (cs_ts <= t_strato) then
      write(*,'(a)') &
        '  cold_start_init: ts <= t_strato — isothermal column above surface.'
      T_at_int(1:pverp-1) = t_strato

    else if (ihz_profile) then
      ! Kasting (1988) IHZ structure (p.476, Figs 1,4,6), inventory-based switch:
      !   surface SATURATED  (esat(Ts) <= inventory; ocean present) → moist
      !       pseudoadiabat from the surface (moist-greenhouse branch, Fig 1a);
      !   surface SUBSATURATED (ocean fully evaporated) → DRY adiabat at the base
      !       with H2O mixing ratio fixed by the inventory, joining a moist
      !       pseudoadiabat aloft where the rising parcel first reaches saturation
      !       (runaway branch, Fig 1b).
      ! Dry-adiabat lapse rate: dT/dp = κ_mix · T/p (exact for ideal gas),
      ! κ_mix = R_mix/cp_mix for the H2O + background mixture at the surface VMR.
      f_vmr = f_vmr_surf      ! constant H2O mole fraction in the dry region

      if (.not. surface_saturated) then
        ! Compute dry-adiabat kappa for the H2O + dry-gas mixture.
        Mw_mix      = f_vmr * SHR_CONST_MWWV + (1._r8 - f_vmr) * mwdry_new
        w_h2o       = f_vmr * SHR_CONST_MWWV / Mw_mix      ! H2O mass fraction
        R_mix       = w_h2o * Rv + (1._r8 - w_h2o) * Rd
        cp_h2o_steam = 1870._r8                              ! J/kg/K, steam ~300-700 K
        cp_mix      = w_h2o * cp_h2o_steam + (1._r8 - w_h2o) * cpdry_new
        kappa_mix   = R_mix / cp_mix
        ! Constant H2O MMR in the dry layer
        h2ommr_dry  = w_h2o
        ! Fixed mass ratio alpha_v = rho_v/rho_n carried through the dry layer
        ! (Kasting A3), used by the non-ideal A11/A12 unsaturated lapse rate.
        alfav_dry   = w_h2o / (1._r8 - w_h2o)

        write(*,'(a,f7.4,a,f6.4,a)') &
          '    IHZ runaway: subsaturated surface, dry base  f_vmr = ', f_vmr, &
          '  κ_mix = ', kappa_mix, ' (dry adiabat from surface)'

        ! Fill interface temperatures: dry adiabat until ws ≤ 1, then moist.
        k_cond      = 1          ! fallback: entire column is dry → moist at TOA
        in_dry_layer = .true.

        do k = pverp-1, 1, -1
          if (.not. in_dry_layer) then
            ! Moist adiabat from the condensation level upward
            if (T_at_int(k+1) <= t_strato) then
              T_at_int(k) = t_strato
              cycle
            end if
            T_lev       = T_at_int(k+1)
            p_lev       = pint(k+1)
            dlogp_layer = log(pint(k) / pint(k+1))
            dlogp_sub   = dlogp_layer / real(NSUB, r8)
            do isub = 1, NSUB
              if (use_nonideal_adiabat .and. T_lev >= IAPWS_TT) then
                ! Kasting A4/A5 non-ideal saturated moist pseudoadiabat (liquid-
                ! vapour SVP; below the triple point fall back to the ideal malr)
                dT_sub = T_lev * steam_dlnTdlnP_sat(T_lev, p_lev, Rd, cpdry_new) &
                         * dlogp_sub
              else
                Gm      = malr(T_lev, p_lev, Rd, exo_g, cpdry_new)
                es_sub  = min(esat(T_lev), 0.99_r8 * p_lev)
                r_sub   = eps_wv * es_sub / (p_lev - es_sub)
                T_v_lev = T_lev * (1._r8 + r_sub / eps_wv) / (1._r8 + r_sub)
                dT_sub  = Gm * Rd * T_v_lev / exo_g * dlogp_sub
              end if
              T_prev_sub = T_lev;  p_prev_sub = p_lev
              T_lev   = T_lev + dT_sub
              p_lev   = p_lev * exp(dlogp_sub)
              if (T_lev <= t_strato) then
                if (p_coldpoint < 0._r8) p_coldpoint = &
                  interp_logp(T_prev_sub, p_prev_sub, T_lev, p_lev, t_strato)
                T_lev = t_strato; exit
              end if
            end do
            T_at_int(k) = T_lev
          else
            if (use_nonideal_adiabat) then
              ! Non-ideal unsaturated (constant-composition) adiabat, Kasting
              ! A11/A12: integrate one log-pressure layer up from the level below.
              T_at_int(k) = T_at_int(k+1) + T_at_int(k+1) * &
                steam_dlnTdlnP_dry(T_at_int(k+1), pint(k+1), alfav_dry, Rd, cpdry_new) * &
                log(pint(k) / pint(k+1))
            else
              ! Dry adiabat: T(p) = Ts · (p/ps)^κ_mix  (exact for ideal gas)
              T_at_int(k) = cs_ts * (pint(k) / ps_use)**kappa_mix
            end if
            if (T_at_int(k) <= t_strato) T_at_int(k) = t_strato
            ! Top of the dry region = where the constant-mixing-ratio parcel first
            ! reaches saturation (Kasting 1988): its vapour partial pressure
            ! pH2O(level) = f_vmr·p  meets/exceeds esat(T).  Because the dry adiabat
            ! is steeper than the SVP curve, this happens at some height above the
            ! (subsaturated) surface and marks the switch to the moist adiabat.
            if (f_vmr * pint(k) >= esat(T_at_int(k))) then
              k_cond       = k
              in_dry_layer = .false.
              write(*,'(a,es10.3,a,f7.2,a)') &
                '    IHZ condensation level: p = ', pint(k)/100._r8, &
                ' hPa  T = ', T_at_int(k), ' K  (switch to moist adiabat)'
            end if
          end if
        end do

      else
        ! Surface saturated (esat(Ts) ≤ inventory): moist pseudoadiabat from the
        ! surface — the moist-greenhouse branch (Kasting 1988, Fig 1a).  No dry base.
        write(*,'(a)') &
          '    IHZ: surface saturated — moist pseudoadiabat from the surface'
        k_cond = 0
        do k = pverp-1, 1, -1
          if (T_at_int(k+1) <= t_strato) then
            T_at_int(k) = t_strato; cycle
          end if
          T_lev       = T_at_int(k+1)
          p_lev       = pint(k+1)
          dlogp_layer = log(pint(k) / pint(k+1))
          dlogp_sub   = dlogp_layer / real(NSUB, r8)
          do isub = 1, NSUB
            if (use_nonideal_adiabat .and. T_lev >= IAPWS_TT) then
              ! Kasting A4/A5 non-ideal saturated moist pseudoadiabat (liquid-
              ! vapour SVP; below the triple point fall back to the ideal malr)
              dT_sub = T_lev * steam_dlnTdlnP_sat(T_lev, p_lev, Rd, cpdry_new) &
                       * dlogp_sub
            else
              Gm      = malr(T_lev, p_lev, Rd, exo_g, cpdry_new)
              es_sub  = min(esat(T_lev), 0.99_r8 * p_lev)
              r_sub   = eps_wv * es_sub / (p_lev - es_sub)
              T_v_lev = T_lev * (1._r8 + r_sub / eps_wv) / (1._r8 + r_sub)
              dT_sub  = Gm * Rd * T_v_lev / exo_g * dlogp_sub
            end if
            T_prev_sub = T_lev;  p_prev_sub = p_lev
            T_lev   = T_lev + dT_sub
            p_lev   = p_lev * exp(dlogp_sub)
            if (T_lev <= t_strato) then
              if (p_coldpoint < 0._r8) p_coldpoint = &
                interp_logp(T_prev_sub, p_prev_sub, T_lev, p_lev, t_strato)
              T_lev = t_strato; exit
            end if
          end do
          T_at_int(k) = T_lev
        end do
      end if

    else
      ! Standard: moist pseudoadiabat from surface, cap at t_strato
      do k = pverp-1, 1, -1
        if (T_at_int(k+1) <= t_strato) then
          T_at_int(k) = t_strato
          cycle
        end if
        T_lev       = T_at_int(k+1)
        p_lev       = pint(k+1)
        dlogp_layer = log(pint(k) / pint(k+1))
        dlogp_sub   = dlogp_layer / real(NSUB, r8)
        do isub = 1, NSUB
          if (use_nonideal_adiabat .and. T_lev >= IAPWS_TT) then
            ! Kasting A4/A5 non-ideal saturated moist pseudoadiabat (liquid-
            ! vapour SVP; below the triple point fall back to the ideal malr)
            dT_sub = T_lev * steam_dlnTdlnP_sat(T_lev, p_lev, Rd, cpdry_new) &
                     * dlogp_sub
          else
            Gm      = malr(T_lev, p_lev, Rd, exo_g, cpdry_new)
            es_sub  = min(esat(T_lev), 0.99_r8 * p_lev)
            r_sub   = eps_wv * es_sub / (p_lev - es_sub)
            T_v_lev = T_lev * (1._r8 + r_sub / eps_wv) / (1._r8 + r_sub)
            dT_sub  = Gm * Rd * T_v_lev / exo_g * dlogp_sub
          end if
          T_prev_sub = T_lev;  p_prev_sub = p_lev
          T_lev   = T_lev + dT_sub
          p_lev   = p_lev * exp(dlogp_sub)
          if (T_lev <= t_strato) then
            if (p_coldpoint < 0._r8) p_coldpoint = &
              interp_logp(T_prev_sub, p_prev_sub, T_lev, p_lev, t_strato)
            T_lev = t_strato
            exit
          end if
        end do
        T_at_int(k) = T_lev
      end do
    end if

    tint(:) = T_at_int(:)
    do k = 1, pver
      tmid(k) = 0.5_r8 * (T_at_int(k) + T_at_int(k+1))
    end do

    write(*,'(a,f6.1,a,f6.1,a)') &
      '    T profile  : ts = ', cs_ts, ' K  →  TOA tint(1) = ', tint(1), ' K'

    ! ---- 4. Moisture: Kasting (1988) cold-trap convention ----
    ! Find the topmost tropospheric layer (first from TOA where T > t_strato).
    k_top_conv = 0
    do k = 1, pver
      if (tmid(k) > t_strato) then
        k_top_conv = k
        exit
      end if
    end do

    ! Stratospheric H2O: qsat(t_strato, p_tropo), uniform in the stratosphere.
    ! p_tropo is the TRUE cold-point pressure where the continuous adiabat crosses
    ! t_strato (interpolated within the integration substep, p_coldpoint), NOT the
    ! snapped interface pint(k_top_conv).  Using the snapped interface made the
    ! cold-trap mixing ratio jump discontinuously each time the integer tropopause
    ! index stepped between layers as Ts swept — the source of the staircase in the
    ! inner-HZ OLR/albedo/Seff sweeps.  Fall back to the snapped interface only if
    ! the cold point was never reached (e.g. cold adiabat captured in the dry layer).
    if (k_top_conv > 0) then
      es_tropo    = esat(t_strato)
      if (p_coldpoint > 0._r8) then
        p_tropo_use = p_coldpoint
      else
        p_tropo_use = pint(k_top_conv)
      end if
      q_cold_trap = eps_wv * es_tropo / (p_tropo_use - es_tropo)
      write(*,'(a,f9.2,a)') &
        '    cold-point p (interp) : ', p_tropo_use/100._r8, ' hPa'
    else
      q_cold_trap = 0._r8   ! isothermal column: no stratospheric H2O
    end if

    do k = 1, pver
      if (tmid(k) > t_strato) then
        ! In the dry layer (subsaturated-surface runaway branch), h2ommr is held
        ! constant at the inventory-set value.  In the moist layer, h2ommr = rh·qsat.
        if (ihz_profile .and. (.not. surface_saturated) .and. k >= k_cond) then
          h2ommr(k) = h2ommr_dry
        else
          es_k      = min(esat(tmid(k)), 0.99_r8 * pmid(k))
          qsat_k    = eps_wv * es_k / (pmid(k) - es_k)
          h2ommr(k) = rh_init * qsat_k
        end if
      else
        h2ommr(k) = rh_init * q_cold_trap
      end if
    end do

    ! ---- 5. Dry-gas MMR profiles (well-mixed, except O3 per o3_profile) ----
    co2mmr(:)  = mmr_dry(1)
    ch4mmr(:)  = mmr_dry(2)
    c2h6mmr(:) = mmr_dry(3)
    h2mmr(:)   = mmr_dry(4)
    n2mmr(:)   = mmr_dry(5)
    select case (trim(o3_profile))
    case ('none')
      o3mmr(:) = 0.0_r8
    case ('earth')
      call set_earth_o3_profile(pmid, mwdry_new, o3mmr)
    case ('rcemip')
      call set_rcemip_o3_profile(pmid, mwdry_new, o3mmr)
    case default  ! 'uniform'
      o3mmr(:) = mmr_dry(6)
    end select
    o2mmr(:)   = mmr_dry(7)

    ! ---- 6. Hydrostatic heights ----
    zint(pverp) = 0._r8
    do k = pver, 1, -1
      zint(k) = zint(k+1) + (Rd / exo_g) * &
                tmid(k) * (1._r8 + h2ommr(k) / eps_wv) / (1._r8 + h2ommr(k)) * &
                log(pint(k+1) / pint(k))
    end do

    ! ---- 7. Surface + scalar state ----
    ts        = cs_ts
    ps        = ps_use
    mwdry_col = mwdry_new
    cpdry_col = cpdry_new
    coszrs    = cs_coszrs
    asdir     = cs_asdir
    asdif     = cs_asdif
    aldir     = cs_aldir
    aldif     = cs_aldif
    msdist    = cfg_msdist

    ! Clear-sky defaults
    cfrc(:)   = 0._r8
    cicewp(:) = 0._r8
    cliqwp(:) = 0._r8
    rei(:)    = 0._r8
    rel(:)    = 0._r8

    ! Refresh derived dry-pressure arrays now that pdel, pint, h2ommr are set.
    call exocol_update_derived()

    write(*,'(a)') '  cold_start_init: initial column ready.'

  end subroutine cold_start_init

  ! -----------------------------------------------------------------------

  subroutine build_pressure_grid(ps_use)
  ! Build the full pint/pmid/pdel arrays by pure log-spacing from p_top to ps
  ! (lowest level ~400 m at pver=70).  The surface-coupled mixed layer
  ! (convadj_surface) handles the near-surface coupling, so no fine near-surface
  ! grid is needed.
    real(r8), intent(in) :: ps_use

    integer :: k

    do k = 1, pverp
      pint(k) = exp(log(p_top) + &
                    real(k-1, r8) / real(pverp-1, r8) * log(ps_use / p_top))
    end do

    ! Midpoints and layer thicknesses
    do k = 1, pver
      pmid(k) = 0.5_r8 * (pint(k) + pint(k+1))
      pdel(k) = pint(k+1) - pint(k)
    end do

  end subroutine build_pressure_grid

  ! -----------------------------------------------------------------------

  pure function interp_logp(T1, p1, T2, p2, Tx) result(px)
  ! Log-pressure interpolation: the pressure at which temperature = Tx along the
  ! segment (T1,p1) -> (T2,p2).  Used to locate the true cold-point pressure
  ! (where the continuous adiabat crosses t_strato) so the cold-trap water is a
  ! smooth function of Ts instead of snapping to a discrete grid interface.
    real(r8), intent(in) :: T1, p1, T2, p2, Tx
    real(r8) :: px, f
    f  = (T1 - Tx) / (T1 - T2)
    px = exp(log(p1) + f * (log(p2) - log(p1)))
  end function interp_logp

end module exocol_coldstart
