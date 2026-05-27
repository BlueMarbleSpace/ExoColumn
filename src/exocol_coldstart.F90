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
!      Euler in log(p):  dT/d(ln p) = Γm · R_d · T / g.
!   4. Moisture: h2ommr(k) = rh_init · qsat(T(k), p(k)) where T > t_strato;
!      zero in the stratosphere.
!   5. Heights: hydrostatic integration  zint(k) = zint(k+1) + R_d T /g · ln(p(k+1)/p(k)).
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
                            n_sfc_layers, dp_sfc_bot, sfc_stretch, &
                            co2_vmr, ch4_vmr, c2h6_vmr,           &
                            h2_vmr,  n2_vmr,  o3_vmr, o2_vmr,     &
                            o3_profile,                            &
                            MW_CO2, MW_CH4, MW_C2H6, MW_H2,       &
                            MW_N2,  MW_O3,  MW_O2,                 &
                            CP_CO2, CP_CH4, CP_C2H6, CP_H2,       &
                            CP_N2,  CP_O3,  CP_O2
  use exocol_convadj, only: esat_cc, malr
  use exocol_ozone,   only: set_earth_o3_profile

  implicit none
  private

  public :: cold_start_init

  ! Earth-like dry-air defaults used when a gas VMR is not specified in
  ! &exocol_composition.  These sum to ~0.99 (the missing ~1% corresponds to
  ! Ar/other inert gases not modelled by ExoRT; not enough to affect mwdry).
  real(r8), parameter :: DEF_CO2  = 4.0e-4_r8
  real(r8), parameter :: DEF_CH4  = 0.0_r8
  real(r8), parameter :: DEF_C2H6 = 0.0_r8
  real(r8), parameter :: DEF_H2   = 0.0_r8
  real(r8), parameter :: DEF_N2   = 0.78_r8
  real(r8), parameter :: DEF_O3   = 0.0_r8
  real(r8), parameter :: DEF_O2   = 0.21_r8

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
    real(r8) :: T_at_int(pverp)
    real(r8) :: T_lev, p_lev, dlogp_layer, dlogp_sub, dT_sub, Gm
    real(r8) :: es_k, qsat_k
    character(len=4) :: dry_name(7)

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

    vmr_sum = sum(vmr_dry)
    if (abs(vmr_sum - 1.0_r8) > 1.0e-2_r8) then
      write(*,'(a,f8.5,a)') &
        '  WARNING: cold-start dry-air VMRs sum to ', vmr_sum, &
        ' (expected ~1.0). Values used as-is.'
    end if

    mwdry_new = sum(vmr_dry * Mw_dry)
    if (mwdry_new <= 0._r8) then
      write(*,'(a)') '  cold_start_init: computed mwdry <= 0. Aborting.'
      stop 'cold_start_init: bad composition'
    end if

    do k = 1, 7
      mmr_dry(k) = vmr_dry(k) * Mw_dry(k) / mwdry_new
    end do

    ! Mass-weighted mean cp from the actual composition.  This replaces the
    ! namelist cpdry so non-Earth compositions (especially H2-rich) get the
    ! correct specific heat automatically.
    cpdry_new = sum(mmr_dry * cp_dry)

    write(*,'(a,f7.3,a)') '    mwdry  = ', mwdry_new, ' g/mol'
    write(*,'(a,f7.1,a)') '    cpdry  = ', cpdry_new, ' J/kg/K (auto-computed from composition)'
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

    ! ---- 3. Temperature: moist adiabat from surface, capped at t_strato ----
    Rd     = SHR_CONST_RGAS / mwdry_new
    eps_wv = SHR_CONST_MWWV / mwdry_new

    T_at_int(pverp) = cs_ts

    if (cs_ts <= t_strato) then
      ! ts at or below stratosphere: column is isothermal at t_strato above
      ! the surface (ts is preserved at pverp as the skin temperature).
      write(*,'(a)') &
        '  cold_start_init: ts <= t_strato — isothermal column above surface.'
      T_at_int(1:pverp-1) = t_strato
    else
      do k = pverp-1, 1, -1
        if (T_at_int(k+1) <= t_strato) then
          T_at_int(k) = t_strato
          cycle
        end if
        T_lev       = T_at_int(k+1)
        p_lev       = pint(k+1)
        dlogp_layer = log(pint(k) / pint(k+1))    ! < 0  (going up)
        dlogp_sub   = dlogp_layer / real(NSUB, r8)
        do isub = 1, NSUB
          Gm     = malr(T_lev, p_lev, Rd, exo_g, cpdry_new)
          dT_sub = Gm * Rd * T_lev / exo_g * dlogp_sub   ! < 0
          T_lev  = T_lev + dT_sub
          p_lev  = p_lev * exp(dlogp_sub)
          if (T_lev <= t_strato) then
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

    ! ---- 4. Moisture (rh_init · qsat in troposphere, dry stratosphere) ----
    do k = 1, pver
      if (tmid(k) > t_strato) then
        es_k      = min(esat_cc(tmid(k)), 0.99_r8 * pmid(k))
        qsat_k    = eps_wv * es_k / (pmid(k) - es_k)
        h2ommr(k) = rh_init * qsat_k
      else
        h2ommr(k) = 0._r8
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
    case default  ! 'uniform'
      o3mmr(:) = mmr_dry(6)
    end select
    o2mmr(:)   = mmr_dry(7)

    ! ---- 6. Hydrostatic heights ----
    zint(pverp) = 0._r8
    do k = pver, 1, -1
      zint(k) = zint(k+1) + (Rd / exo_g) * tmid(k) * log(pint(k+1) / pint(k))
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
  ! Build the full pint/pmid/pdel arrays.
  !
  ! When n_sfc_layers > 0 (hybrid mode):
  !   Bottom n_sfc_layers interfaces are built upward from ps using a
  !   geometric series: pdel(pver) = dp_sfc_bot, pdel(pver-1) = dp*sfc_stretch,
  !   etc.  The remaining pver-n_sfc_layers upper layers are log-spaced from
  !   p_top down to p_match (the top of the surface block).
  !
  ! When n_sfc_layers = 0: pure log-spacing from p_top to ps (legacy).
    real(r8), intent(in) :: ps_use

    integer  :: k, n_upper, n_sfc
    real(r8) :: p_match, dp_k

    n_sfc = n_sfc_layers

    ! Validate: surface layers must not exceed the total column.
    if (n_sfc >= pver) then
      write(*,'(a,i0,a,i0,a)') &
        '  WARNING: n_sfc_layers (', n_sfc, &
        ') >= pver (', pver, ') — reverting to pure log-spacing'
      n_sfc = 0
    end if

    if (n_sfc == 0) then
      ! Pure log-spacing (legacy / n_sfc_layers = 0)
      do k = 1, pverp
        pint(k) = exp(log(p_top) + &
                      real(k-1, r8) / real(pverp-1, r8) * log(ps_use / p_top))
      end do

    else
      ! --- Surface block: geometric layers from ps upward ---
      pint(pverp) = ps_use
      do k = pver, pver - n_sfc + 1, -1
        dp_k   = dp_sfc_bot * sfc_stretch**(pver - k)
        pint(k) = pint(k+1) - dp_k
      end do
      p_match = pint(pver - n_sfc + 1)

      if (p_match <= p_top) then
        write(*,'(a,es10.3,a)') &
          '  WARNING: surface layers span above p_top (p_match = ', p_match, &
          ' Pa) — reduce n_sfc_layers or dp_sfc_bot. Reverting to log-spacing.'
        do k = 1, pverp
          pint(k) = exp(log(p_top) + &
                        real(k-1, r8) / real(pverp-1, r8) * log(ps_use / p_top))
        end do

      else
        ! --- Upper block: log-spaced from p_top to p_match ---
        n_upper = pver - n_sfc
        do k = 1, n_upper + 1
          pint(k) = exp(log(p_top) + &
                        real(k-1, r8) / real(n_upper, r8) * log(p_match / p_top))
        end do
        ! pint(n_upper+1) is written twice (log-spaced then surface block both
        ! set it to p_match); values are identical, no conflict.

      end if
    end if

    ! Midpoints and layer thicknesses
    do k = 1, pver
      pmid(k) = 0.5_r8 * (pint(k) + pint(k+1))
      pdel(k) = pint(k+1) - pint(k)
    end do

  end subroutine build_pressure_grid

end module exocol_coldstart
