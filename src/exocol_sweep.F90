module exocol_sweep
! Flux sweep over surface temperature, solar zenith angle and surface albedo.
!
! Purpose: build radiation lookup tables (OLR and planetary albedo) for
! energy-balance models, which need the TOA fluxes of a prescribed profile as
! a function of those three quantities.  Everything here is about paying the
! ExoRT initialisation once:
!
!   * neither the zenith angle nor the surface albedo feeds back on the
!     cold-start profile — the moist adiabat is built from ts, p_top and the
!     composition alone — so the inner (zenith x albedo) block reuses one
!     column and costs only one aerad_driver call per point (~0.05 s at
!     pver = 200) against ~1.6 s of fixed initialisation;
!   * the optional outer loop over sweep_ts DOES rebuild the column
!     (cold_start_init), but still inside the same process, so a whole
!     temperature axis is covered for one initialisation.  This matters more
!     than the arithmetic suggests: a single-point run carries a ~240 MB
!     working set, and repaying that per grid point makes a large table
!     memory-bandwidth bound rather than compute bound.
!
! Activated by &exocol_sweep::sweep_mode (which forces flux_only); see
! exocol_config for the namelist.  Output is a plain-text table written to
! sweep_outfile: a '#'-prefixed header carrying the run-wide scalars, then one
! record per (Ts, coszrs, albedo) combination, each carrying its own column
! scalars (the total surface pressure and surface humidity change with Ts
! under variable_ps).
!
! The longwave is independent of both swept quantities, so the OLR column is
! constant by construction; it is written on every record anyway so that a
! consumer can assert the invariance as a sanity check.

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pver, pverp
  use exocol_mod,      only: ts, ps, coszrs, asdir, asdif, aldir, aldif, &
                             msdist, mwdry_col, cpdry_col, h2ommr,     &
                             exocol_setgas
  use exocol_config,   only: sweep_coszrs, sweep_albedo, MAX_SWEEP, &
                             sweep_ts, sweep_t_strato,               &
                             cfg_ps => ps, co2_vmr, n2_vmr,          &
                             cfg_ts => ts, cfg_t_strato => t_strato
  use exocol_radiation, only: exocol_rad_tend
  use exocol_coldstart, only: cold_start_init

  implicit none
  private

  public :: run_flux_sweep

contains

  subroutine run_flux_sweep(outfile)
  ! Loop the radiation call over every listed (coszrs, albedo) combination and
  ! write the TOA fluxes.  The column state must already be built and ExoRT
  ! initialised; on return coszrs and the four albedos are left at the last
  ! swept values (the driver's subsequent write_output therefore reports that
  ! last combination, which is intentional — it makes the NetCDF file a
  ! profile record for the sweep, not an independent state).
    character(len=*), intent(in) :: outfile

    real(r8), dimension(pver)  :: LWHR, SWHR
    real(r8), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

    real(r8) :: mu_list(MAX_SWEEP), alb_list(MAX_SWEEP)
    real(r8) :: ts_list(MAX_SWEEP), tstrat_list(MAX_SWEEP)
    integer  :: nmu, nalb, nts, i, j, k, unit, ios
    real(r8) :: olr, sdn, sup, palb
    logical  :: loop_ts

    ! Collect the set entries (sentinel -1 marks unused list slots).
    nmu = 0
    do i = 1, MAX_SWEEP
      if (sweep_coszrs(i) >= 0.0_r8) then
        nmu = nmu + 1
        mu_list(nmu) = sweep_coszrs(i)
      end if
    end do

    nalb = 0
    do i = 1, MAX_SWEEP
      if (sweep_albedo(i) >= 0.0_r8) then
        nalb = nalb + 1
        alb_list(nalb) = sweep_albedo(i)
      end if
    end do

    if (nmu == 0 .or. nalb == 0) then
      write(*,'(a)') '  exocol_sweep: sweep_coszrs or sweep_albedo is empty — nothing to do.'
      return
    end if

    ! Optional outer temperature loop.  With no sweep_ts the column is left
    ! exactly as the driver built it (bit-identical to a single-Ts sweep).
    nts = 0
    do i = 1, MAX_SWEEP
      if (sweep_ts(i) >= 0.0_r8) then
        nts = nts + 1
        ts_list(nts) = sweep_ts(i)
        if (sweep_t_strato(i) >= 0.0_r8) then
          tstrat_list(nts) = sweep_t_strato(i)
        else
          tstrat_list(nts) = cfg_t_strato
        end if
      end if
    end do
    loop_ts = (nts > 0)
    if (.not. loop_ts) then
      nts = 1
      ts_list(1)     = ts
      tstrat_list(1) = cfg_t_strato
    end if

    open(newunit=unit, file=trim(outfile), status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) then
      write(*,'(3a)') '  exocol_sweep: cannot open ', trim(outfile), ' for writing.'
      stop 1
    end if

    ! Header: everything a table builder needs to key this column, including
    ! the total surface pressure actually used (which exceeds the requested dry
    ! pressure when variable_ps adds the surface water vapour on top).
    write(unit,'(a)')             '# ExoColumn flux sweep (Ts x zenith x surface albedo)'
    write(unit,'(a,es16.8)')      '# ps_dry_Pa   = ', cfg_ps
    write(unit,'(a,es16.8)')      '# co2_vmr     = ', co2_vmr
    write(unit,'(a,es16.8)')      '# n2_vmr      = ', n2_vmr
    write(unit,'(a,es16.8)')      '# msdist      = ', msdist
    write(unit,'(a,i0,a,i0,a,i0)') '# nts         = ', nts, '   nmu = ', nmu, &
                                   '   nalb = ', nalb
    write(unit,'(a)')             '# ts_K  ps_total_Pa  mwdry  h2ommr_srf  ' //  &
                                  'coszrs  surfalb  olr_Wm2  swdn_toa_Wm2  ' //  &
                                  'swup_toa_Wm2  palb'

    do k = 1, nts

      if (loop_ts) then
        ! Rebuild the prescribed column at this surface temperature.  Only the
        ! column state changes — the ExoRT initialisation stands.
        cfg_ts       = ts_list(k)
        cfg_t_strato = tstrat_list(k)
        call cold_start_init()
        call exocol_setgas()
      end if

      do i = 1, nmu
        coszrs = mu_list(i)
        do j = 1, nalb
          ! HEXTOR-style grey surface: one broadband albedo in all four ExoRT
          ! channels (visible/near-IR x direct/diffuse).
          asdir = alb_list(j)
          asdif = alb_list(j)
          aldir = alb_list(j)
          aldif = alb_list(j)

          call exocol_rad_tend(LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)

          olr = LWUP(1)
          sdn = SWDN(1)
          sup = SWUP(1)
          if (sdn > 0.0_r8) then
            palb = sup / sdn
          else
            ! No incident beam (coszrs = 0): the planetary albedo is undefined.
            ! Flag it rather than emit 0/0; the table builder must not use it.
            palb = -1.0_r8
          end if

          write(unit,'(f9.3,es16.8,f9.4,es16.8,2f9.5,4es20.10)')                &
                ts, ps, mwdry_col, h2ommr(pver),                                &
                mu_list(i), alb_list(j), olr, sdn, sup, palb
        end do
      end do

    end do

    close(unit)

    write(*,'(a,i0,a,a)') '  exocol_sweep: wrote ', nts*nmu*nalb, &
                          ' flux records to ', trim(outfile)

  end subroutine run_flux_sweep

end module exocol_sweep
