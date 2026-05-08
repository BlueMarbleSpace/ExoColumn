module exocol_io
! Input/output for ExoColumn using the NetCDF-Fortran 90 interface.
!
! read_initial_conditions: populates exocol_mod arrays from RTprofile_in.nc
!   (or any ExoRT-format input file).  Validates that the file's pver matches
!   the compile-time constant from ppgrid.
!
! write_output: writes the equilibrium column state to RTprofile_out.nc in
!   ExoRT-compatible format so existing ExoRT plotting tools work unchanged.
!   Heating rates are written in K/day (our internal unit).

  use netcdf
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pverp
  use exocol_mod

  implicit none
  private

  public :: read_initial_conditions
  public :: write_output

contains

  ! -----------------------------------------------------------------------
  ! Internal helpers
  ! -----------------------------------------------------------------------

  subroutine nc_check(status, context)
    integer,          intent(in) :: status
    character(len=*), intent(in) :: context
    if (status /= nf90_noerr) then
      write(*,'(a,a,a)') 'NetCDF error in ', trim(context), ':'
      write(*,'(a)')     nf90_strerror(status)
      stop 'exocol_io: NetCDF error'
    end if
  end subroutine nc_check

  subroutine get_scalar(ncid, name, val)
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: val
    integer :: varid
    call nc_check(nf90_inq_varid(ncid, name, varid), 'inq_varid:'//name)
    call nc_check(nf90_get_var(ncid, varid, val),    'get_var:'//name)
  end subroutine get_scalar

  subroutine get_1d(ncid, name, arr)
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: arr(:)
    integer :: varid
    call nc_check(nf90_inq_varid(ncid, name, varid), 'inq_varid:'//name)
    call nc_check(nf90_get_var(ncid, varid, arr),    'get_var:'//name)
  end subroutine get_1d

  subroutine get_1d_opt(ncid, name, arr)
  ! Like get_1d but silently leaves arr unchanged if the variable is absent.
  ! Used for gases and cloud fields not present in all input files.
    integer,          intent(in)    :: ncid
    character(len=*), intent(in)    :: name
    real(r8),         intent(inout) :: arr(:)
    integer :: varid, status
    status = nf90_inq_varid(ncid, name, varid)
    if (status == nf90_noerr) then
      call nc_check(nf90_get_var(ncid, varid, arr), 'get_var:'//name)
    end if
  end subroutine get_1d_opt

  subroutine put_scalar(ncid, varid, val, context)
    integer,          intent(in) :: ncid, varid
    real(r8),         intent(in) :: val
    character(len=*), intent(in) :: context
    call nc_check(nf90_put_var(ncid, varid, val), context)
  end subroutine put_scalar

  subroutine put_1d(ncid, varid, arr, context)
    integer,          intent(in) :: ncid, varid
    real(r8),         intent(in) :: arr(:)
    character(len=*), intent(in) :: context
    call nc_check(nf90_put_var(ncid, varid, arr), context)
  end subroutine put_1d

  ! -----------------------------------------------------------------------

  subroutine read_initial_conditions(filename)
  ! Read column state from an ExoRT RTprofile_in.nc file.
  ! Populates all exocol_mod arrays.  cfrc is set to zero (clear sky);
  ! msdist is set to 1.0 (1 AU equivalent).
    character(len=*), intent(in) :: filename

    integer :: ncid, dimid
    integer :: npver_file, npverp_file

    write(*,'(2a)') 'exocol_io: reading ', trim(filename)
    call nc_check(nf90_open(filename, nf90_nowrite, ncid), 'open:'//filename)

    ! Validate grid dimensions against compile-time constants
    call nc_check(nf90_inq_dimid(ncid, 'pver', dimid),          'inq_dimid:pver')
    call nc_check(nf90_inquire_dimension(ncid, dimid, len=npver_file), 'inq_dim:pver')
    if (npver_file /= pver) then
      write(*,'(a,i0,a,i0)') 'ERROR: file pver=', npver_file, &
                              ' but compiled pver=', pver
      stop 'exocol_io: pver mismatch'
    end if

    call nc_check(nf90_inq_dimid(ncid, 'pverp', dimid),             'inq_dimid:pverp')
    call nc_check(nf90_inquire_dimension(ncid, dimid, len=npverp_file), 'inq_dim:pverp')
    if (npverp_file /= pverp) then
      write(*,'(a,i0,a,i0)') 'ERROR: file pverp=', npverp_file, &
                              ' but compiled pverp=', pverp
      stop 'exocol_io: pverp mismatch'
    end if

    ! Scalars
    call get_scalar(ncid, 'ts',     ts)
    call get_scalar(ncid, 'ps',     ps)
    call get_scalar(ncid, 'coszrs', coszrs)
    call get_scalar(ncid, 'mw',     mwdry_col)
    call get_scalar(ncid, 'cp',     cpdry_col)
    call get_scalar(ncid, 'asdir',  asdir)
    call get_scalar(ncid, 'asdif',  asdif)
    call get_scalar(ncid, 'aldir',  aldir)
    call get_scalar(ncid, 'aldif',  aldif)

    ! Layer-midpoint arrays
    call get_1d(ncid, 'tmid',   tmid)
    call get_1d(ncid, 'pmid',   pmid)
    call get_1d(ncid, 'pdel',   pdel)
    call get_1d(ncid, 'h2ommr', h2ommr)
    call get_1d(ncid,     'co2mmr', co2mmr)
    call get_1d(ncid,     'ch4mmr', ch4mmr)
    call get_1d_opt(ncid, 'c2h6mmr',c2h6mmr)  ! absent in most standard profiles
    call get_1d(ncid,     'o2mmr',  o2mmr)
    call get_1d(ncid,     'o3mmr',  o3mmr)
    call get_1d(ncid,     'n2mmr',  n2mmr)
    call get_1d(ncid,     'h2mmr',  h2mmr)
    call get_1d_opt(ncid, 'cicewp', cicewp)    ! absent in clear-sky profiles
    call get_1d_opt(ncid, 'cliqwp', cliqwp)
    call get_1d_opt(ncid, 'rei',    rei)
    call get_1d_opt(ncid, 'rel',    rel)

    ! Interface arrays
    call get_1d(ncid, 'tint', tint)
    call get_1d(ncid, 'pint', pint)
    call get_1d(ncid, 'zint', zint)

    call nc_check(nf90_close(ncid), 'close:'//filename)

    ! Defaults not stored in the input file
    cfrc(:) = 0._r8   ! clear sky
    msdist  = 1._r8   ! 1 AU; override via namelist if needed

    write(*,'(a,i0,a)') 'exocol_io: read ', pver, ' layers OK'

  end subroutine read_initial_conditions

  ! -----------------------------------------------------------------------

  subroutine write_output(filename, LWHR, SWHR, LWUP, LWDN, SWUP, SWDN)
  ! Write equilibrium column state to NetCDF in ExoRT RTprofile_out.nc format.
  ! LWHR / SWHR are in K/day.  All other variables from exocol_mod state.
    character(len=*), intent(in) :: filename
    real(r8), intent(in), dimension(pver)  :: LWHR, SWHR
    real(r8), intent(in), dimension(pverp) :: LWUP, LWDN, SWUP, SWDN

    integer :: ncid
    integer :: pver_dim, pverp_dim, one_dim
    integer :: lwup_id, lwdn_id, swup_id, swdn_id
    integer :: lwhr_id, swhr_id
    integer :: ts_id,  ps_id
    integer :: tmid_id, tint_id, pmid_id, pint_id, zint_id
    integer :: h2o_id, co2_id, ch4_id, c2h6_id
    integer :: o2_id,  o3_id,  n2_id,  h2_id
    integer :: mw_id,  cp_id

    write(*,'(2a)') 'exocol_io: writing ', trim(filename)
    call nc_check(nf90_create(filename, nf90_clobber, ncid), 'create:'//filename)

    ! ---- Dimensions ----
    call nc_check(nf90_def_dim(ncid, 'pver',  pver,  pver_dim),  'def_dim:pver')
    call nc_check(nf90_def_dim(ncid, 'pverp', pverp, pverp_dim), 'def_dim:pverp')
    call nc_check(nf90_def_dim(ncid, 'ONE',   1,     one_dim),   'def_dim:ONE')

    ! ---- Flux variables (pverp) ----
    call nc_check(nf90_def_var(ncid,'LWUP',nf90_double,[pverp_dim],lwup_id),'def_var:LWUP')
    call nc_check(nf90_put_att(ncid,lwup_id,'long_name','LW upwelling flux'),   'att')
    call nc_check(nf90_put_att(ncid,lwup_id,'units','W/m2'),                    'att')

    call nc_check(nf90_def_var(ncid,'LWDN',nf90_double,[pverp_dim],lwdn_id),'def_var:LWDN')
    call nc_check(nf90_put_att(ncid,lwdn_id,'long_name','LW downwelling flux'), 'att')
    call nc_check(nf90_put_att(ncid,lwdn_id,'units','W/m2'),                   'att')

    call nc_check(nf90_def_var(ncid,'SWUP',nf90_double,[pverp_dim],swup_id),'def_var:SWUP')
    call nc_check(nf90_put_att(ncid,swup_id,'long_name','SW upwelling flux'),   'att')
    call nc_check(nf90_put_att(ncid,swup_id,'units','W/m2'),                    'att')

    call nc_check(nf90_def_var(ncid,'SWDN',nf90_double,[pverp_dim],swdn_id),'def_var:SWDN')
    call nc_check(nf90_put_att(ncid,swdn_id,'long_name','SW downwelling flux'), 'att')
    call nc_check(nf90_put_att(ncid,swdn_id,'units','W/m2'),                   'att')

    ! ---- Heating rates (pver, K/day) ----
    call nc_check(nf90_def_var(ncid,'LWHR',nf90_double,[pver_dim],lwhr_id),'def_var:LWHR')
    call nc_check(nf90_put_att(ncid,lwhr_id,'long_name','LW heating rate'), 'att')
    call nc_check(nf90_put_att(ncid,lwhr_id,'units','K/day'),               'att')

    call nc_check(nf90_def_var(ncid,'SWHR',nf90_double,[pver_dim],swhr_id),'def_var:SWHR')
    call nc_check(nf90_put_att(ncid,swhr_id,'long_name','SW heating rate'), 'att')
    call nc_check(nf90_put_att(ncid,swhr_id,'units','K/day'),               'att')

    ! ---- Surface scalars ----
    call nc_check(nf90_def_var(ncid,'ts',nf90_double,[one_dim],ts_id),'def_var:ts')
    call nc_check(nf90_put_att(ncid,ts_id,'long_name','surface temperature'),'att')
    call nc_check(nf90_put_att(ncid,ts_id,'units','K'),                      'att')

    call nc_check(nf90_def_var(ncid,'ps',nf90_double,[one_dim],ps_id),'def_var:ps')
    call nc_check(nf90_put_att(ncid,ps_id,'long_name','surface pressure'),'att')
    call nc_check(nf90_put_att(ncid,ps_id,'units','Pa'),                  'att')

    ! ---- T, P, Z profiles ----
    call nc_check(nf90_def_var(ncid,'tmid',nf90_double,[pver_dim], tmid_id),'def_var:tmid')
    call nc_check(nf90_put_att(ncid,tmid_id,'long_name','midlayer temperatures'),'att')
    call nc_check(nf90_put_att(ncid,tmid_id,'units','K'),                        'att')

    call nc_check(nf90_def_var(ncid,'tint',nf90_double,[pverp_dim],tint_id),'def_var:tint')
    call nc_check(nf90_put_att(ncid,tint_id,'long_name','interface temperatures'),'att')
    call nc_check(nf90_put_att(ncid,tint_id,'units','K'),                         'att')

    call nc_check(nf90_def_var(ncid,'pmid',nf90_double,[pver_dim], pmid_id),'def_var:pmid')
    call nc_check(nf90_put_att(ncid,pmid_id,'long_name','midlayer pressures'),'att')
    call nc_check(nf90_put_att(ncid,pmid_id,'units','Pa'),                    'att')

    call nc_check(nf90_def_var(ncid,'pint',nf90_double,[pverp_dim],pint_id),'def_var:pint')
    call nc_check(nf90_put_att(ncid,pint_id,'long_name','interface pressures'),'att')
    call nc_check(nf90_put_att(ncid,pint_id,'units','Pa'),                     'att')

    call nc_check(nf90_def_var(ncid,'zint',nf90_double,[pverp_dim],zint_id),'def_var:zint')
    call nc_check(nf90_put_att(ncid,zint_id,'long_name','interface heights'),'att')
    call nc_check(nf90_put_att(ncid,zint_id,'units','m'),                   'att')

    ! ---- Gas mixing ratios (pver) ----
    call nc_check(nf90_def_var(ncid,'h2ommr', nf90_double,[pver_dim],h2o_id), 'def_var:h2ommr')
    call nc_check(nf90_put_att(ncid,h2o_id,'long_name','H2O specific humidity'),  'att')
    call nc_check(nf90_put_att(ncid,h2o_id,'units','kg/kg_moist_air'),            'att')

    call nc_check(nf90_def_var(ncid,'co2mmr', nf90_double,[pver_dim],co2_id), 'def_var:co2mmr')
    call nc_check(nf90_put_att(ncid,co2_id,'long_name','CO2 mass mixing ratio'),  'att')
    call nc_check(nf90_put_att(ncid,co2_id,'units','kg/kg_dry_air'),             'att')

    call nc_check(nf90_def_var(ncid,'ch4mmr', nf90_double,[pver_dim],ch4_id), 'def_var:ch4mmr')
    call nc_check(nf90_put_att(ncid,ch4_id,'long_name','CH4 mass mixing ratio'),  'att')
    call nc_check(nf90_put_att(ncid,ch4_id,'units','kg/kg_dry_air'),             'att')

    call nc_check(nf90_def_var(ncid,'c2h6mmr',nf90_double,[pver_dim],c2h6_id),'def_var:c2h6mmr')
    call nc_check(nf90_put_att(ncid,c2h6_id,'long_name','C2H6 mass mixing ratio'), 'att')
    call nc_check(nf90_put_att(ncid,c2h6_id,'units','kg/kg_dry_air'),              'att')

    call nc_check(nf90_def_var(ncid,'o2mmr',  nf90_double,[pver_dim],o2_id),  'def_var:o2mmr')
    call nc_check(nf90_put_att(ncid,o2_id,'long_name','O2 mass mixing ratio'),    'att')
    call nc_check(nf90_put_att(ncid,o2_id,'units','kg/kg_dry_air'),              'att')

    call nc_check(nf90_def_var(ncid,'o3mmr',  nf90_double,[pver_dim],o3_id),  'def_var:o3mmr')
    call nc_check(nf90_put_att(ncid,o3_id,'long_name','O3 mass mixing ratio'),    'att')
    call nc_check(nf90_put_att(ncid,o3_id,'units','kg/kg_dry_air'),              'att')

    call nc_check(nf90_def_var(ncid,'n2mmr',  nf90_double,[pver_dim],n2_id),  'def_var:n2mmr')
    call nc_check(nf90_put_att(ncid,n2_id,'long_name','N2 mass mixing ratio'),    'att')
    call nc_check(nf90_put_att(ncid,n2_id,'units','kg/kg_dry_air'),              'att')

    call nc_check(nf90_def_var(ncid,'h2mmr',  nf90_double,[pver_dim],h2_id),  'def_var:h2mmr')
    call nc_check(nf90_put_att(ncid,h2_id,'long_name','H2 mass mixing ratio'),    'att')
    call nc_check(nf90_put_att(ncid,h2_id,'units','kg/kg_dry_air'),              'att')

    ! ---- Atmospheric constants ----
    call nc_check(nf90_def_var(ncid,'mw',nf90_double,[one_dim],mw_id),'def_var:mw')
    call nc_check(nf90_put_att(ncid,mw_id,'long_name','molecular weight of dry air'),'att')
    call nc_check(nf90_put_att(ncid,mw_id,'units','AMU'),                            'att')

    call nc_check(nf90_def_var(ncid,'cp',nf90_double,[one_dim],cp_id),'def_var:cp')
    call nc_check(nf90_put_att(ncid,cp_id,'long_name','specific heat of dry air'),'att')
    call nc_check(nf90_put_att(ncid,cp_id,'units','J/kg/K'),                      'att')

    call nc_check(nf90_enddef(ncid), 'enddef')

    ! ---- Write data ----
    call put_1d(ncid, lwup_id, LWUP,   'put:LWUP')
    call put_1d(ncid, lwdn_id, LWDN,   'put:LWDN')
    call put_1d(ncid, swup_id, SWUP,   'put:SWUP')
    call put_1d(ncid, swdn_id, SWDN,   'put:SWDN')
    call put_1d(ncid, lwhr_id, LWHR,   'put:LWHR')
    call put_1d(ncid, swhr_id, SWHR,   'put:SWHR')

    call put_scalar(ncid, ts_id, ts,   'put:ts')
    call put_scalar(ncid, ps_id, ps,   'put:ps')

    call put_1d(ncid, tmid_id, tmid,   'put:tmid')
    call put_1d(ncid, tint_id, tint,   'put:tint')
    call put_1d(ncid, pmid_id, pmid,   'put:pmid')
    call put_1d(ncid, pint_id, pint,   'put:pint')
    call put_1d(ncid, zint_id, zint,   'put:zint')

    call put_1d(ncid, h2o_id,  h2ommr, 'put:h2ommr')
    call put_1d(ncid, co2_id,  co2mmr, 'put:co2mmr')
    call put_1d(ncid, ch4_id,  ch4mmr, 'put:ch4mmr')
    call put_1d(ncid, c2h6_id, c2h6mmr,'put:c2h6mmr')
    call put_1d(ncid, o2_id,   o2mmr,  'put:o2mmr')
    call put_1d(ncid, o3_id,   o3mmr,  'put:o3mmr')
    call put_1d(ncid, n2_id,   n2mmr,  'put:n2mmr')
    call put_1d(ncid, h2_id,   h2mmr,  'put:h2mmr')

    call put_scalar(ncid, mw_id, mwdry_col, 'put:mw')
    call put_scalar(ncid, cp_id, cpdry_col, 'put:cp')

    call nc_check(nf90_close(ncid), 'close:'//filename)
    write(*,'(2a)') 'exocol_io: wrote ', trim(filename)

  end subroutine write_output

end module exocol_io
