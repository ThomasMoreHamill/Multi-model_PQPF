SUBROUTINE read_supp_locns_ccpa(cmonth, nx_ccpa, ny_ccpa, &
    nsupp, xlocation_ccpa, ylocation_ccpa, lons_ccpa, lats_ccpa, &
    conusmask_ccpa)

USE netcdf

CHARACTER (len=3), INTENT(IN) :: cmonth

INTEGER, INTENT(IN) :: nx_ccpa, ny_ccpa, nsupp

INTEGER, INTENT(OUT), DIMENSION(nx_ccpa,ny_ccpa,nsupp) :: xlocation_ccpa, ylocation_ccpa
INTEGER(2), INTENT(OUT), DIMENSION(nx_ccpa,ny_ccpa) :: conusmask_ccpa
REAL, INTENT(OUT), DIMENSION(nx_ccpa,ny_ccpa) :: lons_ccpa, lats_ccpa

CHARACTER (len=120) :: infile
CHARACTER (len=20) :: cfield

!infile = '/Users/thamill/refcst2/test/ppn_analog_locns_ccpa9_'//cmonth//'.dat'
!PRINT *, 'reading from ', TRIM(infile)
!OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='unformatted')
!READ (1) nxfin, nyfin, nxain, nyain, nsuppin
!READ (1) xlocation_ccpa
!READ (1) ylocation_ccpa
!READ (1) conusmask_ccpa
!READ (1) lons_ccpa
!READ (1) lats_ccpa
!CLOSE (1)


print *,'nx_ccpa, ny_ccpa, nsupp = ', nx_ccpa, ny_ccpa, nsupp
infile = 'supplemental_locations_eighth_degree_' // cmonth // '_v9.nc'
CALL check (nf90_open(infile,NF90_NOWRITE,netid))
PRINT *,'netid, reading from ',netid, TRIM(infile)

cfield = 'longitudes'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lons_ccpa, &
    start=(/1,1/),count=(/nx_ccpa,ny_ccpa/)))
print *,'done reading longitudes'
    
cfield = 'latitudes'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lats_ccpa, &
    start=(/1,1/),count=(/nx_ccpa,ny_ccpa/)))
print *,'done reading latitudes'

cfield = 'xlocations'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,xlocation_ccpa, &
    start=(/1,1,1/),count=(/nx_ccpa,ny_ccpa,nsupp/)))
print *,'done reading xlocations'
    
cfield = 'ylocations'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,ylocation_ccpa, &
    start=(/1,1,1/),count=(/nx_ccpa,ny_ccpa,nsupp/)))
print *,'done reading ylocations'

cfield = 'conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
print *,'ivar = ', ivar
CALL check (nf90_get_var(netid,ivar,conusmask_ccpa, &
    start=(/1,1/),count=(/nx_ccpa,ny_ccpa/)))
print *,'done reading conusmask'

! --- close file.

CALL check(nf90_close(netid))


RETURN
END SUBROUTINE read_supp_locns_ccpa