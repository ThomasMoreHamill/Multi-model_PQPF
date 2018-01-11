SUBROUTINE load_precipquantiles_ccpa_analonly_v9 (imonth, nxa, nya, npct, &
     nxa_small, nya_small, ioffset, joffset, conusmask, lonsa, latsa, pctv, &
     panal_quantiles, terrain)
!
! purpose: read in the data set for this month of the precipitation quantiles, 
!    both forecast and analyzed, as well as the rank correlation between 
!    forecast and analysis.
! 
USE netcdf

INTEGER, INTENT(IN) :: imonth, nxa, nya, npct, nxa_small, nya_small
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: panal_quantiles
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: lonsa, latsa
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: terrain
REAL, INTENT(OUT), DIMENSION(npct) :: pctv

CHARACTER*80 infile
CHARACTER*20 cfield
CHARACTER*3, DIMENSION(12) :: cmonths

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: panal_quantiles_in

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun',&
	          'Jul','Aug','Sep','Oct','Nov','Dec'/

ALLOCATE (panal_quantiles_in(nxa_small,nya_small,npct))

panal_quantiles(:,:,:) = -99.99 ! initialize expanded array to missing.

! ---- open the file (464,224) dimension

infile = '/Projects/Reforecast2/netcdf/' // &
	'refcstv2_apcp_CCPAgrid_CDF_fhour024_to_048_'// cmonths(imonth)//'.nc'

PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

cfield = 'pct'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,pctv,start=(/1/),count=(/npct/)))

cfield = 'panal_quantiles'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,panal_quantiles_in, &
	start=(/1,1,1/),count=(/nxa_small,nya_small,npct/)))

! --- close file.

CALL check(nf90_close(netid))

! ---- Now map data from the small CCPA domain to the expanded one.
!      from Finan, the information on grid offset:
!      Here are the grid specs for the expanded domain:
!      LL: 22.313N, 232.437E
!      nxlat = 515
!      nylon = 262
!      offset: (1,1) = (21,22) < set in routine above to ioffset, joffset

panal_quantiles(1+ioffset:ioffset+nxa_small,1+joffset:joffset+nya_small,:) = &
     panal_quantiles_in(:,:,:)

! ---- open the terrain file on the expanded 1/8-degree grid;  
!      read in the terrain height and practical_mask and conusmask

infile = 'terrain_expanded_eighthdegree_from_NDFD.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

cfield ='lons'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lonsa,&
  start=(/1,1/),count=(/nxa,nya/)))

cfield ='lats'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,latsa,&
  start=(/1,1/),count=(/nxa,nya/)))

! ---- read in the conus mask (1 if inside conus) , 
!      and lats and lons 

cfield ='conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,conusmask,&
  start=(/1,1/),count=(/nxa,nya/)))

terrain(:,:) = 0.0
cfield ='terrain_height'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,terrain,&
  start=(/1,1/),count=(/nxa,nya/)))

CALL check(nf90_close(netid))

DO ixa = 1, nxa
	DO jya = 1, nya
		IF(conusmask(ixa,jya) .lt. 0) conusmask(ixa,jya) = 0
	END DO
END DO

DEALLOCATE (panal_quantiles_in)

RETURN
END subroutine load_precipquantiles_ccpa_analonly_v9
