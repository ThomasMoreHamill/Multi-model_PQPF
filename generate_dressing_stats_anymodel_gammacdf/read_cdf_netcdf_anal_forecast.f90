SUBROUTINE read_cdf_netcdf_anal_forecast(nxa, nya, npct, &
    nmembers, iyyyymmddhh, cleade, cmodel, data_directory, &
    precip_anal_cdf, ensemble_cdf, thresh)
    
    
    
!    read_cdf_netcdf_anal_forecast(nxa, nya, npct, &
!        nens_cdf, iyyyymmddhh, cleade, cprefix, data_directory, &
!        precip_anal_cdf, ensemble_cdf, thresh)
  
! note trickery here.  CMC ensemble doesn't have members with
! exchangeable errors, so we need to keep track of CDFs member
! by member.  To make code generalizable to other systems 
! whose members have exchangeable errors (NCEP, ECMWF, etc)
! we have the third array dimension of ensemble_cdf be
! set to 1 (on input) for these ensemble systems. 
    
USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, npct, nmembers, iyyyymmddhh
CHARACTER*(*), INTENT(IN) :: cleade, cmodel, data_directory
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: precip_anal_cdf
REAL, INTENT(OUT), DIMENSION(nxa,nya,nmembers, npct) :: ensemble_cdf
REAL, INTENT(OUT), DIMENSION(npct) :: thresh

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_cdf_3d
CHARACTER*3 clead
CHARACTER*256 infilename
CHARACTER*10 cyyyymmddhh_start, cyyyymmddhh_end 
CHARACTER*20 cfield

IF (nmembers .eq. 1) ALLOCATE (ensemble_cdf_3d(nxa,nya,npct))


! ---- determine which CDF file to use based on the current date

clead = cleade
IF (cleade .eq. '12' .or. cleade .eq. '012') clead = '12'
IF (cleade .eq. '24' .or. cleade .eq. '024') clead = '24'
IF (cleade .eq. '36' .or. cleade .eq. '036') clead = '36'
IF (cleade .eq. '48' .or. cleade .eq. '048') clead = '48'
IF (cleade .eq. '60' .or. cleade .eq. '060') clead = '60'
IF (cleade .eq. '72' .or. cleade .eq. '072') clead = '72'
IF (cleade .eq. '84' .or. cleade .eq. '084') clead = '84'
IF (cleade .eq. '96' .or. cleade .eq. '096') clead = '96'

! ---- Eric, this is my kludge.  When I developed my own test environment,
!      I didn't want to keep updating a CDF every day.  So instead I 
!      created a CDF file for every 10th day.  The code below decides
!      which of these to read in.   You replace with your own CDF files.

IF (iyyyymmddhh .lt. 2016013000) THEN
    iyyyymmddhh_reference = 2016013000
ELSE IF (iyyyymmddhh .ge. 2016013100 .and. iyyyymmddhh .lt. 2016021000) THEN
    iyyyymmddhh_reference = 2016021000
ELSE IF (iyyyymmddhh .ge. 2016021100 .and. iyyyymmddhh .lt. 2016022000) THEN
    iyyyymmddhh_reference = 2016022000
ELSE IF (iyyyymmddhh .ge. 2016022100 .and. iyyyymmddhh .lt. 2016022800) THEN
    iyyyymmddhh_reference = 2016022800
ELSE IF (iyyyymmddhh .ge. 2016022900 .and. iyyyymmddhh .lt. 201603100) THEN
    iyyyymmddhh_reference = 201603100
ELSE IF (iyyyymmddhh .ge. 2016031100 .and. iyyyymmddhh .lt. 2016032000) THEN
    iyyyymmddhh_reference = 2016032000
ELSE IF (iyyyymmddhh .ge. 2016032100 .and. iyyyymmddhh .lt. 2016033000) THEN
    iyyyymmddhh_reference = 2016033000                        
ELSE IF (iyyyymmddhh .ge. 2016033100 .and. iyyyymmddhh .lt. 2016041000) THEN
    iyyyymmddhh_reference = 2016041000
ELSE IF (iyyyymmddhh .ge. 2016042000 .and. iyyyymmddhh .lt. 2016043000) THEN
    iyyyymmddhh_reference = 2016042000   
ELSE IF (iyyyymmddhh .ge. 2016043000 .and. iyyyymmddhh .lt. 2016051000) THEN
    iyyyymmddhh_reference = 2016043000 
ELSE IF (iyyyymmddhh .ge. 2016051000 .and. iyyyymmddhh .lt. 2016052000) THEN
    iyyyymmddhh_reference = 2016051000 
ELSE IF (iyyyymmddhh .ge. 2016052000 .and. iyyyymmddhh .lt. 2016053000) THEN
    iyyyymmddhh_reference = 2016052000 
ELSE IF (iyyyymmddhh .ge. 2016053000 .and. iyyyymmddhh .lt. 2016061000) THEN
    iyyyymmddhh_reference = 2016053000 
ELSE IF (iyyyymmddhh .ge. 2016061000 .and. iyyyymmddhh .lt. 2016062000) THEN
    iyyyymmddhh_reference = 2016061000 
ELSE 
    iyyyymmddhh_reference = 2016062000 
ENDIF

READ (cleade,'(i3)') ileade
idaysbefore = - 1 -int(ileade/24)
ihoursbefore = idaysbefore*24
call updat(iyyyymmddhh_reference,ihoursbefore,iyyyymmddhh_end)
call updat(iyyyymmddhh_end,-24*61,iyyyymmddhh_start)
PRINT *,'expecting file dates from ',iyyyymmddhh_start,' to ',iyyyymmddhh_end 
WRITE (cyyyymmddhh_start,'(i10)') iyyyymmddhh_start
WRITE (cyyyymmddhh_end,'(i10)') iyyyymmddhh_end 


infilename = TRIM(data_directory)// TRIM(cmodel) // '_CDF_flead' // &
    TRIM(clead) // '_' // cyyyymmddhh_start // '_to_' // &
    cyyyymmddhh_end // '.nc'
PRINT *, 'reading cdf information from ', TRIM(infilename)

print *,'npct ',npct

! --- read cdf and precip analysis cdf

CALL check (nf90_open(infilename,NF90_NOWRITE,netid))

cfield='panal_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,precip_anal_cdf,&
    start=(/1,1,1/),count=(/nxa,nya,npct/)))
    
cfield='pfcst_CDF'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
IF (TRIM(cmodel) .ne. 'CMC') THEN
    CALL check(nf90_get_var(netid,ivar, ensemble_cdf_3D,&
        start=(/1,1,1/),count=(/nxa,nya,npct/)))
    ensemble_cdf(:,:,1,:) = ensemble_cdf_3d(:,:,:)
ELSE ! now for CMC
    CALL check(nf90_get_var(netid,ivar, ensemble_cdf,&
        start=(/1,1,1,1/),count=(/nxa,nya,nmembers,npct/)))
ENDIF

! ---- read in the precipitation amounts associated with 
!      CDF indices.

cfield='thrval'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar, thresh,&
    start=(/1/),count=(/npct/)))

CALL check(nf90_close(netid))

IF (nmembers .eq. 1) DEALLOCATE(ensemble_cdf_3d)

RETURN
END SUBROUTINE  read_cdf_netcdf_anal_forecast
