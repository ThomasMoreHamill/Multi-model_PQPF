
SUBROUTINE read_precip_climatology_multi_thresh(nxa, &
    nya, nthreshes, pthreshes, pclimo_infile, &
    climo_prob, rlonsa, rlatsa, conusmask)

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, nthreshes
REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes
CHARACTER*(*), INTENT(IN) :: pclimo_infile
REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: climo_prob
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask

CHARACTER*20 cfield
REAL, DIMENSION(nthreshes) :: pthresh_in

! ---- Initialize

conusmask(:,:)=0
climo_prob(:,:,:) = 0.
rlonsa(:,:) = 0.
rlatsa(:,:) = 0.

! ---- Open the file, read in the mask, lat, lon

netid=0
PRINT *,'netid, reading from ',netid, TRIM(pclimo_infile)
CALL check (nf90_open(pclimo_infile, NF90_NOWRITE, netid))

cfield='conusmask'
CALL check(nf90_inq_varid(netid, trim(adjustl(cfield)), ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa,nya/)))

cfield='lonsa'
CALL check(nf90_inq_varid(netid, trim(adjustl(cfield)), ivar))
CALL check(nf90_get_var(netid, ivar, rlonsa,&
    start=(/1,1/), count=(/nxa,nya/)))

cfield='latsa'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, rlatsa,&
    start=(/1,1/), count=(/nxa,nya/)))
    
cfield='pthreshes'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, pthresh_in,&
        start=(/1/), count=(/nthreshes/))) 

! ---- confirm that the climatology to be read in 
!      is calculated at the expected thresholds

DO i = 1, nthresh
    IF (ABS(pthreshes(i)-pthresh_in(i)) .gt. 0.01) THEN
        PRINT *,'The climatologies are not at the expected thresholds.  Stopping.'
        PRINT *,'The program expects ',pthreshes(:)
        PRINT *,'The netCDF file has ',pthresh_in(:)
        PRINT *,'Stopping.' 
        STOP
    ENDIF
END DO 
    
! ---- read in the climatological probability

cfield = 'climo_prob'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, climo_prob,&
    start=(/1,1,1/),count=(/nxa,nya,nthreshes/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))

PRINT *,'done reading climatologies in subroutine read_precip_climatology_multi_thresh'

RETURN
END subroutine read_precip_climatology_multi_thresh
