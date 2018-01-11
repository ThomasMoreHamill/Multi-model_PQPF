
SUBROUTINE read_precip_climatology_local(nxa, nya, pclimo_infile, climo_prob, &
    rlonsa, rlatsa, conusmask)

! ---- read in the climatological probability of analyzed precipitation.   Previously 
!      calculated in compute_climatology_ppn_multithresh.py

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: pclimo_infile
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: climo_prob, rlonsa, rlatsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
INTEGER nthreshes
REAL, ALLOCATABLE, DIMENSION(:) :: pthreshes

CHARACTER*20 cfield

! ---- Initialize

conusmask(:,:) = 0
climo_prob(:,:) = 0.
rlonsa(:,:) = 0.
rlatsa(:,:) = 0.

rthresh = 0.254


! ---- Open the file, use the number of times in the file later to allocate a date array

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
    
cfield='thra'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_inquire_dimension(netid, ivar, cfield, nthreshes)) 

ALLOCATE (pthreshes(nthreshes))
cfield = 'pthreshes'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, pthreshes,&
    start=(/1/), count=(/nthreshes/)))   

PRINT *,'pthreshes = ',pthreshes
DO ithresh = 1, nthreshes
    IF (ABS(pthreshes(ithresh) - rthresh) .lt. 0.01) THEN
        idxtouse = ithresh
        PRINT *,'idxtouse = ',idxtouse
        GOTO 2000
    ENDIF
END DO

1000 PRINT *,'Error in read_precip_climatology_local: did not find pthreshes value to match rthresh'
PRINT *,'pthreshes = ',pthreshes
PRINT *,'rthresh = ',rthresh
PRINT *,'stopping.'
STOP

2000 CONTINUE

cfield = 'climo_prob'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, climo_prob,&
    start=(/1,1,idxtouse/),count=(/nxa,nya,1/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))
PRINT *,'done reading read_precip_climatology_local'

! ---- check data values

probmin = MINVAL(climo_prob*conusmask)
probmax = MAXVAL(climo_prob*conusmask)
PRINT *,'min, max climo_prob = ', probmin, probmax

IF (probmin .lt. 0.0 .or. probmax .gt. 1.0) THEN    
    PRINT *, 'Error in read_precip_climatology_local.  '
    PRINT *, 'Probabilities out of bounds for thresh = ',cthresh,' . Stopping.'
    STOP
ENDIF

rminlon = MINVAL(rlonsa)
IF (rminlon .gt. 0.0) THEN
    PRINT *,'minimum longitude input does not have longitudes below zero.   Fixing by subtracting 360.'
    rlonsa = rlonsa-360.
ENDIF

rminlat = MINVAL(rlatsa)
rmaxlat = MAXVAL(rlatsa)
IF (rminlat .lt. -90.0 .or. rmaxlat .gt. 90.0) THEN
    PRINT *, 'Error in read_precip_climatology_local.  Stopping. '
    PRINT *, 'Latitudes out of bounds; min, max lat = ', rminlat, rmaxlat
    STOP
ENDIF

DEALLOCATE(pthreshes)

RETURN
END subroutine read_precip_climatology_local
