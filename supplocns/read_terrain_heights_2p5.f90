SUBROUTINE read_terrain_heights_2p5(nxa, nya, infile, &
    terrain, lonsa, latsa)
    
USE netcdf
INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: infile
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: terrain, lonsa, latsa

! ---- read in data 

CALL check(nf90_open(infile, NF90_NOWRITE, netid))

CALL check(nf90_inq_varid(netid, "longitude", ivar))
CALL check(nf90_get_var(netid, ivar, lonsa, &
    start=(/1,1/), count=(/nxa,nya/)))
    
CALL check(nf90_inq_varid(netid, "latitude", ivar))
CALL check(nf90_get_var(netid, ivar, latsa, &
    start=(/1,1/), count=(/nxa,nya/)))

CALL check(nf90_inq_varid(netid, "terrain", ivar))
CALL check(nf90_get_var(netid, ivar, terrain, &
    start=(/1,1/), count=(/nxa, nya/)))
    
PRINT *, 'max, min terrain height = ', maxval(terrain), &
    minval(terrain)
        
CALL check(nf90_close(netid))

!  ---- make sure that the data is as expected; if not,
!       change or stop program execution.

rmaxlon = MAXVAL(lonsa)
print *,'lonsa(1:nxa:10, nya/2)', lonsa(1:nxa:10, nya/2)

IF (rmaxlon .gt. 180.) lonsa = lonsa - 360.
print *,'min, max lonsa = ', minval(lonsa), maxval(lonsa)

rminlat = MINVAL(latsa)
rmaxlat = MAXVAL(latsa)
IF (rminlat .lt. 0.0 .or. rmaxlat .gt. 90.0) THEN
    PRINT *, 'Problem detected in read_terrain_heights_2p5: '
    PRINT *, 'Unexpected latitudes read in from ', TRIM(infile)
    PRINT *, 'Latitudes range from ',latsa(1,1),' to ',latsa(1,nya)
    PRINT *, 'Stopping.'
    STOP
ENDIF

! ---- Death Valley is -86 m, and highest peak in CONUS is 4421m (Mt Whitney)
!      though there might be higher peaks in S. Canada.  Check terrain for
!      realistic values and exit if bounds exceeded.

tmin = MINVAL(terrain)
tmax = MAXVAL(terrain)
IF (tmin .lt. -100 .or. tmax .gt. 5500.) THEN
    PRINT *,'Problem detected in read_terrain_heights_2p5: '
    PRINT *,'Unrealistically small or large terrain height values found.'
    PRINT *,'Min, max terrain height = ', tmin, tmax
    PRINT *,'Stopping.'
    STOP
ENDIF

        

RETURN
END SUBROUTINE read_terrain_heights_2p5




