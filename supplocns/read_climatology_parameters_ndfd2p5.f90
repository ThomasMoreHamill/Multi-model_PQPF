
SUBROUTINE read_climatology_parameters_ndfd2p5(nxa, nya, npthresh, &
    infile, conusmask, practical_mask, lonsa, latsa, fraction_zero, &
    cdf, pthreshes)    
    
USE netcdf
INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: infile

REAL, INTENT(IN), DIMENSION(nxa,nya) :: lonsa, latsa
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: fraction_zero
REAL, INTENT(OUT), DIMENSION(nxa,nya,npthresh) :: cdf
REAL, INTENT(OUT), DIMENSION(npthresh) :: pthreshes
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: practical_mask

! ---- read in the data

CALL check(nf90_open(infile, NF90_NOWRITE, netid))

CALL check(nf90_inq_varid(netid, "fraction_zero", ivar))
CALL check(nf90_get_var(netid, ivar, fraction_zero, &
    start=(/1,1/), count=(/nxa,nya/)))
    
CALL check(nf90_inq_varid(netid, "quantile_fpamt", ivar))
CALL check(nf90_get_var(netid, ivar, cdf, &
    start=(/1,1,1/), count=(/nxa,nya,npthresh/)))
    
CALL check(nf90_inq_varid(netid, "pthreshes", ivar))
CALL check(nf90_get_var(netid, ivar, pthreshes, &
    start=(/1/), count=(/npthresh/)))

CALL check(nf90_inq_varid(netid, "validmask", ivar))
CALL check(nf90_get_var(netid, ivar, practical_mask, &
    start=(/1,1/), count=(/nxa, nya/)))
        
CALL check(nf90_inq_varid(netid, "conusmask", ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa, nya/)))

CALL check(nf90_close(netid))

! ---- check for realistic values.

fmin = MINVAL(fraction_zero*conusmask)
fmax = MAXVAL(fraction_zero*conusmask)
IF (fmin .lt. 0.0 .or. fmax .gt. 1.0) THEN
    PRINT *, 'Problem detected in read_climatology_parameters_ndfd2p5: '
    PRINT *, 'Min/max fraction zero over conus not bounded by 0/1'
    PRINT *, 'Min, max = ', fmin, fmax
    DO ixa = 1, nxa
        DO jya = 1, nya
            IF (conusmask(ixa,jya) .eq. 1 .and. fraction_zero(ixa,jya) .lt. 0.0) THEN
                PRINT *,'ixa, jya, lon, lat, fraction_zero = ', &
                ixa, jya, lonsa(ixa,jya), latsa(ixa,jya), fraction_zero(ixa,jya)
            ENDIF
        END DO
    END DO
    PRINT *, 'Stopping.'
    STOP
ENDIF

DO ithresh = 1, npthresh
    qmin = MINVAL(cdf(:,:,ithresh)*conusmask(:,:))
    qmax = MAXVAL(cdf(:,:,ithresh)*conusmask(:,:))
    IF (qmin .lt. 0.0 .or. qmax .gt. 1.0) THEN
        PRINT *, 'Problem detected in read_climatology_parameters_ndfd2p5: '
        PRINT *, 'Min/max cdf over conus not bounded by 0/1 for ithresh = ', ithresh
        PRINT *, 'Min, max = ', qmin, qmax
        PRINT *, 'Stopping.'
        STOP
    ENDIF
END DO

pmin = MINVAL(pthreshes)
pmax = MAXVAL(pthreshes)
IF (pmin .lt. 0.0 .or. pmax .gt. 200.0) THEN
    PRINT *,'Problem detected in read_climatology_parameters_ndfd2p5: '
    PRINT *,'Unexpected value in pthreshes array:'
    PRINT *,'pthreshes = ', phthreshes
    PRINT *,'Stopping.'
    STOP
ENDIF

imaskmin = MINVAL(conusmask)
imaskmax = MAXVAL(conusmask)
IF (imaskmin .ne. 0 .or. imaskmax .ne. 1) THEN
    PRINT *,'Problem detected in read_climatology_parameters_ndfd2p5: '
    PRINT *,'Conusmask should be 1/0'
    PRINT *,'Minval, maxval of conusmask is ', imaskmin, imaskmax
    PRINT *,'Stopping.'
    STOP
ENDIF

imaskmin = MINVAL(practical_mask)
imaskmax = MAXVAL(practical_mask)
IF (imaskmin .ne. 0 .or. imaskmax .ne. 1) THEN
    PRINT *,'Problem detected in read_climatology_parameters_ndfd2p5: '
    PRINT *,'practical_mask should be 1/0'
    PRINT *,'Minval, maxval of practical_mask is ', imaskmin, imaskmax
    PRINT *,'Stopping.'
    STOP
ENDIF


RETURN
END SUBROUTINE read_climatology_parameters_ndfd2p5