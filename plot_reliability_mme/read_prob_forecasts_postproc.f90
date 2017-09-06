! ============================================================================

SUBROUTINE read_prob_forecasts_postproc(nxa, nya, cyyyymmddhh, cmodel, clead_use, &
    rthresh, nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
    prob_forecast, climo_prob)

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: cyyyymmddhh
CHARACTER*(*), INTENT(IN) :: cmodel 
CHARACTER*(*), INTENT(IN) :: clead_use
REAL, INTENT(IN) :: rthresh
REAL, INTENT(OUT), DIMENSION(nthreshes) :: pthreshes

CHARACTER*120 infile
CHARACTER*25 cfield
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlonsa, rlatsa, &
    prob_forecast, climo_prob
LOGICAL iex

infile = '/Users/thamill/precip/ecmwf_data/'//TRIM(cmodel)//'_'//&
    TRIM(clead_use)//'h_IC'//cyyyymmddhh//'.nc'
PRINT *,TRIM(infile)
INQUIRE (file=infile,exist=iex)
IF (iex) THEN
    
    ! ---- open the file

    netid = 0
    CALL check (nf90_open(infile,NF90_NOWRITE,netid))

    ! ---- read in the list of dates/times in yyyymmddhh format associated with
    !      each time index stored in the netcdf file

    cfield ='pthreshes'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,pthreshes,&
    	start=(/1/),count=(/nthreshes/)))
        
    DO ithresh = 1, nthreshes
        IF (ABS(pthreshes(ithresh) - rthresh) .lt. 0.01) GOTO 1000
    END DO
    PRINT *,'Did not find threshold = ', rthresh,' in subroutine verify_relia_bss_anymodel.f90'
    PRINT *,'Stopping.'
    STOP
    1000 CONTINUE    
    PRINT *,'using thresh = ',ithresh,rthresh
        
    cfield ='conusmask'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,conusmask,&
        start=(/1,1/),count=(/nxa,nya/)))
    
    cfield ='rlonsa'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,rlonsa,&
        start=(/1,1/),count=(/nxa,nya/)))

    cfield ='rlatsa'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,rlatsa,&
        start=(/1,1/),count=(/nxa,nya/)))
    	
    cfield ='prob_forecast'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,prob_forecast,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))
    	
    cfield ='climo_prob'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,climo_prob,&
        start=(/1,1,ithresh/),count=(/nxa,nya,1/)))

    CALL check(nf90_close(netid))
ENDIF
RETURN
END SUBROUTINE read_prob_forecasts_postproc
