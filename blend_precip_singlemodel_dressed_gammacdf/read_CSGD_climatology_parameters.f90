SUBROUTINE read_CSGD_climatology_parameters (infile_CSGD_climatology, &
    nxa, nya, iyear, imo, iday, CSGD_climo_mean, CSGD_climo_mu, &
    CSGD_climo_sigma, CSGD_climo_shift)
    
USE netCDF
CHARACTER*(*), INTENT(IN) :: infile_CSGD_climatology
INTEGER, INTENT(IN) :: iyear, imo, iday
INTEGER julday

REAL, DIMENSION(nxa,nya), INTENT(OUT) :: CSGD_climo_mean
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: CSGD_climo_mu
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: CSGD_climo_sigma
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: CSGD_climo_shift

LOGICAL exist
CHARACTER*28 cfield

INQUIRE (file=infile_CSGD_climatology, exist=exist)
IF (exist) THEN

    CALL calc_julday(iyyyy, imonth, iday, julday)
    IF (mod(iyear, 4) .eq. 0 .and. julday .gt. 60) julday = julday-1
    PRINT *, 'reading CSGD information from ',TRIM(infile_CSGD_climatology)

    CALL check (nf90_open(infile_CSGD_climatology,NF90_NOWRITE,netid))

    cfield='mean_cl'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,CSGD_climo_mean,&
        start=(/1,1,julday/),count=(/nxa,nya,1/)))
    
    cfield='mu_cl'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,CSGD_climo_mu,&
        start=(/1,1,julday/),count=(/nxa,nya,1/)))
        
    cfield='sigma_cl'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,CSGD_climo_sigma,&
        start=(/1,1,julday/),count=(/nxa,nya,1/)))     

    cfield='shift_cl'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,CSGD_climo_shift,&
        start=(/1,1,julday/),count=(/nxa,nya,1/)))    

    CALL check(nf90_close(netid))
    
ELSE
    
    PRINT *,'Unable to find file = ',TRIM(infile_CSGD_climatology)
    PRINT *,'Stopping in subroutine read_CSGD_climatology_parameters'
    STOP
    
ENDIF

PRINT *,'CSGD_climo_mean(1:nxa:10,nya/2) = ',CSGD_climo_mean(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_mu(1:nxa:10,nya/2) = ',CSGD_climo_mu(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_sigma(1:nxa:10,nya/2) = ',CSGD_climo_sigma(1:nxa:10,nya/2)
PRINT *,'CSGD_climo_shift(1:nxa:10,nya/2) = ',CSGD_climo_shift(1:nxa:10,nya/2)


RETURN
END SUBROUTINE read_CSGD_climatology_parameters