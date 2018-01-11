SUBROUTINE read_precipitation_analysis(nxa,nya,&
	iyyyymmddhh_end, infile, analysis, istat)

! read in the two 6-hourly precipitation analyses,
! sum, and return that after making sure that
! missing values are set to -99.99

USE netcdf
	
INTEGER, INTENT(IN) :: nxa, nya, iyyyymmddhh_end
CHARACTER*(*) :: infile
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: analysis

INTEGER, INTENT(OUT) :: istat

CHARACTER*20 cfield
	
INTEGER ntimes
INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_anal
REAL, DIMENSION(nxa, nya) :: analysis_early, analysis_late

! ---- open the file

netid = 0
CALL check (nf90_open(infile, NF90_NOWRITE, netid))

! ---- read in the list of dates/times in yyyymmddhh format associated with
!      each time index stored in the netcdf file

!print *,'time'
cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_inquire_dimension(netid,ivar,cfield,ntimes))
ALLOCATE (iyyyymmddhh_anal(ntimes))

!print *,'yyyymmddhh_anal_end'
cfield ='yyyymmddhh_anal_end'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,iyyyymmddhh_anal,&
	start=(/1/),count=(/ntimes/)))

! ---- determine the time slice we want to read in
	
ifound = -99
DO itime = 1, ntimes
	IF (iyyyymmddhh_anal(itime) .eq. iyyyymmddhh_end) THEN
		ifound = 1
		GOTO 6000
	ENDIF
END DO
6000 CONTINUE  

IF (ifound .eq. 1) THEN
	cfield = 'apcp_anal'
    !print *,'apcp_anal'
    CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
	PRINT *,'trying to read analysis, records ',&
        itime, itime-1,' nxa, nya = ', nxa, nya
    CALL check (nf90_get_var(netid,ivar, analysis_late,&
        start=(/1,1,itime/),count=(/nxa,nya,1/)))
    CALL check (nf90_get_var(netid,ivar, analysis_early,&
        start=(/1,1,itime-1/),count=(/nxa,nya,1/)))
    analysis = analysis_late + analysis_early
    DO jya = 1, nya
        DO ixa = 1, nxa
            IF (analysis(ixa,jya) .lt. -99.99) analysis(ixa,jya) = -99.99
        END DO
    END DO   
	istat = 1
ELSE
	analysis(:,:) = -99.99
	istat = -1
ENDIF

! --- close netcdf file.

CALL check(nf90_close(netid))

DEALLOCATE (iyyyymmddhh_anal)

RETURN
END SUBROUTINE read_precipitation_analysis