SUBROUTINE read_precipitation_analysis(nxa,nya,&
	iyyyymmddhh_end, infile, analysis, istat)

USE netcdf
	
INTEGER, INTENT(IN) :: nxa, nya, iyyyymmddhh_end
CHARACTER*(*) :: infile
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: analysis

INTEGER, INTENT(OUT) :: istat

CHARACTER*12 cfield
	
INTEGER ntimes
INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_anal

! ---- open the file

netid = 0
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

! ---- read in the list of dates/times in yyyymmddhh format associated with
!      each time index stored in the netcdf file

cfield ='time'
CALL check(nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_inquire_dimension(netid,ivar,cfield,ntimes))
ALLOCATE (iyyyymmddhh_anal(ntimes))

cfield ='yyyymmddhh'
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
    CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
	PRINT *,'trying to read analysis, record ',itime, ' nxa, nya = ', nxa, nya
    CALL check (nf90_get_var(netid,ivar, analysis,&
         start=(/1,1,itime/),count=(/nxa,nya,1/)))
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