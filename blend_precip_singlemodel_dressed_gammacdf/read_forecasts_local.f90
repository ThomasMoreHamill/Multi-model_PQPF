SUBROUTINE read_forecasts_local (nxa, nya, nens, &
	infile_early, infile_late, ensemble_ccpa)

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya, nens
CHARACTER*(*), INTENT(IN) :: infile_early, infile_late
REAL, INTENT(OUT), DIMENSION(nxa,nya,nens) :: ensemble_ccpa

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ensemble_ccpa_early, ensemble_ccpa_late

CHARACTER*20 cfield

! ---- Initialize

ensemble_ccpa(:,:,:) = -99.99
ALLOCATE(ensemble_ccpa_early(nxa,nya,nens), ensemble_ccpa_late(nxa,nya,nens))

! ---- Open the early file, read in

netid=0
PRINT *,'netid, reading from ',netid, TRIM(infile_early)
CALL check (nf90_open(infile_early,NF90_NOWRITE,netid))

cfield='apcp_fcst_ens'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,ensemble_ccpa_early,&
           start=(/1,1,1/),count=(/nxa,nya,nens/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))


! ---- Open the late file, read in

netid=0
PRINT *,'netid, reading from ',netid, TRIM(infile_late)
CALL check (nf90_open(infile_late,NF90_NOWRITE,netid))

cfield='apcp_fcst_ens'
CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid,ivar,ensemble_ccpa_late,&
           start=(/1,1,1/),count=(/nxa,nya,nens/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))
!PRINT *,'max of earlier forecast = ', maxval(ensemble_ccpa_early)
!PRINT *,'max of later forecast   = ', maxval(ensemble_ccpa_late)

IF (ensemble_ccpa_early(1,1,1) .gt. -98. .and. &
    ensemble_ccpa_late(1,1,1) .gt. -98.) THEN
    ensemble_ccpa(:,:,:) = ensemble_ccpa_late(:,:,:) - ensemble_ccpa_early(:,:,:) 
ELSE
    ensemble_ccpa(:,:,:) = ensemble_ccpa_late(:,:,:) - ensemble_ccpa_early(:,:,:) 
ENDIF

! ---- Reset slightly negative forecasts to zero

DO iens = 1, nens
	DO jya = 1, nya
		DO ixa = 1, nxa
			IF (ensemble_ccpa(ixa,jya,iens) .gt. -2.0 .and. &
			ensemble_ccpa(ixa,jya,iens) .lt. 0) ensemble_ccpa(ixa,jya,iens) = 0.0
		END DO
	END DO
END DO

pmax = maxval(ensemble_ccpa)
IF (pmax .eq. 0.0) ensemble_ccpa(:,:,:) = -99.99

!print *,'ensemble_ccpa(1:nxa:10,nya/2,1) = ', ensemble_ccpa(1:nxa:10,nya/2,1)
!print *,'max of ensemble_ccpa = ', maxval(ensemble_ccpa)
!print *,'done reading in subroutine read_forecasts_local'

DEALLOCATE(ensemble_ccpa_early, ensemble_ccpa_late)

RETURN
END subroutine read_forecasts_local
