
SUBROUTINE read_precip_climatology_local(nxa, nya, pclimo_infile, cthresh, climo_prob, &
    rlonsa, rlatsa, conusmask)

! ---- read in the climatological probability of analyzed precipitation.   Previously 
!      calculated in compute_climatology_ppn_multithresh.py

USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: pclimo_infile
CHARACTER*(*), INTENT(IN) :: cthresh
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: climo_prob, rlonsa, rlatsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask

CHARACTER*20 cfield

! ---- Initialize

conusmask(:,:)=0
climo_prob(:,:) = 0.
rlonsa(:,:) = 0.
rlatsa(:,:) = 0.

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

IF (TRIM(cthresh) .eq. 'POP') THEN
    cfield = 'climo_prob_POP'
ELSE IF (TRIM(cthresh) .eq. '1mm') THEN
    cfield = 'climo_prob_1mm'
ELSE IF (TRIM(cthresh) .eq. '2p5mm') THEN
    cfield = 'climo_prob_2p5mm'
ELSE IF (TRIM(cthresh) .eq. '5mm') THEN
    cfield = 'climo_prob_5mm'
ELSE IF (TRIM(cthresh) .eq. '10mm') THEN
    cfield = 'climo_prob_10mm'
ELSE IF (TRIM(cthresh) .eq. '25mm') THEN
    cfield = 'climo_prob_25mm'   
ELSE IF (TRIM(cthresh) .eq. '50mm') THEN
    cfield = 'climo_prob_50mm'
ELSE
    write(6,*)' **** Invalid threshold value: ',cthresh
    stop
ENDIF

CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check(nf90_get_var(netid, ivar, climo_prob,&
    start=(/1,1/),count=(/nxa,nya/)))

! ---- Close netcdf file.

CALL check(nf90_close(netid))

print *,'done reading read_precip_climatology_local'

RETURN
END subroutine read_precip_climatology_local
