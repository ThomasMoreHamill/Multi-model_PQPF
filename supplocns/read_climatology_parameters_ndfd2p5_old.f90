
SUBROUTINE read_climatology_parameters_ndfd2p5(nxa, nya, infile, &
    fraction_zero, alphahat, betahat, conusmask)    
    
USE netcdf
INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: infile
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: fraction_zero, alphahat, betahat
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya)  :: conusmask

CALL check(nf90_open(infile, NF90_NOWRITE, netid))

CALL check(nf90_inq_varid(netid, "fraction_zero", ivar))
CALL check(nf90_get_var(netid, ivar, fraction_zero, &
    start=(/1,1/), count=(/nxa,nya/)))
    
CALL check(nf90_inq_varid(netid, "alpha", ivar))
CALL check(nf90_get_var(netid, ivar, alphahat, &
    start=(/1,1/), count=(/nxa,nya/)))

CALL check(nf90_inq_varid(netid, "beta", ivar))
CALL check(nf90_get_var(netid, ivar, betahat, &
    start=(/1,1/), count=(/nxa, nya/)))
    
CALL check(nf90_inq_varid(netid, "conusmask", ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa,nya/)))

CALL check(nf90_close(netid))

RETURN
END SUBROUTINE read_climatology_parameters_ndfd2p5