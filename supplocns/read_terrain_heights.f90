SUBROUTINE read_terrain_heights(nxa, nya, infile, terrain, lonsa, latsa, &
    conusmask)
    
USE netcdf
INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: infile
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: terrain, lonsa, latsa
INTEGER*2, DIMENSION(nxa,nya), INTENT(OUT) :: conusmask


CALL check(nf90_open(infile, NF90_NOWRITE, netid))

CALL check(nf90_inq_varid(netid, "lons", ivar))
CALL check(nf90_get_var(netid, ivar, lonsa, &
    start=(/1,1/), count=(/nxa,nya/)))
    
CALL check(nf90_inq_varid(netid, "lats", ivar))
CALL check(nf90_get_var(netid, ivar, latsa, &
    start=(/1,1/), count=(/nxa,nya/)))

CALL check(nf90_inq_varid(netid, "conusmask", ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa, nya/)))

CALL check(nf90_inq_varid(netid, "terrain_height", ivar))
CALL check(nf90_get_var(netid, ivar, terrain, &
    start=(/1,1/), count=(/nxa, nya/)))

CALL check(nf90_close(netid))

RETURN
END SUBROUTINE read_terrain_heights




