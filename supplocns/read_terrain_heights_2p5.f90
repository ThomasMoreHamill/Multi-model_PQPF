SUBROUTINE read_terrain_heights_2p5(nxa, nya, infile, &
    terrain, lonsa, latsa, conusmask)
    
USE netcdf
INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*(*), INTENT(IN) :: infile
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: terrain, lonsa, latsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask

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
    
CALL check(nf90_inq_varid(netid, "conusmask", ivar))
CALL check(nf90_get_var(netid, ivar, conusmask, &
    start=(/1,1/), count=(/nxa, nya/)))

PRINT *, 'max, min terrain height = ', maxval(terrain), &
    minval(terrain)
    
CALL check(nf90_close(netid))

RETURN
END SUBROUTINE read_terrain_heights_2p5




