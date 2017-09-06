SUBROUTINE read_closest_histogram_singlemodel(nmembers, npcatvals, &
    histofile, closest_histogram, precip_histogram_thresholds)

USE netcdf
INTEGER, INTENT(IN) :: nmembers, npcatvals
CHARACTER*(*), INTENT(IN) :: histofile
REAL, INTENT(OUT), DIMENSION(nmembers,npcatvals) :: closest_histogram
REAL, INTENT(OUT), DIMENSION(npcatvals) :: precip_histogram_thresholds

PRINT *,'nmembers, npcatvals  = ', nmembers, npcatvals
PRINT *,TRIM(histofile)

CALL check(nf90_open(TRIM(histofile), NF90_NOWRITE, netid))
CALL check(nf90_inq_varid(netid, "closest_histogram", ivar))
CALL check(nf90_get_var(netid, ivar, closest_histogram, &
    start=(/1,1/), count=(/nmembers,npcatvals/)))
CALL check(nf90_inq_varid(netid, "precip_thresholds", ivar))    
CALL check(nf90_get_var(netid, ivar, precip_histogram_thresholds, &
    start=(/1/), count=(/npcatvals/)))    
CALL check(nf90_close(netid))

RETURN

END SUBROUTINE read_closest_histogram_singlemodel
	
