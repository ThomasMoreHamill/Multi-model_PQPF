SUBROUTINE find_climo_category (nccats, climo_pop_thresholds, climo_prob, iclim)
	
INTEGER, INTENT(IN) :: nccats
REAL, INTENT(IN), DIMENSION(nccats) :: climo_pop_thresholds	
REAL, INTENT(IN) :: climo_prob
INTEGER, INTENT(OUT) :: iclim

DO icat = 1, nccats+1
	IF (icat .eq. 1) THEN
		rlow = 0.0
		rhigh = climo_pop_thresholds(1)
	ELSE IF (icat .eq. nccats+1) THEN
		rlow = climo_pop_thresholds(nccats)
		rhigh = 1.0
	ELSE 
		rlow = climo_pop_thresholds(icat-1)
		rhigh = climo_pop_thresholds(icat)
	ENDIF
	
	IF (climo_prob .ge. rlow .and. climo_prob .lt. rhigh) THEN
		iclim = icat
		GOTO 6000
	ENDIF
END DO

6000 RETURN
END SUBROUTINE find_climo_category
	