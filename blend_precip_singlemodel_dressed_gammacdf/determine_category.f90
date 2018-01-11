SUBROUTINE determine_category(rinput, nvals, cat_dividers, ipcat)

REAL, INTENT(IN) :: rinput
INTEGER, INTENT(IN) :: nvals
REAL, INTENT(IN), DIMENSION(nvals) :: cat_dividers
INTEGER, INTENT(OUT) :: ipcat   
         
!PRINT *,'nvals+1 = ',nvals+1                           
DO icat = 1, nvals+1
    IF (icat .eq. 1) THEN
        plow = 0.0
        phigh = cat_dividers(icat)
    ELSE IF (icat .eq. nvals+1) THEN
        plow = cat_dividers(icat-1)
        phigh = 9999.
    ELSE 
        plow = cat_dividers(icat-1)
        phigh = cat_dividers(icat)
    ENDIF
    !PRINT *,'icat, rinput, plow, phigh = ', icat, rinput, plow, phigh 
    IF (rinput .ge. plow .and. rinput .lt. phigh) THEN
        ipcat = icat
        GOTO 6000
    ENDIF
END DO

6000 CONTINUE
RETURN
END SUBROUTINE determine_category