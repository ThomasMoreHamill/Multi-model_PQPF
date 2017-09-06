SUBROUTINE determine_precip_value_category(npvals, precip_values, &
    ensemble_member, ipvlower)
    
    ! category number closest to but lower than
    
INTEGER, INTENT(IN) :: npvals
REAL, INTENT(IN), DIMENSION(npvals) :: precip_values
REAL, INTENT(IN) :: ensemble_member
INTEGER, INTENT(OUT) :: ipvlower

ipvlower = -99
DO ival = 1, npvals-1
    IF (ensemble_member .ge. precip_values(ival) .and. &
    ensemble_member .lt. precip_values(ival+1)) THEN
        ipvlower = ival
        GOTO 7000
    ENDIF
END DO

PRINT *,'in subroutine determine_precip_value_category, could not find '
PRINT *,'suitable value of lower category.  Ensemble member value = '
PRINT *,ensemble_member,' and precip_values array = ', precip_values
PRINT *,'Stopping'
STOP

7000 RETURN                           
                            
END SUBROUTINE  determine_precip_value_category   