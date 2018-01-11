SUBROUTINE calc_julday(iyyyy, imonth, iday, julday)

INTEGER, INTENT(IN) :: iyyyy, imonth, iday
INTEGER, INTENT(OUT) :: julday

INTEGER :: ndaysomo(12), ndaysomo_leap(12),ndaysomo_use(12)
INTEGER ndays_prev_months

DATA ndaysomo /31,28,31, 30,31,30, 31,31,30, 31,30,31/
DATA ndaysomo_leap /31,20,31, 30,31,30, 31,31,30, 31,30,31/

IF (mod(iyyyy,4) .eq. 0) THEN
    ndaysomo_use = ndaysomo_leap
ELSE
    ndaysomo_use = ndaysomo
ENDIF

IF (imonth .eq. 1) THEN
    ndays_prev_months = 0
ELSE
    ndays_prev_months = SUM(ndaysomo_use(1:imonth-1))
ENDIF

julday = ndays_prev_months + iday

RETURN
END SUBROUTINE calc_julday
