
SUBROUTINE calculate_quantile_for_precip_values(npthreshes, nya, nxa, &
    rthreshes, fraction_zero, alphahat, &
    betahat, conusmask_in, quant_fpamt)

! compile with f2py -c -m calculate_quantile_for_precip_values calculate_quantile_for_precip_values.f90 cumgam.f90 qnorm.f90 pnorm.f90 dnorm.f90 gamma.f90 pgamma.f90 dgamma.f90 qgamma.f90 gamma_inc.f90 error_f.f90 error_fc.f90 rlog.f90 rexp.f90 exparg.f90 gam1.f90 ipmpar.f90 cumnor.f90 r8_swap.f90

INTEGER, INTENT(IN) :: npthreshes, nxa, nya 
REAL*8, INTENT(IN), DIMENSION(nya, nxa) :: alphahat, betahat, fraction_zero
REAL, DIMENSION(npthreshes) :: rthreshes
INTEGER*2, INTENT(IN), DIMENSION(nya, nxa) :: conusmask_in
REAL, INTENT(OUT), DIMENSION(npthreshes, nya, nxa) :: quant_fpamt

REAL*8 ksi, cum, ccum, alpha, beta

! f2py intent(in) npthreshes, nya, nxa
! f2py intent(in) rthreshes
! f2py intent(in) conusmask_in, fraction_zero, alphahat, betahat
! f2py intent(out) quant_fpamt
! f2py depend(npthreshes) rthreshes
! f2py depend(nya,nxa) fraction_zero, alphahat, betahat, conusmask_in
! f2py depend(npthreshes,nya,nxa) quant_fpamt

print *,'calculate_quantile_for_precip_values'
print *,'npthreshes, nya, nxa = ', npthreshes, nya, nxa 

DO ixa = 1, nxa
    DO jya = 1, nya
        alpha = alphahat(jya,ixa)
        beta = betahat(jya,ixa)
        IF (conusmask_in(jya,ixa) .eq. 1 .and. alpha .gt.0) THEN
            DO ithresh = 1, nthreshes
                IF (fraction_zero(jya,ixa) .lt. 1.0) THEN
                    ksi = rthreshes(ithresh) / beta
                    CALL cumgam(ksi, alpha, cum, ccum) ! ccum is exceedance probability
                    quant_fpamt(ithresh,jya,ixa) = fraction_zero(jya,ixa) + &
                        (1.0 - fraction_zero(jya,ixa))*cum
                ELSE
                    quant_fpamt(ithresh,jya,ixa) = 1.0 
                ENDIF   
            END DO                            
        ELSE
            quant_fpamt(:,jya,ixa) = -99.99 
        ENDIF
    END DO ! jya
END DO ! ixa   

RETURN
END SUBROUTINE calculate_quantile_for_precip_values
