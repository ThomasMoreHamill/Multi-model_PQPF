SUBROUTINE ensemble_probs_HXLR (nclim_p1_vals, nclim_vals, nthreshes, nxa, nya, &
    pthreshes, climo_pop_thresholds, gamma_shape_fclimpop, gamma_scale_fclimpop, &
    fraction_zeros_fclimpop, b0_mean, b1_mean, b0_spread, &
    b1_spread, ensmean_pxform, stddev_pxform, conusmask, &
    climo_prob, prob_forecast_HXLR)
    
INTEGER, INTENT(IN) :: nclim_p1_vals, nthreshes, nxa, nya
REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: gamma_shape_fclimpop
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: gamma_scale_fclimpop
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: fraction_zeros_fclimpop
REAL, INTENT(IN), DIMENSION(nclim_vals) :: climo_pop_thresholds
REAL, INTENT(IN) :: b0_mean, b1_mean, b0_spread, b1_spread
REAL, INTENT(IN), DIMENSION(nxa, nya)  :: ensmean_pxform, stddev_pxform
INTEGER*2, INTENT(IN), DIMENSION(nxa, nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: climo_prob

REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: prob_forecast_HXLR


REAL*8 alpha, beta, ksi, cum, ccum

pzero_xform = 0.01**0.33

PRINT *,'b0_mean, b1_mean, b0_spread, b1_spread = ', b0_mean, b1_mean, b0_spread, b1_spread
PRINT *,'max, min ensmean_pxform = ', maxval(ensmean_pxform), minval(ensmean_pxform)
PRINT *,'max, min stddev_pxform = ', maxval(stddev_pxform), minval(stddev_pxform)
PRINT *,'fraction_zeros_fclimpop = ', fraction_zeros_fclimpop
PRINT *,'climo_pop_thresholds = ', climo_pop_thresholds

DO ithresh = 1, nthreshes

    rthresh = pthreshes(ithresh)
    rthresh_0p33 = rthresh**0.33
    !rthresh_0p33 = rthresh
    PRINT *,'processing threshold = ', rthresh

    ! ---- process all grid points inside mask

    DO jya = 1, nya
        DO ixa = 1, nxa
            
            IF (conusmask(ixa,jya) .eq. 1) THEN      
                          
                emean = ensmean_pxform(ixa,jya)
                espread = stddev_pxform(ixa,jya)
                
                IF (emean .lt. pzero_xform) THEN

                    ! ---- keep original dressing functionality and estimate
                    !      postprocessed values from past statistics of
                    !      analyzed precipitation given zero mean forecast
                    !      and given the climatological probability of precip.
                    
                    CALL determine_category(climo_prob(ixa,jya,1), &
                        nclim_vals, climo_pop_thresholds, iclim)
                    !PRINT *, 'climo_prob(ixa,jya,1), iclim = ',climo_prob(ixa,jya,1), iclim
                    alpha = gamma_shape_fclimpop(iclim)
                    beta = gamma_scale_fclimpop(iclim)
                    fcz = fraction_zeros_fclimpop(iclim)
                    ksi = rthresh / beta
                    IF (alpha .ne. alpha) THEN ! happens only when NaN
                        sumprobability = 0.
                    ELSE
                        CALL cumgam(ksi, alpha, cum, ccum)
                        sumprobability = (1.-fcz)*ccum
                    ENDIF
                    
                    !IF (jya .eq. nya/2 .and. ithresh .eq. 1) THEN
                    !    PRINT *,'ixa, iclim, clim  alpha beta fz sumprob ',ixa, iclim, &
                    !        climo_prob(ixa,jya,1), alpha, beta, fraction_zeros_fclimpop(iclim), sumprobability    
                    !    !PRINT *,'fraction_zeros_fclimpop, iclim = ', fraction_zeros_fclimpop, iclim
                    !ENDIF
                    
                ELSE
                    
                    ! ---- use heteroscedastic ensemble logistic regression
                    !      to estimate probabilities.
                    
                    predmean = b0_mean + b1_mean * emean
                    predspread = exp(b0_spread + b1_spread * espread)
                    z = (rthresh_0p33 - predmean) / predspread
                    expz = exp(z)
                    prob_forecast_HXLR(ixa,jya,ithresh) = 1.0 - expz / (1.0 + expz)
            
                ENDIF ! emean .lt. pzero_xform
            ENDIF ! conusmask
        END DO ! ixa
    END DO ! jya
        
    write (6,fmt='(A40,f8.4)') 'HXLR probabilities for thresh = ',rthresh
    write (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
    write (6,fmt='(1X,3(F10.5,1X))')&
        minval(conusmask(:,:)*prob_forecast_HXLR(:,:,ithresh)), &
        maxval(conusmask(:,:)*prob_forecast_HXLR(:,:,ithresh)), &
        SUM(prob_forecast_HXLR(:,:,ithresh))/REAL(nxa*nya)
    
END DO ! ithresh

RETURN
END SUBROUTINE ensemble_probs_HXLR
                    
                    