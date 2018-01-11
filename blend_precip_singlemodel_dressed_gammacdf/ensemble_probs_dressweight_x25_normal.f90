SUBROUTINE ensemble_probs_dressweight_x25_normal (n25, nxa, nya, nens, &
    nmembersx25, npcatvals, nthreshes, pthreshes, &
    ensemble_ccpa_x25, closest_histogram, precip_histogram_thresholds, &
    conusmask, climo_prob, prob_forecast_qmapped, prob_forecast)  

    ! ---- a simplified version of the dressing routine, whereby we 
    !      use the variable weighting of members from closest-member 
    !      histograms, but not the fancified dressing procedure.
    !      Instead we have normal pdfs that follow the procedure from the
    !      previous Hamill et al. (2017) article, with normal dressing
    !      distributions with the amplitude of the dressing noise 
    !      linearly related to the ensemble forecast member's value.
    
INTEGER, INTENT(IN) :: n25, nxa, nya, nens, nmembersx25
INTEGER, INTENT(IN) :: npcatvals
INTEGER, INTENT(IN) :: nthreshes

REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes
REAL, INTENT(IN), DIMENSION(n25,nxa,nya,nens) :: ensemble_ccpa_x25
REAL, INTENT(IN), DIMENSION(nmembersx25,npcatvals) :: closest_histogram
REAL, INTENT(IN), DIMENSION(npcatvals) :: precip_histogram_thresholds

INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: climo_prob
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: prob_forecast_qmapped
REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: prob_forecast

REAL, DIMENSION(n25, nens) :: ensemble
REAL, DIMENSION(nmembersx25) :: ensemble_sorted
REAL*8, DIMENSION(nmembersx25) :: weights, weight_lower, weight_upper
REAL*8 sumprobability, zval, denom, rthresh, cum, ccum

INTEGER, DIMENSION(1) :: ishape

ishape(1) = nmembersx25
prob_forecast(:,:,:) = -99.99
iprint = 1
rmean = 1000.*SUM(ensemble_ccpa_x25)/REAL(n25*nxa*nya*nens)
idum = -1*INT(rmean)
PRINT *,'npcatvals = ', npcatvals
PRINT *,'precip_histogram_thresholds = ', precip_histogram_thresholds

! ---- process all grid points inside mask
    
DO jya = 1, nya
    DO ixa = 1, nxa
    !DO jya = nya/2, nya/2
    !DO ixa = 3*nxa/4, 3*nxa/4
        IF (conusmask(ixa,jya) .eq. 1) THEN
            
            ! ---- sort the ensemble
            
            ensemble(:,:) = ensemble_ccpa_x25(:,ixa,jya,:)       
            ensemble_sorted = RESHAPE(ensemble,ishape)
            emean = SUM(ensemble) / REAL(n25*nens)
            !PRINT *,'ixa, jya, emean = ',ixa, jya, emean
            CALL sort(nmembersx25, ensemble_sorted) 
            
            ! ---- determine based on the ensemble mean which closest
            !      histogram to use.
            
            CALL determine_category(emean, npcatvals, &
                precip_histogram_thresholds, ipcat)
            IF (ipcat .gt. 1) THEN 
                weights(:) = closest_histogram(:,ipcat-1) 
            ELSE
                weights(:) = 1.0 / DBLE(nmembersx25)
            ENDIF
            !PRINT *,'weights(1:10),weights(nmembersx25-10:nmembersx25) = ',&
            !    weights(1:10),weights(nmembersx25-10:nmembersx25)
            
            ! ---- determine exceedance probability for each threshold 
            !      from a weighted linear combination of dressed member
            !      pdfs, here the dressing distributions normal with
            !      magnitude linearly related to the member (quantile-mapped) 
            !      precipitation amount.
            
            DO ithresh = 1, nthreshes
            !DO ithresh = 1, 1
                rthresh = pthreshes(ithresh)
                sumprobability = 0.
                DO imem = 1, nmembersx25
                    IF (ensemble_sorted(imem) .gt. 0.0) THEN
                        denom = 0.15 + 0.15*ensemble_sorted(imem) ! std deviation for this mbr
                        zval = (rthresh - ensemble_sorted(imem)) / denom ! convert to std normal deviate
                        CALL cumnor(zval, cum, ccum)  ! ccum is the exceedance probability
                    ELSE
                        ccum = 0.0
                    ENDIF
                    sumprobability = sumprobability + ccum*DBLE(weights(imem))
                    !PRINT *,'imem, ens, denom, rthresh, zval, ccum, sump = ', &
                    !    imem, ensemble_sorted(imem), denom, &
                    !    rthresh, zval, ccum, sumprobability
                END DO ! imem
                IF (sumprobability .gt. 0.9999) sumprobability = 1.0
                IF (sumprobability .lt. 0.0) sumprobability = 0.0
                prob_forecast(ixa,jya,ithresh) = sumprobability
            END DO
        ENDIF ! conusmask
    END DO ! ixa
END DO ! jya
    
! ---- some diagnostic print statements.

DO ithresh = 1, nthreshes    
    WRITE (6,fmt='(A40,f8.4)') 'Dressed probabilities for thresh = ',pthreshes(ithresh)
    WRITE (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
    WRITE (6,fmt='(1X,3(F10.5,1X))')&
        MINVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        MAXVAL(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        SUM(conusmask(:,:)*prob_forecast(:,:,ithresh))/REAL(nxa*nya)    
END DO ! ithresh = 1, nthreshes

RETURN
END SUBROUTINE ensemble_probs_dressweight_x25_normal