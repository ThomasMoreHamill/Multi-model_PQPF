SUBROUTINE raw_ensemble_probs_singlemodel(nxa, nya, nens, nthreshes, pthreshes, &
    ensemble_ccpa, prob_forecast, ensemble_mean)

INTEGER, INTENT(IN) :: nxa, nya, nens, nthreshes

REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes  ! event threshold amount
REAL, INTENT(IN), DIMENSION(nxa,nya,nens) :: ensemble_ccpa
REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: prob_forecast
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ensemble_mean

REAL, DIMENSION(nens) :: ensemble

!PRINT *,'inside raw_ensemble_probs_singlemodel'
!PRINT *,'pthreshes = ', pthreshes
emax = MAXVAL(ensemble_ccpa)
!print *,'emax = ', emax

IF (emax .le. 0.0) THEN
    PRINT *,'bad data for this day.  Setting all outpust to missing.'
    prob_forecast(:,:,:) = -99.99
ELSE
    DO ithresh = 1, nthreshes 
        rthresh = pthreshes(ithresh)
        DO jya = 1, nya
            DO ixa = 1, nxa
                ensemble(:) = ensemble_ccpa(ixa,jya,:)
                sume = 0.0
                ncount = 0
                ndenom = 0
                DO imem = 1, nens
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. rthresh) ncount = ncount + 1
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. -98.) ndenom = ndenom + 1
                    IF (ensemble_ccpa(ixa,jya,imem) .ge. -98. .and. ithresh .eq. 1) &
                        sume = sume + ensemble_ccpa(ixa,jya,imem) 
                END DO
                IF (ndenom .gt. 0.0) THEN
                    prob_forecast(ixa,jya, ithresh) = REAL(ncount) / REAL(ndenom)
                    IF (ithresh .eq. 1) ensemble_mean(ixa,jya) = sume / REAL(ndenom)
                ELSE
                    prob_forecast(ixa,jya, ithresh) = -99.99
                    IF (ithresh .eq. 1) ensemble_mean(ixa,jya) = -99.99 
                ENDIF
                
                !IF (ixa .eq. 3*nxa/4 .and. jya .eq. nya/2) &
                !    PRINT *,'ithresh rthresh, prob = ', &
                !        ithresh, rthresh, prob_forecast(ixa,jya,ithresh)
                    
            END DO  ! ixa
        END DO  ! jya
        write (6,fmt='(A40,f8.4)') 'Raw probabilities for thresh = ',rthresh
        write (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
        write (6,fmt='(1X,3(F10.5,1X))')minval(prob_forecast(:,:,ithresh)), &
            maxval(prob_forecast(:,:,ithresh)), &
            sum(prob_forecast(:,:,ithresh))/REAL(nxa*nya)
    END DO   ! ithresh

    
    !PRINT *,'RAW POP min, max, mean =', minval(prob_forecast(:,:,1)), &
    !    maxval(prob_forecast(:,:,1)), sum(prob_forecast(:,:,1))/(nya*nxa) 
    !PRINT *, 'prob_forecast(1:nxa:10,nya/2,1) = ', prob_forecast(1:nxa:10,nya/2,1)
    !PRINT *, 'prob_forecast(1:nxa:10,nya/2,3) = ', prob_forecast(1:nxa:10,nya/2,3)
    
ENDIF

RETURN
END SUBROUTINE raw_ensemble_probs_singlemodel
