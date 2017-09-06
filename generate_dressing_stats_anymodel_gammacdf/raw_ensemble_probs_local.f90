SUBROUTINE raw_ensemble_probs_local(nxa, nya, nens_cmc, nens_ncep, nens_ecmwf, &
    rthresh, cmc_ensemble_ccpa, ncep_ensemble_ccpa, ecmwf_ensemble_ccpa, &
    prob_forecast_raw, prob_forecast_raw_CMC, prob_forecast_raw_NCEP, &
    prob_forecast_raw_ECMWF, ensemble_mean)

! --- 

INTEGER, INTENT(IN) :: nxa, nya, nens_cmc, nens_ncep, nens_ecmwf 

REAL, INTENT(IN) :: rthresh  ! event threshold amount

REAL, INTENT(IN), DIMENSION(nxa,nya,nens_cmc) :: cmc_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ncep) :: ncep_ensemble_ccpa
REAL, INTENT(IN), DIMENSION(nxa,nya,nens_ecmwf) :: ecmwf_ensemble_ccpa

REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_NCEP
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_CMC
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw_ECMWF
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: ensemble_mean
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: raw_ensemble

idum = -23434 ! seed for random number generator

! ---- determine the number of members in the downscaled multi-model ensemble, and then allocate
!      arrays to hold that data.  In order to give the deterministic/control forecasts from the 
!      various centers more weight, we'll let them be downscaled multiple times (nmult_ecmwf and
!      nmult_other) with different analysis data used for the downscaling each time.  Also allocate
!      array that will hold a list of indices of the analysis data that is closest to this particular
!      forecast (idate_list)

nens = nens_cmc + nens_ncep + nens_ecmwf
ALLOCATE (raw_ensemble(nxa,nya,nens))

! ---- copy the input deterministic and ensemble data into the work array.  Keep track
!      of how many members have good data (not -99.99)

nens_good = 1
nens_cmc_good = 0
nens_ncep_good = 0
nens_ecmwf_good = 0

cmax = MAXVAL(cmc_ensemble_ccpa)
rnmax = MAXVAL(ncep_ensemble_ccpa)
emax = MAXVAL(ecmwf_ensemble_ccpa)
IF (cmax .le. 0.0 .or. rnmax .le. 0.0 .or. emax .le. 0.0) THEN
    PRINT *,'bad data for this day.  Setting all outpust to missing.'
    prob_forecast_raw(:,:) = -99.99
    prob_forecast_raw_NCEP(:,:) = -99.99
    prob_forecast_raw_CMC(:,:) = -99.99
    prob_forecast_raw_ECMWF(:,:) = -99.99
    ensemble_mean(:,:) = -99.99
ELSE
    PRINT *,'copying in ncep data'
    DO i = 1, nens_ncep
        IF (ncep_ensemble_ccpa(nxa/2,nya/2,i) .ge. -98.) THEN
            raw_ensemble(:,:,nens_good) = ncep_ensemble_ccpa(:,:,i)
            nens_ncep_good = nens_ncep_good + 1
            nens_good = nens_good + 1
        ENDIF
    END DO

    PRINT *,'copying in cmc data'
    DO i = 1, nens_cmc
        IF (cmc_ensemble_ccpa(nxa/2,nya/2,i) .ge. -98.) THEN
            raw_ensemble(:,:,nens_good) = cmc_ensemble_ccpa(:,:,i)
            nens_cmc_good = nens_cmc_good + 1
            nens_good = nens_good + 1
        ENDIF
    END DO

    PRINT *,'copying in ecmwf data'
    DO i = 1, nens_ecmwf
        IF (ecmwf_ensemble_ccpa(nxa/2,nya/2,i) .ge. -98.) THEN
            raw_ensemble(:,:,nens_good) = ecmwf_ensemble_ccpa(:,:,i)
            nens_ecmwf_good = nens_ecmwf_good + 1
            nens_good = nens_good + 1
        ENDIF
    END DO
    nens_good = nens_good - 1

    PRINT *,nens_good,' valid ensemble members out of ',nens
    PRINT *, 'rthresh = ', rthresh

    ! ---- now get the probabilities from ensemble relative frequency

    ! ---- NCEP raw probability

    IF (nens_ncep_good .gt. 0) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                ncount_ncep = 0
                ndenom = 0
                DO imem = 1, nens_ncep
                    IF (ncep_ensemble_ccpa(ixa,jya,imem) .ge. rthresh) &
                        ncount_ncep = ncount_ncep + 1
                    IF (ncep_ensemble_ccpa(ixa,jya,imem) .ge. -98.) &
                        ndenom = ndenom + 1
                END DO
                prob_forecast_raw_NCEP(ixa,jya) = REAL(ncount_ncep) / REAL(ndenom)
            END DO
        END DO
    ELSE
        prob_forecast_raw_NCEP(:,:) = -99.99
    ENDIF

    ! ---- CMC raw probability

    IF (nens_cmc_good .gt. 0) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                ncount_cmc = 0
                ndenom = 0
                DO imem = 1, nens_cmc
                    IF (cmc_ensemble_ccpa(ixa,jya,imem) .ge. rthresh) &
                        ncount_cmc = ncount_cmc + 1
                    IF (cmc_ensemble_ccpa(ixa,jya,imem) .ge. -98.) &
                        ndenom = ndenom + 1
                END DO
                prob_forecast_raw_CMC(ixa,jya) = REAL(ncount_cmc) / REAL(ndenom)
            END DO
        END DO
    ELSE
        prob_forecast_raw_CMC(:,:) = -99.99
    ENDIF

    ! ---- ECMWF raw probability

    IF (nens_ecmwf_good .gt. 0) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                ncount_ecmwf = 0
                ndenom = 0
                DO imem = 1, nens_ecmwf
                    IF (ecmwf_ensemble_ccpa(ixa,jya,imem) .ge. rthresh) &
                        ncount_ecmwf = ncount_ecmwf + 1
                    IF (ecmwf_ensemble_ccpa(ixa,jya,imem) .ge. -98.) &
                        ndenom = ndenom + 1
                END DO
                prob_forecast_raw_ECMWF(ixa,jya) = REAL(ncount_ecmwf) / REAL(ndenom)
            END DO
        END DO
    ELSE
        prob_forecast_raw_ECMWF(:,:) = -99.99
    ENDIF

    ! ---- MME raw probability

    IF (nens_good .gt. 0) THEN
        DO jya = 1, nya
            DO ixa = 1, nxa
                ncount_mme = 0
                ndenom = 0
                DO imem = 1, nens_good
                    IF (raw_ensemble(ixa,jya,imem) .ge. rthresh) &
                        ncount_mme = ncount_mme + 1
                    IF (raw_ensemble(ixa,jya,imem) .ge. -98) &
                        ndenom = ndenom + 1
                END DO
                prob_forecast_raw(ixa,jya) = REAL(ncount_mme) / REAL(ndenom)
                IF (ixa .eq. nxa/2 .and. jya .eq. nya/2) PRINT *,'ncount_mme, ndenom = ', ncount_mme, ndenom
            END DO
        END DO
    ELSE
        prob_forecast_raw(:,:) = -99.99
    ENDIF

    ensemble_mean = SUM(raw_ensemble(:,:,1:nens_good), 3) / REAL(nens_good)

    PRINT *, 'number of good NCEP, CMC, ECMWF, MME members = ', & 
        nens_ncep_good, nens_cmc_good, nens_ecmwf_good, nens_good

    PRINT *,'RAW min, max, mean =', minval(prob_forecast_raw), &
        maxval(prob_forecast_raw), sum(prob_forecast_raw)/(nya*nxa)
    PRINT *,'NCEP min, max, mean =', minval(prob_forecast_raw_NCEP), &
        maxval(prob_forecast_raw_NCEP),sum(prob_forecast_raw_NCEP)/(nya*nxa)
    PRINT *,'CMC min, max, mean =', minval(prob_forecast_raw_CMC), &
        maxval(prob_forecast_raw_CMC),sum(prob_forecast_raw_CMC)/(nya*nxa)
    PRINT *,'ECMWF min, max, mean =', minval(prob_forecast_raw_ECMWF), &
        maxval(prob_forecast_raw_ECMWF),sum(prob_forecast_raw_ECMWF)/(nya*nxa)
        
ENDIF

! ---- deallocate the arrays

DEALLOCATE (raw_ensemble)

RETURN
END SUBROUTINE raw_ensemble_probs_local
