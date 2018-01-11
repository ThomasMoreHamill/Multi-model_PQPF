SUBROUTINE  control_ensemble_quantile_mapping_x25 (nxa, nya, npct, &
    nstride, nmembers, nmemcdf, n25, thresh, conusmask, &
    precip_anal_cdf, ensemble_cdf, ensemble_ccpa, &
    ensemble_ccpa_x25, ensmean, stddev, POP)
    
! controls the quantile mapping of ensemble forecast to the analysis
! using cdfs of each.    

INTEGER, INTENT(IN) :: nxa, nya ! grid dimensions 
INTEGER, INTENT(IN) :: npct ! number of precip thresholds where CDF tallied
INTEGER, INTENT(IN) :: nstride ! stride length when skipping grid pts
INTEGER, INTENT(IN) :: nmembers ! ensemble size
INTEGER, INTENT(IN) :: n25
REAL, INTENT(IN), DIMENSION(npct) :: thresh ! threshold amts for CDF
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask ! conus mask
REAL, INTENT(IN), DIMENSION(nxa,nya,npct):: precip_anal_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya,nmemcdf,npct):: ensemble_cdf
    ! note that Canadian system has different forecast biases assoc'd with diff mbrs
    ! hence nmemcdf = nmembers for Canadian, 1 for other systems
REAL, INTENT(IN), DIMENSION(nxa,nya,nmembers) :: ensemble_ccpa
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: ensmean, stddev, POP


REAL, DIMENSION(nxa,nya) :: ensemble

REAL, INTENT(OUT), DIMENSION(n25,nxa,nya,nmembers)  :: ensemble_ccpa_x25

REAL, DIMENSION(n25,nxa,nya) :: forecast_x25 ! work array
REAL, DIMENSION(nxa,nya,npct):: work_cdf

PRINT *,'Empirical quantile mapping and dressing of ensemble, # members = ',nmembers
imemtot = 1
DO imem = 1, nmembers
    
    ! ---- with canadian ensemble, quantile map each member with a forecast
    !      cdf appropriate to that member, else use same forecast cdf for all.
    
    IF (nmemcdf .eq. nmembers) THEN
        imemcdf = imem
    ELSE
        imemcdf = 1
    ENDIF
	
    work_cdf(:,:,:) = ensemble_cdf(:,:,imemcdf,:)
    
    CALL cdf_correct_x25_local(nxa, nya, npct, nstride, n25, &
        thresh, conusmask, work_cdf, precip_anal_cdf, ensemble_ccpa(1,1,imem), &
        forecast_x25)
    ensemble(:,:) = ensemble_ccpa(:,:,imem)
    PRINT *,'imem, max(ensemble_ccpa), max(forecast_x25)', &
        maxval(ensemble), maxval(forecast_x25)
    ensemble_ccpa_x25(:,:,:,imem) = forecast_x25(:,:,:)
	imemtot = imemtot + 1
END DO

CALL calc_mean_spread(n25, nxa, nya, nmembers, &
    ensemble_ccpa_x25, conusmask, ensmean, stddev, POP)

RETURN
END SUBROUTINE control_ensemble_quantile_mapping_x25
