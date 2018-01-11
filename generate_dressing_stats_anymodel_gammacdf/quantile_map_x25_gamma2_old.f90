
SUBROUTINE quantile_map_x25_gamma2(nxa, nya, nstride, n25, nens_qmap, &
    conusmask, gamma_shape_qmap_forecast, &
    gamma_scale_qmap_forecast, fraction_zero_qmap_forecast, &
    gamma_shape_qmap_analysis, gamma_scale_qmap_analysis, &
    fraction_zero_qmap_analysis, forecast, forecast_x25)

! --- Perform the quantile mapping; here we use 25 (24 + original) 
!     surrounding grid points forecasts as well to increase sample  
!     size and account for position error

INTEGER, INTENT(IN) :: nxa, nya, nstride, n25, nens_qmap
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_shape_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_scale_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: fraction_zero_qmap_forecast
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_shape_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: gamma_scale_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa, nya) :: fraction_zero_qmap_analysis
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast

REAL, INTENT(OUT), DIMENSION(n25,nxa,nya) :: forecast_x25

REAL :: weight, weight2, rmean
DOUBLE PRECISION alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION alpha_anal, beta_anal, fz_anal
DOUBLE PRECISION fcst
DOUBLE PRECISION precip_qmapped
DOUBLE PRECISION qgamma
INTEGER :: idum,itest,jtest, jyn, ixn

rmean = 1000.*SUM(forecast)/REAL(nxa*nya)
idum = -1*INT(rmean)
forecast_x25(:,:,:) = -99.99 ! set to missing

! ---- Not for points inside conus mask, apply the procedure of doing quantile
!      mapping using the forecast at (i,j), but also surrounding locations.
!      also add random number to the input forecast quantile to also account for the
!      tendency of the ensemble to be over-certain of its amount.

weight = 0.0
DO jya = 1, nya
	DO ixa = 1, nxa
        
		IF (conusmask(ixa,jya) .le. 0) THEN
         	! If outside the conus, then make the output forecast array simply replicates
         	! of the input forecast array
         	forecast_x25(:,ixa,jya) = forecast(ixa,jya)
      	ELSE
        
         	ktr = 0

         	! ---- Loop thru the 25 grid points with (ixa,jya) in center

         	DO jyn = jya-2*nstride, jya+2*nstride, nstride
            	DO ixn = ixa-2*nstride, ixa+2*nstride, nstride
					
                    ! ---- deal with the possibility that the stencil point is outside 
                    !      domain.  In that case, use the value on the nearest border point
                    
                    IF (ixn .lt. 1) THEN
                        ixs = 1
                    ELSE IF (ixn .gt. nxa) THEN
                        ixs = nxa
                    ELSE
                        ixs = ixn
                    ENDIF
                    IF (jyn .lt. 1) THEN
                        jys = 1
                    ELSE IF (jyn .gt. nya) THEN
                        jys = nya
                    ELSE
                        jys = jyn
                    ENDIF
                    
                    ktr = ktr + 1
                    IF (conusmask(ixs,jys) .gt. 0) THEN

                        fcst = forecast(ixs,jys)
                        
						IF (forecast(ixs,jys) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = 0.0
                        ELSE

                     		! ---- Perform quantile mapping.  First determine the quantile of 
                            !      forecast distribution associated with today's forecast value.
                            
                            alpha_fcst = gamma_shape_qmap_forecast(ixs,jys)
                            beta_fcst = gamma_scale_qmap_forecast(ixs,jys)
                            fz_fcst = fraction_zero_qmap_forecast(ixs,jys)
                            alpha_anal = gamma_shape_qmap_analysis(ixa,jya)
                            beta_anal = gamma_scale_qmap_analysis(ixa,jya)
                            fz_anal = fraction_zero_qmap_analysis(ixa,jya)
                            CALL gamma_quantile_map(fcst, alpha_fcst, &
                                beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                fz_anal, precip_qmapped)
                            forecast_x25(ktr,ixa,jya) = precip_qmapped 
						ENDIF
					
                    ELSE 
                        
                        ! ---- stencil point is either outside the area where quantile mapping possible.
                        !      As a substitute, take the quantile-mapped value at center point and
                        !      add a small amount of random noise.
                        
                        fcst = forecast(ixa,jya) 
						IF (forecast(ixa,jya) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = 0.0
                        ELSE
                        
                            alpha_fcst = gamma_shape_qmap_forecast(ixa,jya)
                            beta_fcst = gamma_scale_qmap_forecast(ixa,jya)
                            fz_fcst = fraction_zero_qmap_forecast(ixa,jya)
                            alpha_anal = gamma_shape_qmap_analysis(ixa,jya)
                            beta_anal = gamma_scale_qmap_analysis(ixa,jya)
                            fz_anal = fraction_zero_qmap_analysis(ixa,jya)
                            CALL gamma_quantile_map(fcst, alpha_fcst, &
                                beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                fz_anal, precip_qmapped)
                            forecast_x25(ktr,ixa,jya) = precip_qmapped
                            r = 1 + 0.1*(ran1(idum) - 0.5)
                            forecast_x25(ktr,ixa,jya) = forecast_x25(ktr,ixa,jya)*r 
                        ENDIF      
                        
                 	END IF ! jyn>1, jyn < nya, etc.
                    
            	END DO ! ixn
         	END DO ! jyn
      	END IF ! conusmask
   	END DO ! ixa
END DO ! jya

RETURN
END SUBROUTINE quantile_map_x25_gamma2

! ==================================================================

SUBROUTINE gamma_quantile_map (fcst, alpha_fcst, beta_fcst, fz_fcst, &
    alpha_anal, beta_anal, fz_anal, qmapped_fcst)
    
DOUBLE PRECISION, INTENT(IN) :: fcst, alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION, INTENT(IN) :: alpha_anal, beta_anal, fz_anal 
DOUBLE PRECISION, INTENT(OUT) :: qmapped_fcst 
    
DOUBLE PRECISION ksi, cum, ccum, nonexceedance_prob, qgamma
LOGICAL tootrue
LOGICAL toofalse
    
! ---- determine the cumulative probability of being less than 
!      or equal to the forecast amount using the CDF defined by
!      the forecast gamma distribution parameters and the forecast
!      fraction zero  

tootrue = .TRUE.
toofalse = .FALSE.

ksi = fcst / beta_fcst
CALL cumgam(ksi, alpha_fcst, cum, ccum)
nonexceedance_prob = fz_fcst + (1.-fz_fcst)*cum

! ---- Next, using Michael Scheuerer's code, use the quantile 
!      function to retrieve the analyzed value associates with
!      the same cumulative probability.

IF (nonexceedance_prob .gt. fz_anal) THEN
    cum = (nonexceedance_prob - fz_anal) / (1.-fz_anal)
    qmapped_fcst = qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse)
    IF (qmapped_fcst .lt. 0.0) qmapped_fcst = 0.0

ELSE ! map to zero
    qmapped_fcst = 0.0

ENDIF

RETURN
END SUBROUTINE gamma_quantile_map
    

