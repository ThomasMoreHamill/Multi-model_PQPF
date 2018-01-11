
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

REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99

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

        !PRINT *,'processing ixa, jya, conusmask ',ixa,jya,conusmask(ixa,jya)
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
                    
                    !PRINT *,'ixs, jys = ',ixs,jys, conusmask(ixs,jys)
                    
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
                            !PRINT *,'alpha_fcst, beta_fcst, fz_fcst = ', &
                            !    alpha_fcst, beta_fcst, fz_fcst
                            !PRINT *,'alpha_anal, beta_anal, fz_anal = ', &
                            !    alpha_anal, beta_anal, fz_anal  
                            IF (alpha_fcst .ne. alpha_fcst .or. alpha_anal .ne. alpha_anal) THEN
                                forecast_x25(ktr,ixa,jya) = fcst  
                            ELSE   
                                CALL gamma_quantile_map(fcst, alpha_fcst, &
                                    beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                    fz_anal, precip_qmapped)       
                                forecast_x25(ktr,ixa,jya) = precip_qmapped
                            ENDIF
                            forecast_x25(ktr,ixa,jya) = precip_qmapped 
						ENDIF
					
                    ELSE 
                        
                        ! ---- stencil point is either outside the area where quantile mapping possible.
                        !      As a substitute, take the quantile-mapped value at center point and
                        !      add a small amount of random noise.
                        
                        fcst = forecast(ixa,jya) 
                        !print *,'fcst = ',fcst
						IF (forecast(ixa,jya) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = 0.0
                        ELSE
                        
                            alpha_fcst = gamma_shape_qmap_forecast(ixa,jya)
                            beta_fcst = gamma_scale_qmap_forecast(ixa,jya)
                            fz_fcst = fraction_zero_qmap_forecast(ixa,jya)
                            alpha_anal = gamma_shape_qmap_analysis(ixa,jya)
                            beta_anal = gamma_scale_qmap_analysis(ixa,jya)
                            fz_anal = fraction_zero_qmap_analysis(ixa,jya)
                            
                            !PRINT *,'alpha_fcst, beta_fcst, fz_fcst = ', &
                            !    alpha_fcst, beta_fcst, fz_fcst
                            !PRINT *,'alpha_anal, beta_anal, fz_anal = ', &
                            !    alpha_anal, beta_anal, fz_anal   
                            IF (alpha_fcst .ne. alpha_fcst .or. alpha_anal .ne. alpha_anal) THEN
                                forecast_x25(ktr,ixa,jya) = fcst  
                            ELSE 
                                CALL gamma_quantile_map(fcst, alpha_fcst, &
                                    beta_fcst, fz_fcst, alpha_anal, beta_anal, &
                                    fz_anal, precip_qmapped)
                                forecast_x25(ktr,ixa,jya) = precip_qmapped
                            ENDIF
                            r = 1 + 0.1*(ran1(idum) - 0.5)
                            forecast_x25(ktr,ixa,jya) = forecast_x25(ktr,ixa,jya)*r 
                        ENDIF  
                        
                 	END IF ! jyn>1, jyn < nya, etc.
                    
                    !PRINT *,'forecast_x25(ktr,ixa,jya) = ', forecast_x25(ktr,ixa,jya)
                    
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

REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99
    
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

!print *, 'fcst = ', fcst
IF (nonexceedance_prob .gt. fz_anal) THEN
    cum = (nonexceedance_prob - fz_anal) / (1.-fz_anal)
    
    IF (cum .gt. 0.95) THEN 

        ! ---- with minimal training data, estimation of the correction
        !      at the tails of the distribution may be error prone.
        !      Accordingly, if we are above the 95th percentile of the 
        !      distribution, we will apply an alternative procedure to
        !      direct quantile mapping.  We will estimate forecast and 
        !      analyzed values associated with the 95th, 96, 97, 98, and
        !      99th percentiles of the distribution, and use this
        !      data to form a regression relationship.  The actual
        !      correction applied will be a regression correction based
        !      on this data.   First step here is to get the 95-99th
        !      percentiles, analyzed and forecast        

        CALL gamma_get_95_to_99(alpha_fcst, beta_fcst, fz_fcst, &
            alpha_anal, beta_anal, fz_anal, a95_to_a99, f95_to_f99)

        !      Determine the regression slope associated with the correction
        !      following Schuerer's method discussed in the article
        !      http://journals.ametsoc.org/doi/pdf/10.1175/MWR-D-15-0061.1
        !      Apply that assuming regression intercept
        !      is set by the analyzed value at 95th percentile.

        IF (SUM((f95_to_f99(2:5) - f95_to_f99(1))**2) .gt. 0.0) THEN
       	    slope = SUM((a95_to_a99(2:5) - a95_to_a99(1)) * &
    		    (f95_to_f99(2:5) - f95_to_f99(1))) / &
    	        SUM((f95_to_f99(2:5) - f95_to_f99(1))**2)
            qmapped_fcst = a95_to_a99(1) + (fcst-f95_to_f99(1))*slope
            
            IF (qmapped_fcst .gt. 25.) THEN            
                rmaxmult = 1 + 0.5*exp(-(qmapped_fcst-25.)**2/625.)
                IF (fcst .gt. 25.0) rmaxmult = MIN(rmaxmult, qmapped_fcst/fcst)
                qmapped_fcst = fcst*rmaxmult
            ENDIF
            
            IF (qmapped_fcst .lt. 0.0) qmapped_fcst = 0.0
        ELSE
            qmapped_fcst = fcst
        ENDIF 
    ELSE
        qmapped_fcst = qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse)
        IF (qmapped_fcst .lt. 0.0) qmapped_fcst = 0.0
    ENDIF
    IF (qmapped_fcst .ne. qmapped_fcst) qmapped_fcst = fcst

ELSE ! map to zero
    qmapped_fcst = 0.0

ENDIF
!print *,'qmapped_fcst = ', qmapped_fcst

RETURN
END SUBROUTINE gamma_quantile_map

! ==================================================================================

SUBROUTINE gamma_get_95_to_99(alpha_fcst, beta_fcst, fz_fcst, &
    alpha_anal, beta_anal, fz_anal, a95_to_a99, f95_to_f99)
    
DOUBLE PRECISION, INTENT(IN) :: alpha_fcst, beta_fcst, fz_fcst
DOUBLE PRECISION, INTENT(IN) :: alpha_anal, beta_anal, fz_anal 

DOUBLE PRECISION ksi, cum, ccum, nonexceedance_prob, qgamma
LOGICAL tootrue
LOGICAL toofalse

REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99
    
tootrue = .TRUE.
toofalse = .FALSE.
    
cum = 0.95   
a95_to_a99(1) = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))
f95_to_f99(1) = REAL(qgamma(cum,alpha_fcst,beta_fcst,tootrue,toofalse))

cum = 0.96   
a95_to_a99(2) = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))
f95_to_f99(2) = REAL(qgamma(cum,alpha_fcst,beta_fcst,tootrue,toofalse))

cum = 0.97   
a95_to_a99(3) = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))
f95_to_f99(3) = REAL(qgamma(cum,alpha_fcst,beta_fcst,tootrue,toofalse))

cum = 0.98   
a95_to_a99(4) = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))
f95_to_f99(4) = REAL(qgamma(cum,alpha_fcst,beta_fcst,tootrue,toofalse))

cum = 0.99   
a95_to_a99(5) = REAL(qgamma(cum,alpha_anal,beta_anal,tootrue,toofalse))
f95_to_f99(5) = REAL(qgamma(cum,alpha_fcst,beta_fcst,tootrue,toofalse))

RETURN
END SUBROUTINE gamma_get_95_to_99

    

