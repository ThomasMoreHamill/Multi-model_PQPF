SUBROUTINE quantile_map (npct, thresh, forecast, fcdf, acdf, &
    forecast_interpolated)

! --- perform the quantile mapping

INTEGER, INTENT(IN) :: npct
REAL, INTENT(IN), DIMENSION(npct) :: thresh 
REAL*8, INTENT(IN), DIMENSION(npct) :: fcdf, acdf
REAL, INTENT(IN) :: forecast
REAL, INTENT(OUT) :: forecast_interpolated

REAL, DIMENSION(5) :: a95_to_a99, f95_to_f99


DO ipct = 1, npct-1
    
    IF (ipct .eq. 1) THEN
        tbelow = 0.0
        tabove = thresh(ipct)
        fbelow = 0.0
        fabove = fcdf(ipct)
    ELSE 
        tbelow = thresh(ipct)
        tabove = thresh(ipct+1)
        fbelow = fcdf(ipct)
        fabove = fcdf(ipct+1)
    ENDIF
    
	IF (forecast .ge. tbelow .and. forecast .lt. tabove) THEN

 		!  ---- Determine the percentile in the 
		!       CDF associated with this forecast amount
		
	 	weight = (forecast-tbelow) / (tabove -  tbelow)
		cdf_interpolated = fbelow*(1.-weight) + fabove*weight

		! ---- For non-extreme values, set the new forecast to be the analysis 
	 	! 	   value associated with the same quantile of the cdf.  

	 	IF (cdf_interpolated .lt. 0.95)THEN
  			weight2 = 0.0
  			IF (cdf_interpolated .lt. acdf(1)) THEN
     		    forecast_interpolated = 0.0
     		    GOTO 3000
  		  	END IF

  		    DO icdf = 1,npct-1
     			IF (cdf_interpolated .ge. acdf(icdf) .and. &
        		cdf_interpolated .lt. acdf(icdf+1)) THEN
        			weight2 = (cdf_interpolated-acdf(icdf)) / &
                        (acdf(icdf+1)-acdf(icdf))
        		    forecast_interpolated = &
					    thresh(icdf)*(1.-weight2) + thresh(icdf+1)*weight2
        		    GOTO 3000
     		    END IF
  		 	END DO
	 	ELSE

  			! cdf_interpolated >= 0.95; apply Scheuerer regression analysis approach
  		  	! from appendix A of Scheuerer and Hamill, MWR, 143, 4578-4596. The 
  		  	! underlying rationale is that quantile mapping produces potentially
  		  	! especially unrealistic values at the extreme high percentiles, so the
  		  	! regression approach should diminish this tendency.

  		  	! Find the forecast and analyzed values associated with the 95th
  		  	! thru 99th percentiles of the distribution.
  
  		  	CALL get_95_to_99(npct, acdf, fcdf, thresh, a95_to_a99, f95_to_f99)

  		  	! Determine the regression slope associated with the correction
  		  	! following Schuerer's method.  Apply that assuming regression intercept
  		  	! is set by the analyzed value at 95th percentile.

  		  	IF (SUM((f95_to_f99(2:5) - f95_to_f99(1))**2) .gt. 0.0) THEN
     		   	slope = SUM((a95_to_a99(2:5) - a95_to_a99(1)) * &
					(f95_to_f99(2:5) - f95_to_f99(1))) / &
           		SUM((f95_to_f99(2:5) - f95_to_f99(1))**2)
     		    forecast_interpolated = a95_to_a99(1) + &
					(forecast-f95_to_f99(1))*slope

                IF (forecast_interpolated .gt. 25.) THEN            
                    rmaxmult = 1 + 0.5*exp(-(forecast_interpolated-25.)**2/625.)
                    IF (forecast .gt. 25.0) rmaxmult = MIN(rmaxmult, forecast_interpolated/forecast)
                    forecast_interpolated = forecast*rmaxmult
                ENDIF
                    
                    
  			ELSE
     		    forecast_interpolated = forecast
  		  	ENDIF 
            
            
            
            
            

            
            
            

  		  	GOTO 3000

	 	ENDIF ! cdf_interpoloated < 0.95
	ENDIF  ! forecast >= thresh(ipct), < thresh(ipct)+1
END DO ! ipct

3000 CONTINUE
RETURN
END SUBROUTINE quantile_map
