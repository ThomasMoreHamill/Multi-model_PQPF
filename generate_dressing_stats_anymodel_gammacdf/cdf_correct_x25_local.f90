SUBROUTINE cdf_correct_x25_local(nxa, nya, npct, nstride, n25, thresh, &
    conusmask, forecast_cdf, analysis_cdf, forecast, forecast_x25)

! --- Perform the CDF bias correction; here we use 8 surrounding grid points forecasts as well
!     to increase sample size and account for position error

INTEGER, INTENT(IN) :: nxa, nya, npct, nstride, n25
REAL, INTENT(IN), DIMENSION(npct) :: thresh
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa,nya,npct) :: forecast_cdf, analysis_cdf
REAL, INTENT(IN), DIMENSION(nxa,nya) :: forecast
REAL, INTENT(OUT), DIMENSION(n25,nxa,nya) :: forecast_x25

REAL :: weight, weight2, rmean
INTEGER :: idum,itest,jtest, jyn, ixn

REAL*8, DIMENSION(npct) :: fcdf, acdf ! ---- set random seed from mean forecast


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
		
        !PRINT *, '********** ixa, jya, conusmask = ', ixa, jya, conusmask(ixa,jya)
		IF (conusmask(ixa,jya) .le. 0) THEN
         	! If outside the conus, then make the output forecast array simply replicates
         	! of the input forecast array
         	forecast_x25(:,ixa,jya) = forecast(ixa,jya)
      	ELSE
         	ktr = 0

         	! Loop thru the 25 grid points with  (ixa,jya) in center

         	DO jyn = jya-2*nstride, jya+2*nstride, nstride
            	DO ixn = ixa-2*nstride, ixa+2*nstride, nstride
					
					ktr = ktr + 1

               	 	! Only process if this point is legitimately within the domain
                    
               	 	IF ( (jyn .ge. 1 .and. jyn .le. nya) .and. &
			   	 	(ixn .ge. 1 .and. ixn .lt. nxa)) THEN


						IF (forecast(ixn,jyn) .eq. 0.0) THEN
							forecast_x25(ktr,ixa,jya) = forecast(ixn,jyn)
							GOTO 3000
						ENDIF
						
                  		IF (conusmask(ixn,jyn) .eq. 1 .and. forecast(ixn,jyn) .ge. 0.0) THEN

                     		! ---- Perform quantile mapping

                  		  	fcdf(:) = forecast_cdf(ixn,jyn,:)
                  		  	acdf(:) = analysis_cdf(ixa,jya,:)
                
                            CALL quantile_map (npct, thresh, forecast(ixn,jyn), fcdf, acdf, &
                                forecast_x25(ktr,ixa,jya))        
                                
                                    
                  	 	ELSE ! conusmask, forecast>0
							forecast_x25(ktr,ixa,jya) = forecast(ixn,jyn)
						ENDIF
						3000 CONTINUE
                        
                        IF (forecast_x25(ktr,ixa,jya) .lt. -99.96) THEN
                            PRINT *,'forecast(ixn,jyn) = ', forecast(ixn,jyn), ixn, jyn
                            stop
                        ENDIF
					
                    ELSE	! jyn>1, jyn < nya, etc.
                        
                        ! --- fill this forecast_x25 value with a slightly perturbed 
                        !     quantile-mapped center grid point value
                        !     so that all points have some sort of value for each of the 
                        !     9 elements surrounding the grid point of interest.
                        
                  		IF (conusmask(ixa,jya) .eq. 1 .and. forecast(ixa,jya) .ge. 0.0) THEN
              		  	    fcdf(:) = forecast_cdf(ixa,jya,:)
              		  	    acdf(:) = analysis_cdf(ixa,jya,:)
                            CALL quantile_map (npct, thresh, forecast(ixa,jya), fcdf, acdf, &
                                forecast_x25(ktr,ixa,jya))
                            r = 1 + 0.1*(ran1(idum) - 0.5)
                            forecast_x25(ktr,ixa,jya) = forecast_x25(ktr,ixa,jya)*r    
                        ELSE
                            r = 1 + 0.1*(ran1(idum) - 0.5)
                            forecast_x25(ktr,ixa,jya) = forecast(ixa,jya)*r    
                        ENDIF                 
                        
                 	END IF ! jyn>1, jyn < nya, etc.
            	END DO ! ixn
         	END DO ! jyn
      	END IF ! conusmask
   	END DO ! ixa
END DO ! jya

RETURN
END SUBROUTINE cdf_correct_x25_local
