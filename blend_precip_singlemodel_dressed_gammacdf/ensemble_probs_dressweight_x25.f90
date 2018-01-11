SUBROUTINE ensemble_probs_dressweight_x25 (n25, nxa, nya, nens, &
    nmembersx25, npcatvals, npvals, nlo_int_hi_vals, &
    nclim_p1_vals, nclim_vals, nthreshes, pthreshes, precip_values, &
    ensemble_ccpa_x25, closest_histogram, precip_histogram_thresholds, &
    gamma_shapes, gamma_scales, fraction_zeros, gamma_shape_fclimpop, &
    gamma_scale_fclimpop, fraction_zeros_fclimpop, climo_pop_thresholds, &
    conusmask, climo_prob, prob_forecast_qmapped, prob_forecast)  
    
INTEGER, INTENT(IN) :: n25, nxa, nya, nens, nmembersx25
INTEGER, INTENT(IN) :: npcatvals, npvals, nlo_int_hi_vals
INTEGER, INTENT(IN) :: nclim_p1_vals, nclim_vals, nthreshes



REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes
REAL, INTENT(IN), DIMENSION(npvals) :: precip_values
REAL, INTENT(IN), DIMENSION(n25,nxa,nya,nens) :: ensemble_ccpa_x25
REAL, INTENT(IN), DIMENSION(nmembersx25,npcatvals) :: closest_histogram
REAL, INTENT(IN), DIMENSION(npcatvals) :: precip_histogram_thresholds
REAL, INTENT(IN), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: gamma_shapes ! shape parameters for 
    ! Gamma distributions
REAL, INTENT(IN), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: gamma_scales ! scale parameters for 
    ! Gamma distributions
REAL, INTENT(IN), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: fraction_zeros ! fraction of 
    ! samples with zero precip.
    !npcatvals = 3 : closest_histogram stats are stratified by ens-mean amount.
    !nlo_int_hi_vals = index for lowest, intermediate, highest member of closest histogram
    
    
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: gamma_shape_fclimpop ! when ens mean = 0., 
    ! shape parameter of nonzero dressed values.
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: gamma_scale_fclimpop ! when ens mean = 0., 
    ! scale parameter of nonzero dressed values.
REAL, INTENT(IN), DIMENSION(nclim_p1_vals) :: fraction_zeros_fclimpop ! when ens mean = 0., 
    ! fraction of samples with zero dressed values
REAL, INTENT(IN), DIMENSION(nclim_vals) :: climo_pop_thresholds 
    ! climatological probability thresholds 
    ! between elements of gamma_shape_fclimpop, etc.
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: climo_prob
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: prob_forecast_qmapped
REAL, INTENT(OUT), DIMENSION(nxa,nya,nthreshes) :: prob_forecast

REAL, DIMENSION(n25, nens) :: ensemble
REAL, DIMENSION(nmembersx25) :: ensemble_sorted
REAL*8, DIMENSION(nmembersx25) :: weights, weight_lower, weight_upper
REAL*8 sumprobability
REAL*8 alpha, beta, cum, ccum, ksi, zval
INTEGER, DIMENSION(1) :: ishape

ishape(1) = nmembersx25
prob_forecast(:,:,:) = -99.99
iprint = 1

! ---- process all precipitation thresholds

!PRINT *,'precip_histogram_thresholds(1) = ',precip_histogram_thresholds(1)

DO ithresh = 1, nthreshes
    
    rthresh = pthreshes(ithresh)

    ! ---- process all grid points inside mask
    
    DO jya = 1, nya
        DO ixa = 1, nxa

    !DO jya = 52, 52
    !    DO ixa = 123, 123

            !IF (iprint .eq. 1) &
            !PRINT *,'====== PROCESSING ixa, jya, mask, rthresh = ', &
            !    ixa, jya, conusmask(ixa,jya), rthresh
            IF (conusmask(ixa,jya) .eq. 1) THEN
                sumprobability = 0.
                ensemble(:,:) = ensemble_ccpa_x25(:,ixa,jya,:)       
                ensemble_sorted = RESHAPE(ensemble,ishape)
                emean = SUM(ensemble) / REAL(n25*nens)
                CALL sort(nmembersx25, ensemble_sorted)   
                emin = MINVAL(ensemble_sorted) 
                emax = MINVAL(ensemble_sorted)             
                CALL determine_category(emean, npcatvals, &
                    precip_histogram_thresholds, ipcat)
                ipcat = ipcat-1 ! reduce index b/c we don't care about gamma dist for near zero
                !IF (iprint .eq. 1) PRINT *, ' ixa, jya, emean, ipcat = ', ixa, jya, emean, ipcat 
                !IF (ixa .eq. 352 .and. jya .eq. 2) THEN
                !PRINT *,'ensemble_sorted = ', ensemble_sorted
                !ENDIF

                IF (emean .lt. precip_histogram_thresholds(1)) THEN    
                    
                    ! NEVER OCCURS, SO TURN OFF
                    
                    ! --- we have near zero ens-mean precip.  We will estimate
                    !     probabilities from the _fclimpop arrays, which
                    !     indicate the fraction zero and gamma parameters as 
                    !     a function of the climatological POP.
                    
                    CALL determine_category(climo_prob(ixa,jya,1), &
                        nclim_vals, climo_pop_thresholds, iclim)
                    alpha = gamma_shape_fclimpop(iclim) 
                    beta = gamma_scale_fclimpop(iclim) 
                    fz = fraction_zeros_fclimpop(iclim) 
                    ksi = rthresh / beta
                    IF (alpha .ne. alpha) THEN ! happens only when NaN
                        sumprobability = 0.
                    ELSE
                        CALL cumgam(ksi, alpha, cum, ccum)
                        sumprobability = (1.-fz)*ccum
                    ENDIF
                    
                ELSE ! emean .ge. precip_histogram_thresholds(1)
                    
                    ! --- we have nonzero mean precipitation.  We will form
                    !     a probability estimate of exceeding the threshold of 
                    !     interest from a weighted linear combination of
                    !     exceedance probabilities estimated from the Gamma
                    !     distribution parameters appropriate to that sorted
                    !     ensemble member number.
                    
                    weights(:) = closest_histogram(:,ipcat)
                    !print *,'weights = ', weights
                    !print *,'nmembersx25  = ', nmembersx25
                    
                    DO imem = 1, nmembersx25
                        
                        ! ---- determine the bounding lower precipitation value
                        !      in the precip_values array, and the weight
                        !      to apply to this lower value
                        
                        !IF (iprint .eq. 1)&
                        !   PRINT *,'-- processing imem = ',imem,&
                        !    ' quantile mapped value = ',ensemble_sorted(imem)                        
                        IF (imem .eq. 1) THEN
                            ilomidhi = 1  ! index for lowest sorted member
                        ELSE IF (imem  .eq. nmembersx25) THEN
                            ilomidhi = 3  ! index for highest
                        ELSE
                            ilomidhi = 2  ! index for intermediate
                        ENDIF
                        
                        CALL determine_precip_value_category(npvals, precip_values, &
                            ensemble_sorted(imem), ipvlower)
                        ipvupper = ipvlower + 1
                        !IF (iprint .eq. 1) &
                            !PRINT *,'   ipvlower, ipcat, ilomidhi, ens = ',&
                            !ipvlower, ipcat, ilomidhi, ensemble_sorted(imem)
                        !IF (iprint .eq. 1) &
                            !PRINT *,'   precip_values lo/up',precip_values(ipvlower),&    
                            !precip_values(ipvupper)
                        weight_gamma_to_lower = 1.0 - &
                            (ensemble_sorted(imem) - precip_values(ipvlower)) / &
                            (precip_values(ipvupper) - precip_values(ipvlower)) 
                        IF (weight_gamma_to_lower .lt. 0.0) weight_gamma_to_lower = 0.0
                        IF (weight_gamma_to_lower .gt. 1.0) weight_gamma_to_lower = 1.0
                        weight_gamma_to_upper = 1.0 - weight_gamma_to_lower
                        !IF (iprint .eq. 1) &
                        !    PRINT *,'   weight_gamma_to_lower = ', weight_gamma_to_lower
                    
                        ! ---- linearly interpolate the fraction zero
                        !      to a value appropriate for this precip amount.
                                            
                        fz = weight_gamma_to_lower * fraction_zeros(ipvlower, ipcat, ilomidhi) + &
                            weight_gamma_to_upper * fraction_zeros(ipvupper, ipcat, ilomidhi)
                        IF (fz .lt. 0.0) THEN  ! -99.99, missing; use nearest value
                            DO ip = ipvlower, 1, -1
                                IF (fraction_zeros(ip,ipcat,ilomidhi) .ge. 0.0) THEN
                                    fz = fraction_zeros(ip,ipcat,ilomidhi)
                                    GOTO 1001
                                ENDIF
                            END DO
                            
                            DO ip = ipvupper, npvals, 1
                                IF (fraction_zeros(ip,ipcat,ilomidhi) .ge. 0.0) THEN
                                    fz = fraction_zeros(ip,ipcat,ilomidhi)
                                    GOTO 1001
                                ENDIF
                            END DO
                            1001 CONTINUE    
                        ENDIF

                        !IF (iprint .eq. 1) &
                        !    PRINT *,'   ipvlower, ipcat, ilomidhi, fz = ', &
                        !    ipvlower, ipcat, ilomidhi, fz 
                            
                        IF (ilomidhi .eq. 1 .or. ilomidhi .eq. 3) THEN
                            
                            ! ---- linearly interpolate the gamma distribution parameters 
                            !      to a value appropriate for this precip amount.
                            
                            alpha = weight_gamma_to_lower * gamma_shapes(ipvlower, ipcat, ilomidhi) + &
                                weight_gamma_to_upper * gamma_shapes(ipvupper, ipcat, ilomidhi)
                            beta = weight_gamma_to_lower * gamma_scales(ipvlower, ipcat, ilomidhi) + &
                                weight_gamma_to_upper * gamma_scales(ipvupper, ipcat, ilomidhi) 
                            ccum = 0.0
                            IF (alpha .gt. 0.0) THEN  ! valid dressing stats exist
                                IF (fz .lt. 1.0) THEN
                                    ksi = rthresh / beta
                                    CALL cumgam(ksi, alpha, cum, ccum)  ! ccum is exceedance probability
                                    !CALL determine_exceedance_probability_normal(nmembersx25, imem, &
                                    !    weights, ensemble_sorted(imem), rthresh, ccum)
                                ENDIF
                                IF (ensemble_sorted(1) .eq. 0) THEN
                                    zero_fzsave = fz
                                    zero_ccum = ccum
                                ENDIF
                            ELSE
                                ! ---- for lack of training data to set Gamma distribution parameters, 
                                !      put a kernel of probability density centered on
                                !      the quantile-mapped value, with a small std dev.
                                
                                CALL determine_exceedance_probability_normal(nmembersx25, imem, &
                                    weights, ensemble_sorted(imem), rthresh, ccum) 
                            ENDIF
                        ELSE ! ilomidhi =  2; intermediate rank 
                            
                            IF (ensemble_sorted(imem) .gt. 0.0) THEN
                                CALL determine_exceedance_probability_normal(nmembersx25, imem, &
                                    weights, ensemble_sorted(imem), rthresh, ccum) 
                            ELSE
                                fz = 1.0
                                ccum = 0.0
                                !fz = zero_fzsave
                                !ccum = zero_ccum
                            ENDIF                                
                        ENDIF
                        sumprobability = sumprobability + weights(imem)*(1.-fz)*ccum
                        !IF (iprint .eq. 1) &
                        !    PRINT 206,'  sumprob, ens, ccum, fz, wgt = ',&
                        !    sumprobability, ensemble_sorted(imem), ccum, fz, &
                        !    weights(imem)
                        !206 format(a37,5(f10.6,1x))
                        
                    END DO ! imem
                    IF (sumprobability .gt. 0.9999) sumprobability = 1.0
                    IF (sumprobability .lt. 0.0) sumprobability = 0.0
                    !IF (iprint .eq. 1) &
                    !    PRINT *,'   sumprobability = ', sumprobability
                    
                ENDIF ! use climatological probability or gamma/normal distributions.
                prob_forecast(ixa,jya,ithresh) = sumprobability
            ELSE
                prob_forecast(ixa,jya,ithresh) = prob_forecast_qmapped(ixa,jya,ithresh)
            ENDIF ! conusmask
            
            !IF (prob_forecast(ixa,jya,ithresh) .lt. 0.0) THEN
            !    PRINT *,'ixa, jya, ithresh, prob_forecast ',&
            !        ixa,jya,ithresh,prob_forecast(ixa,jya,ithresh)
            !ENDIF
            
        END DO ! ixa
    END DO ! jya
    
    !PRINT *,'prob_forecast(3*nxa/4,nya/2,:) = ', prob_forecast(3*nxa/4,nya/2,:) 
    !stop
    
    write (6,fmt='(A40,f8.4)') 'Dressed probabilities for thresh = ',rthresh
    write (6,fmt='(3(A10,1X))')'MIN','MAX','MEAN'
    write (6,fmt='(1X,3(F10.5,1X))')&
        minval(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        maxval(conusmask(:,:)*prob_forecast(:,:,ithresh)), &
        SUM(prob_forecast(:,:,ithresh))/REAL(nxa*nya)
        !sum(conusmask(:,:)*prob_forecast(:,:,ithresh))/REAL(SUM(conusmask))
    
END DO ! ithresh = 1, nthreshes

!PRINT *,'prob_forecast(3*nxa/4,nya/2,:) = ', prob_forecast(3*nxa/4,nya/2,:)

RETURN
END SUBROUTINE ensemble_probs_dressweight_x25                     
         
! --------------------------------------------------------------------------
       
SUBROUTINE determine_exceedance_probability_normal(nmembersx25, imem, weights, &
    ensemble_member, rthresh, ccum)               
INTEGER, INTENT(IN) :: nmembersx25, imem
REAL*8, DIMENSION(nmembersx25) :: weights
REAL, INTENT(IN) :: ensemble_member, rthresh
REAL*8 zval, cum, ccum
                
! ---- for interior ranks in the sorted ensemble or for a sample where
!      for lack of training data we do not have Gamma distribution 
!      information, we assume that a narrow Gaussian probability 
!      distribution is appropriate, centered on the quantile-mapped value.   
!      The standard deviation of the Gaussian distribution is a product 
!      of three factors:
!      (a) the inverse of the ensemble size.  A larger ensemble
!          with more members will have members more tightly squeezed
!          together, and hence the distribution of analyzed around the 
!          best-member forecast will be smaller;
!      (b) the ratio of the weight (closest_histogram) value for this 
!          sorted member divided by the minimum weight value. In
!          this way, near the extreme ranks where the weight is larger,
!          this means effectively that there is a bit more distance
!          to the next-closest member, and that should be accounted for
!          in the dressing distribution; and 
!      (c) the value of the quantile-mapped member itself.  Generally
!          uncertainty is proportional to amount 

wmin = MINVAL(weights)
chratio = weights(imem) / wmin
!denom = (1./REAL(nmembersx25)) * chratio * ensemble_member
denom = ensemble_member / (REAL(nmembersx25)/20.)

IF (denom .lt. 0.001) denom = 0.001
zval = (rthresh - ensemble_member) / denom ! convert to std normal deviate
CALL cumnor(zval, cum, ccum)  ! ccum is the exceedance probability
RETURN
END SUBROUTINE determine_exceedance_probability_normal