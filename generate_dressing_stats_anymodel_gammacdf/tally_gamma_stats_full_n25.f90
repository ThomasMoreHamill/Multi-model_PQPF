SUBROUTINE tally_gamma_stats_full_n25 (n25, nxa, nya, nmembers, &
    n_climocats, nthreshes, nout_thresh, thresh_light, &
    thresh_mod, thresh_high, output_threshes, &
    ensemble_x25, analysis, conusmask, outfilename_nc, &
    gamma_threshes, climo_prob, climo_pop_thresholds, istat)

    ! this subroutine is the heart of the generate_dressing_stats_anymodel_gammacdf.x
    ! Here is where we actually compare ensemble forecasts for a given prediction 
    ! system to the verifying analyses, tallying the information we need
    ! to generate (later, over many case days) the closest-histogram statistics
    ! the fraction of samples with zero precipitation, and the gamma shape and
    ! scale parameters.  These last three are function of the precipitation amount
    ! and have dependence on whether the member is the lowest sorted member,
    ! an intermediate, or highest member.   There is also a dependence for both
    ! the closest histogram and the gamma-distribution statistics on the 
    ! ensemble-mean precipitation amount.   Data is saved to a netcdf file
    ! for this particular model and forecast initial time and forecast lead.
    !
    ! Gamma-distribution parameters are not estimated here.  Rather the sums
    ! necessary (Wilks 2011 text, Statistical Methods in the Atmospheric Sciences,
    ! 3rd Ed, eq. 4.40) are generated so that at a later step we can generate
    ! ln(xbar) and (1/n)sum(ln(xi)), and from that the Gamma shape and scale
    ! parameters.
    
    ! coded by: Tom Hamill, ESRL/PSD, tom.hamill@noaa.gov (303) 497-3060

USE netcdf
    
INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
INTEGER, INTENT(IN) :: n_climocats ! the number of categories for climatology of POP
INTEGER, INTENT(IN) :: nthreshes ! dimension for discretization of precip amt
INTEGER, INTENT(IN) :: nout_thresh
REAL, INTENT(IN) :: thresh_light, thresh_mod, thresh_high ! for breaking up closest histogram etc
REAL, INTENT(IN), DIMENSION(nout_thresh) :: output_threshes
REAL, INTENT(IN), DIMENSION (n25, nxa, nya, nmembers) :: ensemble_x25 ! qmapped ensemble
REAL, INTENT(IN), DIMENSION (nxa, nya) :: analysis ! analyzed precip amt.
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask 
CHARACTER*(*), INTENT(IN) :: outfilename_nc
REAL, INTENT(IN), DIMENSION (nxa, nya) :: climo_prob ! climatological probability of POP.
REAL, INTENT(IN), DIMENSION(n_climocats-1) :: climo_pop_thresholds ! used for setting
    ! gamma distribution parameters in the eventuality that all ensemble members are zero.
    ! Then we permit a dependence of the nonzero amounts on the precipitation climatology.
REAL, INTENT(IN), DIMENSION(nthreshes) :: gamma_threshes ! values where gamma distribution
    ! sums are to be defined.

REAL, INTENT(OUT) :: istat

! ---- Here are the variables we intend to populate with the forecast, analysis information.
!      dimension 3 below is for zero, light, mod, heavy precip, with thresholds between
!      them set by inputs thresh_light, thresh_mod, thresh_high.
!      dimension 3 below is for lowest sorted member(1), intermediate(2), highest(3) 

REAL*8, DIMENSION(nthreshes,3,3) :: gamma_sum 
REAL*8, DIMENSION(nthreshes,3,3) :: gamma_ln_sum
INTEGER, DIMENSION(nmembers*n25, 3) :: closest_histogram
REAL*8, DIMENSION(nthreshes, 3, 3) :: nzeros
REAL*8, DIMENSION(nthreshes, 3, 3) :: npositive
REAL*8, DIMENSION(nthreshes, 3, 3) :: sum_of_weights
REAL*8, DIMENSION(n_climocats) :: gamma_sum_ensmeanzero_fclim
REAL*8, DIMENSION(n_climocats) :: gamma_ln_sum_ensmeanzero_fclim
REAL*8, DIMENSION(n_climocats) :: nzeros_fclim
REAL*8, DIMENSION(n_climocats) :: npositive_fclim
INTEGER, DIMENSION(n_climocats,nout_thresh) :: exceed_yes, exceed_no

INTEGER, DIMENSION(nthreshes) :: ilower_bound, iupper_bound 
    ! range of gamma_threshes indices we will increment sums 
    ! at when presented with a gamma_threshes index of ipclosest.
REAL*8, DIMENSION(nthreshes, nthreshes) :: weightfn

REAL rdum 
INTEGER idum

integer :: dimid_3d(3), dimid_2d(2), dimid_2da(2), dimid_3da(3)

print *, 'subroutine tally_gamma_stats_full_n25'

istat = 1

! ----- initialize to zero

gamma_sum(:,:,:) = 0.0
gamma_ln_sum(:,:,:) = 0.0
closest_histogram(:,:) = 0.0
nzeros(:,:,:) = 0.
npositive(:,:,:) = 0.
sum_of_weights(:,:,:) = 0.0
gamma_sum_ensmeanzero_fclim(:) = 0.
gamma_ln_sum_ensmeanzero_fclim(:) = 0.
nzeros_fclim(:) = 0.
npositive_fclim(:) = 0.
exceed_no(:,:) = 0
exceed_yes(:,:) = 0

! ---- define bounds for gamma_threshes summation and weighting function

CALL define_indices_weights(nthreshes, gamma_threshes, &
    ilower_bound, iupper_bound, weightfn)

! ---- set a random number seed for the closest rank histogram

rdum = 0.0
DO ixa = 1, nxa
	rdum = rdum + ensemble_x25(1,ixa,nya/2,1)
END DO
idum = NINT(rdum)

! ---- only process if input analysis, ensemble have valid data ....

rm = MINVAL(ensemble_x25)
rma = MAXVAL(analysis)
IF (rm .lt. -98. .and. rma .lt. -98.) THEN
	istat = -1
	PRINT *,'identified bad ensemble or analysis data.  rm, rma = ',rm, rma
	PRINT *,'setting istat to -1 in tally_gamma_stats_full_x25 and exiting subroutine.'
ELSE
    
    ! --- loop thru and process all points in the CONUS
    
	DO jya = 1, nya
		DO ixa = 1, nxa
			
			IF (conusmask(ixa, jya) .eq. 1 .and. analysis(ixa,jya) .ge. 0.0) THEN
            
				! --- determine which member is the closest to the
				!     analyzed and how many members have values lower 
				!     than or equal to the analyzed value
			
				rsum = 0.0
				iktr = 0
				DO i25 = 1, n25
					DO imem = 1, nmembers
						IF (ensemble_x25(i25,ixa,jya,imem) .ge. 0.0) THEN
							rsum = rsum + ensemble_x25(i25,ixa,jya,imem)
							iktr = iktr+1
						ENDIF
					END DO
				END DO
				emean = rsum / REAL(iktr)
                a = analysis(ixa, jya)

                ! ---- we will tally stats one of two ways below.  For situations
                !      where the precipitation is heavier than a threshold that
                !      is slightly greater than zero, we keep track of closest
                !      histogram, fraction zero, and Gamma alpha and beta parameters.
                !      In situations where the mean precip is zero or really 
                !      close to zero, we will tally the fraction zero, and the 
                !      Gamma alpha and beta parameters as a function of the 
                !      climatological probability of precipitation.
    
                IF (emean .ge. thresh_light) THEN ! nonzero precip
                    
                    ! ---- determine what the category is associated with this 
                    !      (quantile-mapped) ensemble-mean precipitation amount
                    
                    IF (emean .ge. thresh_light .and. emean .lt. thresh_mod) THEN
                        ipcat = 1
                    ELSE IF (emean .ge. thresh_mod .and. emean .lt. thresh_high) THEN
                        ipcat = 2
                    ELSE IF (emean .ge. thresh_high) THEN
                        ipcat = 3
                    ELSE
                        ipcat = -99
                    ENDIF  
                    
                    ! ----- PART A:  Increment closest histogram ----------------
                    
    				! ---- find the ensemble forecast member that is the closest
                    !      to the analyzed

    				i25closest = 1
    				imemclosest = 1
    				rclosest = 9999.
                    eclosest = 0.
    				DO i25 = 1, n25
    					DO imem = 1, nmembers
    						e = ensemble_x25(i25,ixa,jya,imem)
    						diff = ABS(a - e)
    						IF (diff .lt. rclosest .and. e .gt. -99) THEN
    							rclosest = diff
    							eclosest = e
    							i25closest = i25
    							imemclosest = imem
    						ENDIF
    					END DO
    				END DO

    				! ---- determine how many other members are lower than the
    				!      closest member, and how many are equal,  in order
                    !      to determine the rank of the closest in the sorted
                    !      ensemble.
				
    				ibelow = 0
    				iequal = 0
    				DO i25 = 1, n25
    					DO imem = 1, nmembers
    						e = ensemble_x25(i25,ixa,jya,imem)
    						IF (i25 .eq. i25closest .and. imem .eq. imemclosest) THEN
                                CONTINUE
                            ELSE
    							IF (e .lt. eclosest) ibelow = ibelow + 1
    							IF (e .eq. eclosest) iequal = iequal + 1
    						ENDIF
    					END DO
    				END DO

    				! --- determine the closest_histogram rank, including a randomization procedure
                    !     in the case of ties of ensemble and analysis (common when both = 0)
				
    				IF (iequal .eq. 0) THEN
    					iclosest = ibelow + 1			
    				ELSE
    					r = ran3(idum) * iequal
    					ir = INT(r)
    					IF (ir .gt. iequal) ir = iequal
    					iclosest = ibelow + ir + 1
    				ENDIF                   
                    
                    ! ---- increment closest histogram for nonzero ens mean precip 
                    
                    closest_histogram(iclosest,ipcat) = closest_histogram(iclosest,ipcat) + 1
                    
                    ! ----- PART B:  Increment gamma distribution statistics ----------------
                    
                    ! --- determine index for whether this is lowest, intermediate, or highest
                    !     ranked closest member.  We'll tally gamma-distribution statistics
                    !     separately in each circumstance.

            	    IF (iclosest .eq. 1) THEN
                        ilomidhi = 1
            	    ELSE IF (iclosest .gt. 1 .and. iclosest .lt. nmembers*n25) THEN
                        ilomidhi = 2
            	    ELSE IF (iclosest .eq. nmembers*n25) THEN
            		    ilomidhi = 3
            	    ENDIF

                    ! ---- determine the gamma amount threshold (ipclosest) that is closest
                    !      to this best-member ensemble member (eclosest)

                    ipclosest = 1
                    dmin = 9999.
                    DO icat = 1, nthreshes
                	    diff = ABS(gamma_threshes(icat) - eclosest)
                	    IF (diff .lt. dmin) THEN
                		    dmin = diff
                		    ipclosest = icat
                	    ENDIF
                    END DO

                    ! --- recenter the analyzed value on the precipitation amount associated with 
                    !     the nearest value in the gamma_threshes array

                    pertval = gamma_threshes(ipclosest) + (a-eclosest)
                    IF (pertval .lt. 0.) pertval = 0.
                    
                    IF (pertval .gt. 0) THEN
            
                        ! ---- because of the limited sample size, we are going to allow gamma_threshes
                        !      that are in the neighborhood of eclosest to be updated with the 
                        !      difference between the analysis and eclosest.   We'll scale that
                        !      difference up or down a bit so that precip amounts larger than
                        !      eclosest will have a larger perturbation associated with them.
                    
                        DO ith = ilower_bound(ipclosest), iupper_bound(ipclosest)
                            IF (gamma_threshes(ipclosest) .ne. 0.0) THEN
                                scaling = gamma_threshes(ith) / gamma_threshes(ipclosest)  
                                pertval2 = gamma_threshes(ith) + (a-eclosest)*scaling
                            ELSE
                                scaling = 1.0  
                                pertval2 = a-eclosest
                            ENDIF

                            sum_of_weights(ith,ipcat,ilomidhi) = sum_of_weights(ith,ipcat,ilomidhi) + &
                                weightfn(ipclosest,ith)                    
                    
                            ! ---- increment the summations we'll need to calculate gamma parameters
                            !      and fraction of samples with zero precipitation.
                    
            		        gamma_sum(ith,ipcat,ilomidhi) = &
                                gamma_sum(ith,ipcat,ilomidhi) + weightfn(ipclosest,ith)*pertval2
            		        gamma_ln_sum(ith,ipcat,ilomidhi) = &
                                gamma_ln_sum(ith,ipcat,ilomidhi) + weightfn(ipclosest,ith)*alog(pertval2)
            		        npositive(ith,ipcat,ilomidhi) = &
                                npositive(ith,ipcat,ilomidhi) + weightfn(ipclosest,ith)
                        END DO 
                    ELSE
                        DO ith = ilower_bound(ipclosest), iupper_bound(ipclosest)
                	        nzeros(ith,ipcat,ilomidhi) = nzeros(ith,ipcat,ilomidhi) + weightfn(ipclosest,ith)
                        END DO
                    ENDIF    
                    
                ELSE ! treat as zero precip
                
                    ! ---- for situations where all members' quantile-mapped ensemble forecasts
                    !      are very near or equal to zero, we want to then permit 
                    !      estimating the dressing distribution of possible nonzero
                    !      amounts as a function of the climatological probability of
                    !      nonzero precip.  The assumption is that this probability 
                    !      distribution will vary from climatologically dry to wet regions.  

                    cprob = climo_prob(ixa,jya)
                    CALL find_climo_category (n_climocats-1, climo_pop_thresholds, cprob, iclim)
                    
                    IF (a .gt. 0.0) THEN  
                        gamma_sum_ensmeanzero_fclim(iclim) = &
                            gamma_sum_ensmeanzero_fclim(iclim) + a
                        gamma_ln_sum_ensmeanzero_fclim(iclim) = &
                            gamma_ln_sum_ensmeanzero_fclim(iclim) + alog(a)
                        npositive_fclim(iclim) = npositive_fclim(iclim) + 1. 
                        
                        DO ithresh = 1, nout_thresh
                            IF (a .gt. output_threshes(ithresh)) THEN
                                exceed_yes(iclim,ithresh) = exceed_yes(iclim,ithresh) + 1
                            ELSE
                                exceed_no(iclim,ithresh) = exceed_no(iclim,ithresh) + 1
                            ENDIF
                        END DO  

                    ELSE
                        nzeros_fclim(iclim) = nzeros_fclim(iclim) + 1. 
                    ENDIF
                    

                    
                ENDIF  ! emean .ge. thresh_light	
            ENDIF ! conusmask
		END DO ! ixa
	END DO ! jya

    ! ---- print out the closest histogram sample

    PRINT *,'closest_histogram'
    DO i25 = 1,n25
        istart = 1 + (i25-1)*20
        iend = istart + 19
        PRINT 200, closest_histogram(istart:iend,3)
        200 FORMAT(20(i3,1x))
    END DO
    
    PRINT *,'output_threshes = ', output_threshes
    DO iclim = 1,n_climocats
        PRINT *,'POP category index = ',iclim
        PRINT 208, 'NO  ',exceed_no(iclim,:)
        PRINT 208, 'YES ',exceed_yes(iclim,:)
        208 FORMAT(a4,1x,7(i9,1x))
    END DO
    
    ! ====================== WRITE DAILY OUTPUT TO NETCDF FILE ==========================
    
    ! ---- Create the netCDF file.
    
    print *,'writing to ',TRIM(outfilename_nc)
    CALL check( nf90_create(TRIM(outfilename_nc), NF90_CLOBBER, ncid) )
    
    ! ---- Define the array dimensions. NetCDF will hand back an ID for each.
  
    PRINT *,'array dimensions'
    CALL check( nf90_def_dim(ncid, "nhist", nmembers*n25, nhist_dimid) )
    CALL check( nf90_def_dim(ncid, "nthreshes", nthreshes, nthreshes_dimid) )
    CALL check( nf90_def_dim(ncid, "nccats", n_climocats, nccats_dimid) )
    CALL check( nf90_def_dim(ncid, "nccats_minus1", n_climocats-1, nccats_minus1_dimid) )
    CALL check( nf90_def_dim(ncid, "nscalar", 1, nscalar_dimid) )
    CALL check( nf90_def_dim(ncid, "nprecipcats", 3, nprecipcats_dimid) )
    CALL check( nf90_def_dim(ncid, "nlomidhirank", 3, nlomidhirank_dimid) )
    CALL check( nf90_def_dim(ncid, "nout_thresh", nout_thresh, nout_thresh_dimid) )
    dimid_3d =  (/ nthreshes_dimid, nprecipcats_dimid, nlomidhirank_dimid /)
    dimid_2d =  (/ nhist_dimid, nprecipcats_dimid /)
    dimid_2da = (/ nccats_dimid, nout_thresh_dimid /)

    ! ---- Define the variables and associated IDs

    PRINT *,'defining variables'
    CALL check( nf90_def_var(ncid, "closest_histogram", &
        NF90_INT, dimid_2d, nhist_varid) )
    CALL check( nf90_def_var(ncid, "gamma_sum", &
        NF90_DOUBLE, dimid_3d, ngsum_varid) )
    CALL check( nf90_def_var(ncid, "gamma_ln_sum", &
        NF90_DOUBLE, dimid_3d, nglnsum_varid) )
    CALL check( nf90_def_var(ncid, "weights", &
        NF90_DOUBLE, dimid_3d, nweights_varid) )
    CALL check( nf90_def_var(ncid, "nzeros", &
        NF90_DOUBLE, dimid_3d, nzeros_varid) )
    CALL check( nf90_def_var(ncid, "npositive", &
        NF90_DOUBLE, dimid_3d, npositive_varid) )    
    CALL check( nf90_def_var(ncid, "gamma_sum_ensmeanzero_fclim", &
        NF90_DOUBLE, nccats_dimid, ngamma_sum_ensmeanzero_fclim_varid) )    
    CALL check( nf90_def_var(ncid, "gamma_ln_sum_ensmeanzero_fclim", &
        NF90_DOUBLE, nccats_dimid, ngamma_ln_sum_ensmeanzero_fclim_varid) ) 
    CALL check( nf90_def_var(ncid, "nzeros_fclim", &
        NF90_DOUBLE, nccats_dimid, nzeros_fclim_varid) )    
    CALL check( nf90_def_var(ncid, "npositive_fclim", &
        NF90_DOUBLE, nccats_dimid, npositive_fclim_varid) )     
    CALL check( nf90_def_var(ncid, "exceed_yes_fclimPOP", &
        NF90_INT, dimid_2da, nexceed_yes_fclimPOP_varid) )    
    CALL check( nf90_def_var(ncid, "exceed_no_fclimPOP", &
        NF90_INT, dimid_2da, nexceed_no_fclimPOP_varid) )                          
            
    CALL check( nf90_def_var(ncid, "gamma_threshes", &
        NF90_FLOAT, nthreshes_dimid, ngamma_threshes_varid) )    
    CALL check( nf90_def_var(ncid, "thresh_light", &
        NF90_FLOAT, nscalar_dimid, nthresh_light_varid) ) 
    CALL check( nf90_def_var(ncid, "thresh_mod", &
        NF90_FLOAT, nscalar_dimid, nthresh_mod_varid) )            
    CALL check( nf90_def_var(ncid, "thresh_high", &
        NF90_FLOAT, nscalar_dimid, nthresh_high_varid) )                 
    CALL check( nf90_def_var(ncid, "climo_pop_thresholds", &
        NF90_FLOAT, nccats_minus1_dimid, nclimo_pop_thresh_varid) )    
    CALL check( nf90_def_var(ncid, "output_threshes", &
        NF90_FLOAT, nout_thresh_dimid, noutput_threshes_varid) )         
        
    ! --- End define mode. This tells netCDF we are done defining metadata.

    CALL check( nf90_enddef(ncid) )

    ! ---- write the data.  

    PRINT *,'writing data'
    CALL check( nf90_put_var(ncid, nhist_varid, closest_histogram))
    CALL check( nf90_put_var(ncid, ngsum_varid, gamma_sum))
    CALL check( nf90_put_var(ncid, nglnsum_varid, gamma_ln_sum))
    CALL check( nf90_put_var(ncid, nweights_varid, sum_of_weights))
    CALL check( nf90_put_var(ncid, nzeros_varid, nzeros))
    CALL check( nf90_put_var(ncid, npositive_varid, npositive))  
    CALL check( nf90_put_var(ncid, ngamma_sum_ensmeanzero_fclim_varid, &
        gamma_sum_ensmeanzero_fclim))
    CALL check( nf90_put_var(ncid, ngamma_ln_sum_ensmeanzero_fclim_varid, &
        gamma_ln_sum_ensmeanzero_fclim)) 
    CALL check( nf90_put_var(ncid, nzeros_fclim_varid, nzeros_fclim))
    CALL check( nf90_put_var(ncid, npositive_fclim_varid, npositive_fclim))  
    CALL check( nf90_put_var(ncid, ngamma_threshes_varid, gamma_threshes))
    CALL check( nf90_put_var(ncid, nthresh_light_varid, thresh_light))
    CALL check( nf90_put_var(ncid, nthresh_mod_varid, thresh_mod))    
    CALL check( nf90_put_var(ncid, nthresh_high_varid, thresh_high))
    CALL check( nf90_put_var(ncid, nclimo_pop_thresh_varid, &
        climo_pop_thresholds))
    CALL check( nf90_put_var(ncid, noutput_threshes_varid, &
        output_threshes))    
    CALL check(nf90_put_var(ncid, nexceed_yes_fclimPOP_varid, exceed_yes))    
    CALL check(nf90_put_var(ncid, nexceed_no_fclimPOP_varid, exceed_no))        
    
    ! ---- Close the file. This frees up any internal netCDF resources
    !      associated with the file, and flushes any buffers.
    
    CALL check( nf90_close(ncid) )
    PRINT *, "*** SUCCESS writing netcdf file "
    
ENDIF
				
RETURN
END SUBROUTINE tally_gamma_stats_full_n25

! ===================================================================

SUBROUTINE define_indices_weights(nthreshes, gamma_threshes, &
    ilower_bound, iupper_bound, weightfn)
    
INTEGER, INTENT(IN) :: nthreshes
REAL, INTENT(IN), DIMENSION(nthreshes) :: gamma_threshes
INTEGER, INTENT(OUT), DIMENSION(nthreshes) :: ilower_bound, iupper_bound
REAL*8, INTENT(OUT), DIMENSION(nthreshes, nthreshes) :: weightfn
REAL, DIMENSION(nthreshes) :: gamma_threshes_alog


gamma_threshes_alog(1) = 0.0
gamma_threshes_alog(2:nthreshes) = ALOG(gamma_threshes(2:nthreshes))

! ---- process each index

weightfn(:,:) = 0.0
weightfn(1,1) = 1.0  ! effectively no weight for 0 on other positive

ilower_bound(:) = 1
iupper_bound(:) = 1

DO ithresh = 2, nthreshes

    ! --- define indices that bound where we will generate the weights.
    
    gmin = gamma_threshes(ithresh)/2.
    gmax = gamma_threshes(ithresh)*2.
    DO it2 = 1, nthreshes
        IF (gamma_threshes(it2) .ge. gmin) THEN
            ilower_bound(ithresh) = it2
            GOTO 777
        ENDIF
    END DO
    777 CONTINUE
    
    DO it2 = nthreshes,1,-1
        IF (gamma_threshes(it2) .le. gmax) THEN
            iupper_bound(ithresh) = it2
            GOTO 888
        ENDIF
    END DO
    888 CONTINUE
    
    ! --- define the weights for indices inside bounds. 
    
    DO iw1 = ilower_bound(ithresh), iupper_bound(ithresh)
        weightfn(ithresh,iw1) = 1.0 - &
            ABS(gamma_threshes_alog(iw1)-gamma_threshes_alog(ithresh)) / 0.6931
        IF (weightfn(ithresh,iw1) .lt. 0.0) weightfn(ithresh,iw1) = 0.0
    END DO 
    
END DO

RETURN
END SUBROUTINE define_indices_weights
    
       
