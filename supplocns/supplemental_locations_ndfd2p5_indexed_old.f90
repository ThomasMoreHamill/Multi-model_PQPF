PROGRAM supplemental_locations_ndfd2p5_indexed

! purpose:  For every (i,j) in on the 2.5-km NDFD Climatology-Calibrated
!    Precipitation Analysis (CCPA; see Hou et al. 2014, available at
!    http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-11-0140.1)  grid 
!    with valid data, determine a set of supplemental data locations.
!    Base these supplemental locations on the similarity of analyzed CDFs
!    (via gamma-distribution alpha, beta, and fraction zero parameters) 
!    as well as the local terrain heights and gradients.  Generate the 
!    supplemental locations separately for every month of the year. We want 
!    grid points that will have somewhat independent data, so make 
!    sure that the set of supplemental locations are not too close together.
! 
USE netcdf

PARAMETER (nyears = 15) ! 2002-2016 here for the CCPA data
PARAMETER (nxa = 2145) ! NDFD 2.5 km grid x dimensions for expanded grid.
PARAMETER (nya = 1597) ! NDFD 2.5 km grid y dimensions for expanded grid
PARAMETER (npoints = nxa*nya)
PARAMETER (nmonths = 12) ! # months of the year
PARAMETER (nsupplemental = 50) ! max number of CCPA grid supplemental locations  
    ! to be stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 25) ! user-defined minimum separation between 
    ! analysis grid locations and supp location (in km)
PARAMETER (maxseparation = 1000) ! user-defined maximum separation between  
    ! analysis gridlocation and supplemental location (in grid pts)

PARAMETER (alpha_coeff = 0.15) ! coefficient to apply to Gamma dist alpha parameter
PARAMETER (beta_coeff = 0.15) ! coefficient to apply to Gamma dist beta parameter
PARAMETER (fz_coeff = 0.15) ! coefficient to apply to fraction zero parameter
PARAMETER (terht_coeff = 0.25) ! coefficient to apply to terrain height differences
PARAMETER (gradx_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (grady_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (distpenalty = .0001) ! coefficient to apply to distance  
PARAMETER (earth_radius_meters = 6370000.)

INTEGER, DIMENSION(2) :: mindices
INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(npoints) :: conusmask_1d  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(nxa,nya) :: maskout  ! array used to block out nearby points !
                           ! from further consideration
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: xlocation_supp ! for each forecast point,  
                           ! a list of which other forecast points
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: ylocation_supp ! have the closest climatologies 
                           ! and forecast-obs relationship
INTEGER minx
INTEGER miny
INTEGER, DIMENSION(nxa,nya) :: i_index ! i coordinate on grid
INTEGER, DIMENSION(nxa,nya) :: j_index ! j coordinate on grid

REAL, DIMENSION(nxa,nya) :: alphahat ! estimated Gamma dist shape parameter read in from file
REAL, DIMENSION(nxa,nya) :: alphahat_mean  ! local spatial mean of alpha parameter
REAL, DIMENSION(nxa,nya) :: alphahat_stddev ! local spatial std dev of alpha parameter
REAL, DIMENSION(nxa,nya) :: betahat ! estimated Gamma dist scale parameter read in from file
REAL, DIMENSION(nxa,nya) :: betahat_mean ! local spatial mean of beta parameter
REAL, DIMENSION(nxa,nya) :: betahat_stddev ! local spatial std dev of beta parameter

REAL, DIMENSION(nxa,nya) :: difference ! total penalty difference between (ixa,jya) 
!  and potential supplemental (ixn,jyn)
!REAL, DIMENSION(nxa,nya) :: difference_alpha ! Gamma dist alpha parameter penalty
!REAL, DIMENSION(nxa,nya) :: difference_beta ! Gamma dist beta parameter penalty
!REAL, DIMENSION(nxa,nya) :: difference_fz ! fraction zero penalty
!REAL, DIMENSION(nxa,nya) :: difference_terht  ! terrain height penalty
!REAL, DIMENSION(nxa,nya) :: difference_gradx  ! E-W terrain height gradient penalty
!REAL, DIMENSION(nxa,nya) :: difference_grady  ! N-S terrain height gradient penalty
!REAL, DIMENSION(nxa,nya) :: difference_dist ! physical distance penalty
REAL, DIMENSION(nxa,nya) :: difference_init_large
REAL, DIMENSION(nxa,nya) :: fraction_zero ! estimated fraction zero parameter from file
REAL, DIMENSION(nxa,nya) :: fz_mean ! local spatial mean of fraction zero parameter
REAL, DIMENSION(nxa,nya) :: fz_stddev ! local spatial std dev of fraction zero parameter
REAL, DIMENSION(nxa,nya) :: gradx_mean ! local spatial mean of E-W terrain gradient
REAL, DIMENSION(nxa,nya) :: gradx_stddev ! local spatial std dev of E-W terrain gradient
REAL, DIMENSION(nxa,nya) :: grady_mean ! local spatial mean of N-S terrain gradient
REAL, DIMENSION(nxa,nya) :: grady_stddev ! local spatial std dev of N-S terrain gradient
REAL, DIMENSION(nxa,nya) :: lonsa ! analysis grid longitudes (negative is W lon)
REAL, DIMENSION(nxa,nya) :: latsa ! analysis grid latitudes in degrees
REAL, DIMENSION(nxa,nya,nsupplemental) :: penalty  ! penalty function 
    ! associated with the supp location
REAL, DIMENSION(nxa,nya) :: terrain ! terrain elevation
REAL, DIMENSION(nxa,nya) :: terrain_gradx ! terrain E-W gradient
REAL, DIMENSION(nxa,nya) :: terrain_grady ! terrain N-S gradient
REAL, DIMENSION(nxa,nya) :: ter_mean
REAL, DIMENSION(nxa,nya) :: ter_stddev

INTEGER*2, DIMENSION(nxa, nya) :: iconsider

REAL time1, time2

REAL diffmin

CHARACTER*120 infile, outfile, data_directory
CHARACTER*3, DIMENSION(12) :: cmonths
CHARACTER*3 csupp
CHARACTER*2 cleadb, cleade

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun', &
             'Jul','Aug','Sep','Oct','Nov','Dec'/
             
DATA data_directory /'/Projects/Reforecast2/netcdf/NationalBlend/'/

CALL getarg(1, cleadb) ! 00 for 00 UTC, or 12
CALL getarg(2, cleade) ! 12 for 12 UTC, or 00

IF (nsupplemental .ge. 100) THEN
   WRITE (csupp,'(i3)') nsupplemental
ELSE
   WRITE (csupp,'(i2)') nsupplemental
   csupp = '0' // csupp(1:2)
ENDIF
print *,'csupp = ',csupp

! ---- for each month, load the quantile 
!      information and rank correlation previously
!      computed, and then determine analogs

xlocation = -99  ! initialize all locations to missing value.
ylocation = -99
xlocation_supp = -99
ylocation_supp = -99

! ---- initialize 

penalty = 0.
penalty_alpha = 0.
penalty_beta = 0.
penalty_fz = 0.
penalty_ter = 0.
penalty_gradx = 0.
penalty_grady = 0.
penalty_dist = 0.
DO jya = 1, nya
    DO ixa = 1, nxa
        i_index(ixa,jya) = ixa
        j_index(ixa,jya) = jya
    END DO
END DO
difference_init_large = 999999999.

! ---- read in terrain facet information for smoothing 
!      at short, intermediate, and larger scales

infile = TRIM(data_directory) // 'blend.precip_const.2p5.nc'
PRINT *, 'calling read_terrain_heights_2p5'
PRINT *, TRIM(infile)
CALL read_terrain_heights_2p5(nxa, nya, infile, terrain, &
    lonsa, latsa, conusmask)
npoints_in_conus = SUM(INT(conusmask))
ktr = 0
DO jya = 1, nya
    DO ixa = 1, nxa
        ktr = ktr+1
        conusmask_1d(ktr) = conusmask(ixa,jya)
    END DO
END DO     

! ---- calculate terrain gradients

PRINT *, 'calling calculate_terrain_gradients'
CALL calculate_terrain_gradients(nxa, nya, earth_radius_meters, &
    terrain, latsa, terrain_gradx, terrain_grady)
    
! ---- loop over months

!DO imonth = 4, 4
DO imonth = 1,12

    ! ---- read in the parameters of the Gamma distribution for 
    !      climatology, and fraction zero.

    infile = TRIM(data_directory) // &
        'climatology_gamma_parameters_ndfd2p5_' // &
        cleadb // '_to_' // cleade // '_' // &
        cmonths(imonth) //'.nc'
    PRINT *, TRIM(infile)
    CALL read_climatology_parameters_ndfd2p5(nxa, nya, infile, &
        fraction_zero, alphahat, betahat)  
    
             
    ! ---- loop thru grid points and precompute statistics for mean, std dev
    !      of alpha, beta, fraction zero, and terrain.  These will be used to
    !      scale the various pieces of information used to calculate 
    !      supplemental location penalties.

    PRINT *,'pre-processing to determine local standard devs of CDF ',&
        '& terrain parameters'
    DO ixa = 1, nxa
        IF (MOD(ixa,100) .eq. 0) PRINT *,'processing  ixa = ',ixa,' of ',nxa

        ! --- only bother searching in a box +/- 10 grid points to 
        !     limit computations.

        ixmin = MAX(1,ixa-10)
        ixmax = MIN(nxa,ixa+10)
        DO jya = 1, nya
            jymin = MAX(1,jya-10)
            jymax = MIN(nya,jya+10)
            sum_fz2 = 0.0
            sum_fz  = 0.0
            sum_alpha2 = 0.0
            sum_alpha  = 0.0
            sum_beta2 = 0.0
            sum_beta = 0.0
            sum_gradx = 0.0
            sum_gradx2 = 0.0
            sum_grady = 0.0
            sum_grady2 = 0.0
            sum_terht = 0.0
            sum_terht2 = 0.0
            nsamps  = 0
            IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point has valid data
                fz_here = fraction_zero(ixa,jya)
                alpha_here = alphahat(ixa,jya)
                beta_here = betahat(ixa,jya)
                gradx_here = terrain_gradx(ixa,jya)
                grady_here = terrain_grady(ixa,jya)
                ter_here = terrain(ixa,jya)
                DO ix2 = ixmin, ixmax
                    DO jy2 = jymin, jymax
                        IF (conusmask(ix2,jy2).gt. 0) THEN
                            fz_there = fraction_zero(ix2,jy2)
                            alpha_there = alphahat(ix2,jy2)
                            beta_there = betahat(ix2,jy2)
                            gradx_there = terrain_gradx(ix2,jy2)
                            grady_there = terrain_grady(ix2,jy2)
                            ter_there = terrain(ix2,jy2)
					 
                            ! --- quantify the difference between forecast distributions at 
                            !     the two grid points
					 
                            sum_fz2 = sum_fz2 + fz_there**2
                            sum_fz  = sum_fz + fz_there
                            sum_alpha2 = sum_alpha2 + alpha_there**2
                            sum_alpha  = sum_alpha + alpha_there
                            sum_beta2 = sum_beta2 + beta_there**2
                            sum_beta = sum_beta + beta_there
                            sum_gradx = sum_gradx + gradx_there
                            sum_gradx2 = sum_gradx2 + gradx_there**2
                            sum_grady = sum_grady + grady_there
                            sum_grady2 = sum_grady2 + grady_there**2
                            sum_terht = sum_terht + ter_there
                            sum_terht2 = sum_terht2 + ter_there**2
                            nsamps  = nsamps  + 1
                        ENDIF
                    END DO
                END DO

                ! --- now calculate std deviation by shortcut (Wilks 2006, eq. 3.20)

                alphahat_mean(ixa,jya) = sum_alpha / REAL(nsamps)
                alphahat_stddev(ixa,jya) = SQRT((sum_alpha2-REAL(nsamps)*alphahat_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
                betahat_mean(ixa,jya) = sum_beta / REAL(nsamps)
                betahat_stddev(ixa,jya) = SQRT((sum_beta2-REAL(nsamps)*betahat_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
                fz_mean(ixa,jya) = sum_fz / REAL(nsamps)
                fz_stddev(ixa,jya) = SQRT((sum_fz2-REAL(nsamps)*fz_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
                gradx_mean(ixa,jya) = sum_gradx / REAL(nsamps)
                gradx_stddev(ixa,jya) = SQRT((sum_gradx2-REAL(nsamps)*gradx_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
                grady_mean(ixa,jya) = sum_grady / REAL(nsamps)
                grady_stddev(ixa,jya) = SQRT((sum_grady2-REAL(nsamps)*grady_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
                ter_mean(ixa,jya) = sum_terht  / REAL(nsamps)
                ter_stddev(ixa,jya) = SQRT((sum_terht2-REAL(nsamps)*ter_mean(ixa,jya)**2) / &
                    REAL(nsamps-1))
            ELSE
                alphahat_mean(ixa,jya) = -99.99
                alphahat_stddev(ixa,jya) = -99.99
                betahat_mean(ixa,jya) = -99.99
                betahat_stddev(ixa,jya) = -99.99
                fz_mean(ixa,jya) = -99.99
                fz_stddev(ixa,jya) = -99.99
                gradx_mean(ixa,jya) = -99.99
                gradx_stddev(ixa,jya) = -99.99
                grady_mean(ixa,jya) = -99.99
                grady_stddev(ixa,jya) = -99.99
                ter_mean(ixa,jya) = -99.99
                ter_stddev(ixa,jya) = -99.99
            ENDIF

        END DO ! jya
    END DO    ! ixa
    
    !PRINT *,'i    lon    alpham alphastd  betam betastd fzmean fzstd  gradxmean gradxstd  termean    terstd'
    !DO ixa = 1, nxa
    !    PRINT 200, ixa, lonsa(ixa,nya/2), alphahat_mean(ixa,nya/2), alphahat_stddev(ixa,nya/2), &
    !        betahat_mean(ixa,nya/2), betahat_stddev(ixa,nya/2), &
    !        fz_mean(ixa,nya/2), fz_stddev(ixa,nya/2), gradx_mean(ixa,nya/2), &
    !        gradx_stddev(ixa,nya/2), ter_mean(ixa,nya/2), ter_stddev(ixa,nya/2)
    !    200 FORMAT (i3,1x,f8.2,8(1x,f7.3),2(1x,f11.1))
    !END DO

    stdmax = MAXVAL(ter_stddev)
    stdmax_sqrt = stdmax**0.5
    gradxmax = MAXVAL(gradx_stddev)
    gradxmax_sqrt = SQRT(gradxmax)
    gradymax = MAXVAL(grady_stddev)
    gradymax_sqrt = SQRT(gradymax)
    
    ! =======================================================================================
    ! ---- find locations of supplemental locations that are the best fit for each grid point
    ! =======================================================================================
    
    CALL cpu_time(time1)
    PRINT *, 'finding supplemental locations'
    nproc = 0

    ktrpoint = 0
    DO ixa = 1, nxa
    !DO ixa = nxa/3, nxa/3
        CALL cpu_time(time2)
        IF (MOD(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa, &
	  	    ' time elapsed, delta: ',time2, time2-time1,'nproc = ', &
            nproc,' total pts = ',SUM(REAL(conusmask))
        DO jya = 1, nya
        !DO jya = nya/2, nya/2

            IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point inside area with CCPA data
                
                maskout(:,:) = 1-conusmask(:,:)
                difference(:,:) = difference_init_large(:,:) 
    	        !difference_alpha(:,:) = difference_init_large(:,:) 
                !difference_beta(:,:) = difference_init_large(:,:) 
                !difference_fz(:,:) = difference_init_large(:,:) 
                !difference_terht(:,:) = difference_init_large(:,:) 
                !difference_gradx(:,:) = difference_init_large(:,:) 
                !difference_grady(:,:) = difference_init_large(:,:) 
                !difference_dist(:,:) = difference_init_large(:,:)
                
                ! ---- rather than examining all points for their closeness to the 
                !      current grid point's alpha, beta, fraction zero, and terrain
                !      values, let's thin down the number of grid points to examine
                !      to a smaller set that have similar values
                
                ktrpoint = ixa + (jya-1)*nxa
                
                iconsider = conusmask
                alphahat_min = MAX(0.001,alphahat(ixa,jya) - alphahat_stddev(ixa,jya))
                alphahat_max = alphahat(ixa,jya) + alphahat_stddev(ixa,jya)
                betahat_min = MAX(0.001,betahat(ixa,jya) - betahat_stddev(ixa,jya))
                betahat_max = betahat(ixa,jya) + betahat_stddev(ixa,jya)
                fraction_zero_min = MAX(0.001,fraction_zero(ixa,jya) - fz_stddev(ixa,jya))
                fraction_zero_max = fraction_zero(ixa,jya) + fz_stddev(ixa,jya)

                nproc = nproc+1
                imin = max(1,ixa-maxseparation)
                imax = min(nxa,ixa+maxseparation)
                jmin = max(1,jya-maxseparation)
                jmax = min(nya,jya+maxseparation)
                
                fz_here = fraction_zero(ixa,jya)
                alpha_here = alphahat(ixa,jya)
                beta_here = betahat(ixa,jya)
                gradx_here = terrain_gradx(ixa,jya)
                grady_here = terrain_grady(ixa,jya)
                ter_here = terrain(ixa,jya)

                terht_coeff2 = terht_coeff * SQRT(ter_stddev(ixa,jya)) / stdmax_sqrt
                gradx_coeff2 = gradx_coeff * SQRT(gradx_stddev(ixa,jya)) / gradxmax_sqrt
                grady_coeff2 = grady_coeff * SQRT(grady_stddev(ixa,jya)) / gradymax_sqrt

                nconsider = 0
                DO ix2 = imin, imax
                    DO jy2 = jmin, jmax
                        
                        IF (conusmask(ix2,jy2) .gt. 0) THEN
                            IF (alphahat(ix2,jy2) .ge. alphahat_min .and. &
                            alphahat(ix2,jy2) .le. alphahat_max .and. & 
                            betahat(ix2,jy2) .ge. betahat_min .and. &
                            betahat(ix2,jy2) .le. betahat_max .and. &
                            fraction_zero(ix2,jy2) .ge. fraction_zero_min .and. &
                            fraction_zero(ix2,jy2) .le. fraction_zero_max) THEN

                                ! ---- Also, only consider this grid point if inside area with data, 
                                !      and is less than max separation, computed in 1/8-degree 
                                !      grid points.  Since there is less trustworthy data outside 
                                !      CONUS (lots of problems over ocean, over Canada at this  
                                !      point, supplement data only with points from inside CONUS

                                CALL haversine(latsa(ixa,jya), lonsa(ixa,jya), &
                                    latsa(ix2,jy2), lonsa(ix2,jy2), dist)

                                IF (dist .le. maxseparation) THEN
                            
                                    nconsider = nconsider + 1
                                    fz_there = fraction_zero(ix2,jy2)
                                    alpha_there = alphahat(ix2,jy2)
                                    beta_there = betahat(ix2,jy2)
                                    gradx_there = terrain_gradx(ix2,jy2)
                                    grady_there = terrain_grady(ix2,jy2)
                                    ter_there = terrain(ix2,jy2)

                                    ! --- quantify the difference between fraction zero, alpha, beta 

                                    fz_diff = ABS(fz_here - fz_there) / fz_stddev(ixa,jya) 
                                    alpha_diff = ABS(alpha_here - alpha_there) / alphahat_stddev(ixa,jya)
                                    beta_diff = ABS(beta_here - beta_there) / betahat_stddev(ixa,jya)
                                    gradx_diff = ABS(gradx_here - gradx_there) / gradx_stddev(ixa,jya)
                                    grady_diff = ABS(grady_here - grady_there) / grady_stddev(ixa,jya)
                                    terr_diff = ABS(ter_here - ter_there) / ter_stddev(ixa,jya)

                                    ! --- calculate differences between this potential supplemental location
                                    !     and the grid point of interest for each of the various factors 
                                    !     we choose to weight 
                            
                                    !difference_alpha(ix2,jy2) = alpha_coeff*alpha_diff 
                                    !difference_beta(ix2,jy2) = beta_coeff*beta_diff
                                    !difference_fz(ix2,jy2) = fz_coeff*fz_diff
                                    !difference_terht(ix2,jy2) = terht_coeff2*terr_diff
                                    !difference_gradx(ix2,jy2) = gradx_coeff2*gradx_diff
                                    !difference_grady(ix2,jy2) = grady_coeff2*grady_diff
                                    !difference_dist(ix2,jy2) = dist*distpenalty
                                    !difference(ix2,jy2) = difference_alpha(ix2,jy2) + &
                                    !    difference_beta(ix2,jy2) + difference_fz(ix2,jy2) + &
                                    !    difference_terht(ix2,jy2) + difference_gradx(ix2,jy2) + &
                                    !    difference_grady(ix2,jy2) + difference_dist(ix2,jy2)
                                    difference(ix2,jy2) = alpha_coeff*alpha_diff  + &
                                        beta_coeff*beta_diff + fz_coeff*fz_diff + &
                                        terht_coeff2*terr_diff + gradx_coeff2*gradx_diff + &
                                        grady_coeff2*grady_diff + dist*distpenalty                    

                                ENDIF ! dist .le. maxseparation
                            ELSE
                                maskout(ix2,jy2) = 1 ! don't consider this grid point
                            ENDIF ! alphahat(ix2,jy2) .ge. alphahat_min, etc
                        ENDIF ! conusmask(ix2,jy2) 
                    END DO ! jy2
                END DO ! ix2
                !PRINT *,'considered ', nconsider,' out of ',(imax-imin+1)*(jmax-jmin+1)

                DO isupp = 1, nsupplemental  ! number of supplemental locations

                    IF (isupp .eq. 1) THEN ! the first supplemental location
                                     !  is simply that orginal grid point
                                     
                        xlocation_supp(ixa,jya,isupp) = ixa
                        ylocation_supp(ixa,jya,isupp) = jya
                        
                        ! ---- don't consider points right around grid pt of interest,
                        !      too strong a correlation (want quasi-independent samples)

                        diffmin = 0.
                        CALL mask_around_thisgridpt(ixa, jya, nxa, nya, &
                            minseparation, maskout)

                    ELSE

				        ! ---- now find the analysis grid point with the next 
				        !      closest similarity. Set this as supplemental location 
				        !      point and then mask around this to eliminate nearby
				        !      points from any future consideration

                        mindices = MINLOC(difference)
                        minx = mindices(1)
                        miny = mindices(2)
                        
                        !minx = -99
                        !miny = -99
                        !diffmin = 999999.
                        !DO ix2 = imin, imax   !1, nxa
                        !    DO jy2 = jmin, jmax  !1, nya
                        !        IF (difference(ix2,jy2) .lt. diffmin .and. &
                        !        maskout(ix2,jy2) .eq. 0) THEN
                        !            diffmin = difference(ix2,jy2)
                        !            minx = ix2
                        !            miny = jy2
                        !        END IF
                        !    END DO ! jy2
                        !END DO  ! ix2
               
                        ! ---- finally, define the supplemental 
				        !      location on the forecast grid to be the 
                        !      grid point that had this closest fit as defined above.
                  
                        xlocation_supp(ixa,jya,isupp) = minx
                        ylocation_supp(ixa,jya,isupp) = miny

                        ! ---- make sure no other grid points very close to  
                        !      this new supplemental location are considered 
				        !      as another potential supplemental location.

                        CALL mask_around_thisgridpt(minx, miny, nxa, nya, &
                            minseparation, maskout)

                    ENDIF ! isupp > 1
                    
                    !PRINT *,'isupp, x, y = ',isupp, xlocation_supp(ixa,jya,isupp), ylocation_supp(ixa,jya,isupp)

                END DO ! isupp
            ENDIF   ! conusmask
        END DO  ! jya
    END DO ! ixa

    ! --- save supplemental locations and penalty informationfor this month

    outfile = TRIM(data_directory) // 'supplemental_locations_CONUS_ndfd2p5_'//&
        cleadb //'_to_'//cleade//'_'//cmonths(imonth) // '.nc'
    CALL write_supp_locns_to_netcdf(outfile, nxa, nya, nsupplemental, &
        xlocation_supp, ylocation_supp, conusmask, lonsa, latsa)

END DO  ! imonth
PRINT *,'Done'
CONTAINS                        

! =============================================================================

SUBROUTINE mask_around_thisgridpt(minx, miny, nx, ny, minseparation, maskout)
! set maskout array around grid point of interest to 1 
!    to indicate nonconsideration of this in the future.
INTEGER, INTENT(IN) :: minx, miny, nx, ny, minseparation
INTEGER*2, INTENT(INOUT), DIMENSION(nx,ny) :: maskout
imin_m = MAX(1,minx-minseparation)
imax_m = MIN(nx,minx+minseparation)
jmin_m = MAX(1,miny-minseparation)
jmax_m = MIN(ny,miny+minseparation)
!$OMP PARALLEL DO DEFAULT(SHARED) COLLAPSE(2) PRIVATE(i,j,dist2,dist)
DO i = imin_m, imax_m
   DO j = jmin_m, jmax_m
      dist2 = REAL((minx-i)**2 + (miny-j)**2)
      dist = SQRT(dist2)
      IF (dist .le. minseparation) maskout(i,j) = 1
   END DO
END DO
!$OMP END PARALLEL DO
RETURN
END SUBROUTINE mask_around_thisgridpt

! ==============================================================================

SUBROUTINE haversine (deglat1, deglon1, deglat2, deglon2, dist)
! great circle distance -- adapted from Matlab 
real,intent(in) :: deglat1,deglon1,deglat2,deglon2
real :: a,c2,dist,dlat,dlon,lat1,lat2
real,parameter :: radius = 6372.8 
 
CALL to_radian(deglat2-deglat1,dlat)
CALL to_radian(deglon2-deglon1,dlon)
CALL to_radian(deglat1,lat1)
CALL to_radian(deglat2,lat2)
a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
c2 = 2*asin(sqrt(a))
dist = radius*c2
END SUBROUTINE haversine

! ======================================================================
 
SUBROUTINE to_radian(degree, rad)
REAL, INTENT(IN) :: degree
REAL, PARAMETER :: deg_to_rad = atan(1.0)/45.
    ! exploit intrinsic atan to generate pi/180 runtime constant
REAL :: rad

rad = degree*deg_to_rad
END SUBROUTINE to_radian

! ======================================================================

SUBROUTINE write_supp_locns_to_netcdf(outfile, nxa, nya, nsupplemental, &
    xlocation_supp, ylocation_supp, conusmask, lonsa, latsa)

USE netcdf  ! use Unidata netCDF fortran interface library

CHARACTER*(*), INTENT(IN) :: outfile
INTEGER, INTENT(IN) :: nxa, nya, nsupplemental
INTEGER, INTENT(IN), DIMENSION(nxa, nya, nsupplemental) :: &
    xlocation_supp, ylocation_supp
INTEGER*2, DIMENSION(nxa,nya), INTENT(IN) :: conusmask
REAL, DIMENSION(nxa,nya), INTENT(IN) :: lonsa, latsa
INTEGER :: dimid_3d(3), dimid_2d(2)

print *,'writing to ',TRIM(outfile)
CALL check( nf90_create(TRIM(outfile), NF90_CLOBBER, ncid) )

! ---- Define the array dimensions. NetCDF will hand back an ID for each.

PRINT *,'array dimensions'
CALL check( nf90_def_dim(ncid, "nxa", nxa, nxa_dimid) )
CALL check( nf90_def_dim(ncid, "nya", nya, nya_dimid) )
CALL check( nf90_def_dim(ncid, "nsupp", nsupplemental, nsupp_dimid) )
dimid_2d =  (/ nxa_dimid, nya_dimid /)
dimid_3d =  (/ nxa_dimid, nya_dimid, nsupp_dimid /)

! ---- Define the variables and associated IDs

PRINT *,'defining variables'
CALL check( nf90_def_var(ncid, "conusmask", &
    NF90_SHORT, dimid_2d, nconusmask_varid) )
CALL check( nf90_def_var(ncid, "lons", &
    NF90_FLOAT, dimid_2d, nlons_varid) )    
CALL check( nf90_def_var(ncid, "lats", &
    NF90_FLOAT, dimid_2d, nlats_varid) ) 
CALL check( nf90_def_var(ncid, "xlocation_supp", &
    NF90_INT, dimid_3d, nxsupp_varid) )
CALL check( nf90_def_var(ncid, "ylocation_supp", &
    NF90_INT, dimid_3d, nysupp_varid) )  
CALL check( nf90_enddef(ncid) )

! ---- write the data.

PRINT *,'writing data'
CALL check( nf90_put_var(ncid, nconusmask_varid, conusmask))
CALL check( nf90_put_var(ncid, nlons_varid, lonsa))
CALL check( nf90_put_var(ncid, nlats_varid, latsa))
print *,'writing supplemental location data'
CALL check( nf90_put_var(ncid, nxsupp_varid, xlocation_supp))
CALL check( nf90_put_var(ncid, nysupp_varid, ylocation_supp))

! ---- Close the file. This frees up any internal netCDF resources
!      associated with the file, and flushes any buffers.

CALL check( nf90_close(ncid) )
PRINT *, "*** SUCCESS writing netcdf file "

RETURN
END SUBROUTINE write_supp_locns_to_netcdf


! ==================================================================

END PROGRAM supplemental_locations_ndfd2p5_indexed