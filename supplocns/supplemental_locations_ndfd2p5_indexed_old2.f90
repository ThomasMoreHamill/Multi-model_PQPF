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
PARAMETER (nxa = 2345) ! NDFD 2.5 km grid x dimensions for expanded grid.
PARAMETER (nya = 1597) ! NDFD 2.5 km grid y dimensions for expanded grid
PARAMETER (npoints = nxa*nya)
PARAMETER (nmonths = 12) ! # months of the year
PARAMETER (nsupplemental = 50) ! max number of CCPA grid supplemental locations  
    ! to be stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 25) ! user-defined minimum separation between 
    ! analysis grid locations and supp location (in grid pts)
!PARAMETER (maxseparation = 1000) ! user-defined maximum separation between  
!    ! analysis gridlocation and supplemental location (in km)
PARAMETER (maxseparation = 700) ! user-defined maximum separation between  
    ! analysis gridlocation and supplemental location (in grid pts)


PARAMETER (alpha_coeff = 0.15) ! coefficient to apply to Gamma dist alpha parameter
PARAMETER (beta_coeff = 0.15) ! coefficient to apply to Gamma dist beta parameter
PARAMETER (fz_coeff = 0.15) ! coefficient to apply to fraction zero parameter
PARAMETER (terht_coeff = 0.25) ! coefficient to apply to terrain height differences
PARAMETER (gradx_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (grady_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (distpenalty = .001) ! coefficient to apply to distance  
PARAMETER (earth_radius_meters = 6370000.)

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(nxa,nya) :: maskout  ! array used to block out nearby points !
                           ! from further consideration
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: xlocation_supp ! for each forecast point,  
                           ! a list of which other forecast points
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: ylocation_supp ! have the closest climatologies 
                           ! and forecast-obs relationship

LOGICAL, DIMENSION(nxa,nya) :: test_diffs_here

INTEGER mindex(1)
INTEGER, DIMENSION(nxa,nya) :: i_index ! i coordinate on grid
INTEGER, DIMENSION(nxa,nya) :: j_index ! j coordinate on grid

REAL, DIMENSION(nxa,nya) :: alphahat ! estimated Gamma dist shape parameter read in from file
REAL*8, DIMENSION(nxa,nya) :: alphahat_mean  ! local spatial mean of alpha parameter
REAL*8, DIMENSION(nxa,nya) :: alphahat_stddev ! local spatial std dev of alpha parameter
REAL, DIMENSION(nxa,nya) :: betahat ! estimated Gamma dist scale parameter read in from file
REAL*8, DIMENSION(nxa,nya) :: betahat_mean ! local spatial mean of beta parameter
REAL*8, DIMENSION(nxa,nya) :: betahat_stddev ! local spatial std dev of beta parameter

REAL, ALLOCATABLE, DIMENSION(:) :: difference ! total penalty difference between (ixa,jya) 
REAL, DIMENSION(nxa,nya) :: dist
!  and potential supplemental (ixn,jyn)
REAL, DIMENSION(nxa,nya) :: fraction_zero ! estimated fraction zero parameter from file
REAL*8, DIMENSION(nxa,nya) :: fz_mean ! local spatial mean of fraction zero parameter
REAL*8, DIMENSION(nxa,nya) :: fz_stddev ! local spatial std dev of fraction zero parameter
REAL*8, DIMENSION(nxa,nya) :: gradx_mean ! local spatial mean of E-W terrain gradient
REAL*8, DIMENSION(nxa,nya) :: gradx_stddev ! local spatial std dev of E-W terrain gradient
REAL*8, DIMENSION(nxa,nya) :: grady_mean ! local spatial mean of N-S terrain gradient
REAL*8, DIMENSION(nxa,nya) :: grady_stddev ! local spatial std dev of N-S terrain gradient
REAL, DIMENSION(nxa,nya) :: lonsa ! analysis grid longitudes (negative is W lon)
REAL, DIMENSION(nxa,nya) :: latsa ! analysis grid latitudes in degrees
REAL, DIMENSION(nxa,nya,nsupplemental) :: penalty  ! penalty function 
    ! associated with the supp location
REAL, DIMENSION(nxa,nya) :: terrain ! terrain elevation
REAL, DIMENSION(nxa,nya) :: terrain_gradx ! terrain E-W gradient
REAL, DIMENSION(nxa,nya) :: terrain_grady ! terrain N-S gradient
REAL*8, DIMENSION(nxa,nya) :: ter_mean
REAL*8, DIMENSION(nxa,nya) :: ter_stddev

REAL*8 :: sum_fz2 
REAL*8 :: sum_fz
REAL*8 :: sum_alpha2
REAL*8 :: sum_alpha 
REAL*8 :: sum_beta2 
REAL*8 :: sum_beta
REAL*8 :: sum_gradx 
REAL*8 :: sum_gradx2
REAL*8 :: sum_grady 
REAL*8 :: sum_grady2 
REAL*8 :: sum_terht 
REAL*8 :: sum_terht2

REAL*8 :: fz_here 
REAL*8 :: alpha_here
REAL*8 :: beta_here 
REAL*8 :: gradx_here
REAL*8 :: grady_here
REAL*8 :: ter_here

REAL*8 :: fz_there 
REAL*8 :: alpha_there
REAL*8 :: beta_there 
REAL*8 :: gradx_there
REAL*8 :: grady_there
REAL*8 :: ter_there

INTEGER*2, DIMENSION(nxa, nya) :: iconsider
INTEGER, ALLOCATABLE, DIMENSION(:) :: i_indices, j_indices


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
DO jya = 1, nya
    DO ixa = 1, nxa
        i_index(ixa,jya) = ixa
        j_index(ixa,jya) = jya
    END DO
END DO

! ---- read in terrain facet information for smoothing 
!      at short, intermediate, and larger scales

infile = TRIM(data_directory) // 'blend.precip_const.2p5.nc'
PRINT *, 'calling read_terrain_heights_2p5'
PRINT *, TRIM(infile)
CALL read_terrain_heights_2p5(nxa, nya, infile, terrain, &
    lonsa, latsa)
npoints_in_conus = SUM(INT(conusmask))    

! ---- calculate terrain gradients

PRINT *, 'calling calculate_terrain_gradients'
CALL calculate_terrain_gradients(nxa, nya, earth_radius_meters, &
    terrain, latsa, terrain_gradx, terrain_grady)
    
! ---- loop over months

!DO imonth = 4, 4
DO imonth = 2,12

    ! ---- read in the parameters of the Gamma distribution for 
    !      climatology, and fraction zero.

    infile = TRIM(data_directory) // &
        'climatology_gamma_parameters_ndfd2p5_' // &
        cleadb // '_to_' // cleade // '_' // &
        cmonths(imonth) //'.nc'
    PRINT *, TRIM(infile)
    CALL read_climatology_parameters_ndfd2p5(nxa, nya, infile, &
        fraction_zero, alphahat, betahat, conusmask)  
    
    PRINT *,'1336  62 conusmask(ixa,jya), fraction_zero(ixa,jya) = ', conusmask(1336,62), fraction_zero(1336,62)
 
             
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
                            !IF (ixa .eq. 1336 .and. jya .eq. 62) THEN
                            !    PRINT *,'ix2, jy2, fraction_zero(ix2,jy2) = ', ix2, jy2, fraction_zero(ix2,jy2)
                            !ENDIF
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
                !IF (ixa .eq. 1336 .and. jya .eq. 62) THEN
                !    PRINT *,'ixa, jya, nsamps, sum_fz2, fz_mean, fz_stddev = ',&
                !        ixa, jya,  nsamps, sum_fz2, fz_mean(ixa,jya), fz_stddev(ixa,jya)
                !ENDIF    
                    
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
    
    !stop
    
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
    !DO ixa = 1, nxa
    DO ixa = 470,470
        CALL cpu_time(time2)
        IF (MOD(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa, &
	  	    ' time elapsed, delta: ',time2, time2-time1,'nproc = ', &
            nproc,' total pts = ',SUM(REAL(conusmask))
        !DO jya = 1, nya
        DO jya = 936, 936

            IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point inside area with CCPA data
                
                PRINT *,'processing ixa, jya = ', ixa, jya
                PRINT *,'lon, lat = ', lonsa(470,936), latsa(470,936)
                maskout(:,:) = 1-conusmask(:,:)
                
                
                
! where I left off.   I think I found a grid point with very unusual alpha, beta
! parameters, but where the underlying cdf isn't really that odd.   I probably need
! to fall back to evaluation of CDFs at predefined quantiles or precip amounts.
! perhaps this can be made more efficient by integer-izing the calculations.                
                
                
                ! ---- rather than examining all points for their closeness to the 
                !      current grid point's alpha, beta, fraction zero, and terrain
                !      values, let's thin down the number of grid points to examine
                !      to a smaller set that have similar values
                
                ktrpoint = ixa + (jya-1)*nxa
                
                iconsider = conusmask
                alphahat_min = MAX(0.00001,alphahat(ixa,jya) - 1.5*alphahat_stddev(ixa,jya))
                alphahat_max = alphahat(ixa,jya) + 1.5*alphahat_stddev(ixa,jya)
                betahat_min = MAX(0.00001,betahat(ixa,jya) - 1.5*betahat_stddev(ixa,jya))
                betahat_max = betahat(ixa,jya) + 1.5*betahat_stddev(ixa,jya)
                fraction_zero_min = MAX(0.00001,fraction_zero(ixa,jya) - 1.5*fz_stddev(ixa,jya))
                fraction_zero_max = MIN(fraction_zero(ixa,jya) + 1.5*fz_stddev(ixa,jya), 1.0)
                PRINT *,'alphahat_min, alphahat, alphahat_max = ', alphahat_min, alphahat(ixa,jya), alphahat_max
                PRINT *,'betahat_min, betahat, betahat_max = ', betahat_min, betahat(ixa,jya), betahat_max
                PRINT *,'fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max =  ',&
                    fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max
                
                nproc = nproc+1
                imin = max(1,ixa-maxseparation)
                imax = min(nxa,ixa+maxseparation)
                jmin = max(1,jya-maxseparation)
                jmax = min(nya,jya+maxseparation)
                PRINT *,'imin, ixa, imax = ', imin, ixa, imax
                PRINT *,'jmin, jya, jmax = ', jmin, jya, jmax
                
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
                dist(:,:) = 999999999.
                !CALL cpu_time(time2)
                !PRINT *,'before distance calculations ', time2
                DO ix2 = imin, imax
                    DO jy2 = jmin, jmax
                        d2 = (ixa-ix2)**2 + (jya-jy2)**2
                        dist(ix2,jy2)= SQRT(d2)
                        !CALL haversine(latsa(ixa,jya), lonsa(ixa,jya), &
                        !    latsa(ix2,jy2), lonsa(ix2,jy2), dist(ix2,jy2))
                    END DO
                END DO
                !CALL cpu_time(time2)
                !PRINT *,'after distance calculations ', time2
                        
                        
                PRINT *,'alphahat_min, alphahat, alphahat_max = ', alphahat_min, alphahat(ixa,jya), alphahat_max
                PRINT *,'betahat_min, betahat, betahat_max = ', betahat_min, betahat(ixa,jya), betahat_max
                PRINT *,'fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max =  ',&
                    fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max
                    
                test_diffs_here = conusmask > 0 .and. dist < maxseparation .and. &
                alphahat >= alphahat_min .and. alphahat <= alphahat_max .and. &
                betahat >= betahat_min .and. betahat <= betahat_max .and. &
                fraction_zero >= fraction_zero_min .and. &
                fraction_zero <= fraction_zero_max 
                numtrue = count(test_diffs_here)
                PRINT *,'test_diffs_here(470, 936) = ', test_diffs_here(470, 936) 
                PRINT *,'numtrue = ',numtrue
                ALLOCATE(difference(numtrue), i_indices(numtrue), j_indices(numtrue))
                
                
                !CALL cpu_time(time2)
                !PRINT *,'after logical test ', time2
                
                ktrtrue = 0
                DO ix2 = imin, imax
                    DO jy2 = jmin, jmax
                        IF (test_diffs_here(ix2,jy2)) THEN 
                            
                            ktrtrue = ktrtrue + 1
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
                            
                            difference(ktrtrue) = alpha_coeff*alpha_diff  + &
                                beta_coeff*beta_diff + fz_coeff*fz_diff + &
                                terht_coeff2*terr_diff + gradx_coeff2*gradx_diff + &
                                grady_coeff2*grady_diff + dist(ix2,jy2)*distpenalty   
                            i_indices(ktrtrue) = ix2
                            j_indices(ktrtrue) = jy2                 
                            
                        ENDIF ! test_diffs_here
                    END DO ! jy2
                END DO ! ix2
                print *,'ktrtrue = ', ktrtrue
                PRINT *,'difference(1:ktrtrue) = ',difference(1:ktrtrue)
                
                
                !CALL cpu_time(time2)
                !PRINT *,'after calculating differences ', time2
                !PRINT *,'considered ', nconsider,' out of ',(imax-imin+1)*(jmax-jmin+1)

                DO isupp = 1, nsupplemental  ! number of supplemental locations

                    !PRINT *, 'isupp, nsupplemental = ', isupp, nsupplemental
                    IF (isupp .eq. 1) THEN ! the first supplemental location
                                     !  is simply that orginal grid point
                                     
                        xlocation_supp(ixa,jya,isupp) = ixa
                        ylocation_supp(ixa,jya,isupp) = jya
                        
                        ! ---- don't consider points right around grid pt of interest,
                        !      too strong a correlation (want quasi-independent samples)

                        diffmin = 0.
                        CALL mask_around_thissample(ixa, jya, minseparation, &
                            numtrue, difference, i_indices, j_indices)

                    ELSE

				        ! ---- now find the analysis grid point with the next 
				        !      closest similarity. Set this as supplemental location 
				        !      point and then mask around this to eliminate nearby
				        !      points from any future consideration

                        mindex = MINLOC(difference)
                        IF (difference(mindex(1)) .lt. 1000.) THEN
                            minx = i_indices(mindex(1))
                            miny = j_indices(mindex(1))
               
                            ! ---- finally, define the supplemental 
				            !      location on the forecast grid to be the 
                            !      grid point that had this closest fit as defined above.
                  
                            xlocation_supp(ixa,jya,isupp) = minx
                            ylocation_supp(ixa,jya,isupp) = miny

                            ! ---- make sure no other grid points very close to  
                            !      this new supplemental location are considered 
				            !      as another potential supplemental location.

                            CALL mask_around_thissample(minx, miny, minseparation, &
                                numtrue, difference, i_indices, j_indices)
                        ELSE
                            xlocation_supp(ixa,jya,isupp:nsupplemental) = -99
                            ylocation_supp(ixa,jya,isupp:nsupplemental) = -99
                            GOTO 3000
                        ENDIF
                                
                    ENDIF ! isupp > 1
                    
                    !PRINT *,'   x,y supplemental = ',xlocation_supp(ixa,jya,isupp), &
                    !    ylocation_supp(ixa,jya,isupp)
                    !PRINT *,'    difference(1:ktrtrue) = ',difference(1:ktrtrue)
                    
                    CALL cpu_time(time2)
                    
                    !PRINT *,'isupp, x, y = ',isupp, xlocation_supp(ixa,jya,isupp), ylocation_supp(ixa,jya,isupp)

                END DO ! isupp
                3000 CONTINUE
                !CALL cpu_time(time2)
                !PRINT *,'after supplemental locations ', time2
                DEALLOCATE(difference, i_indices, j_indices)
            ENDIF   ! conusmask
            
        END DO  ! jya
    END DO ! ixa



    stop
    ! --- save supplemental locations and penalty informationfor this month

    outfile = TRIM(data_directory) // 'supplemental_locations_CONUS_ndfd2p5_'//&
        cleadb //'_to_'//cleade//'_'//cmonths(imonth) // '.nc'
    CALL write_supp_locns_to_netcdf(outfile, nxa, nya, nsupplemental, &
        xlocation_supp, ylocation_supp, conusmask, lonsa, latsa)




END DO  ! imonth
PRINT *,'Done'
CONTAINS                        

! =============================================================================

SUBROUTINE mask_around_thissample(minx, miny, minseparation, &
    numtrue, difference, i_indices, j_indices)
! reset difference to very high value if this sample is close to the input (minx, miny)
! location
INTEGER, INTENT(IN) :: minx, miny, minseparation, numtrue
REAL, INTENT(INOUT), DIMENSION(numtrue) :: difference
INTEGER, INTENT(INOUT), DIMENSION(numtrue) :: i_indices, j_indices

rminseparation2 = minseparation**2
DO itrue = 1, numtrue
    dist2 = (i_indices(itrue)-minx)**2 + (j_indices(itrue)-miny)**2
    IF (dist2 .lt. rminseparation2) difference(itrue) = 99999999.
END DO

RETURN
END SUBROUTINE mask_around_thissample

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