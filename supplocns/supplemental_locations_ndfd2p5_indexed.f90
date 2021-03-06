PROGRAM supplemental_locations_ndfd2p5_indexed

! purpose:  For every (i,j) in on the 2.5-km NDFD Climatology-Calibrated
!    Precipitation Analysis (CCPA; see Hou et al. 2014, available at
!    http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-11-0140.1)  grid 
!    with valid data, determine a set of supplemental data locations.
!    Base these supplemental locations on the similarity of analyzed CDFs
!    (via gamma-distribution cdf differences, and fraction zero parameters) 
!    as well as the local terrain heights and gradients.  Generate the 
!    supplemental locations separately for every month of the year. We want 
!    grid points that will have somewhat independent data, so make 
!    sure that the set of supplemental locations are not too close together.
! 
USE netcdf

PARAMETER (nyears = 15) ! 2002-2016 here for the CCPA data
PARAMETER (nxa = 2345) ! NDFD 2.5 km grid x dimensions for expanded grid.
PARAMETER (nya = 1597) ! NDFD 2.5 km grid y dimensions for expanded grid
PARAMETER (npthresh = 8) ! number of precip thresholds CDF evaluated at
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


PARAMETER (cdf_coeff = 0.30) ! coefficient to apply to Gamma CDF diffferences
PARAMETER (fz_coeff = 0.15) ! coefficient to apply to fraction zero parameter
PARAMETER (terht_coeff = 0.25) ! coefficient to apply to terrain height differences
PARAMETER (gradx_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (grady_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
PARAMETER (distpenalty = .001) ! coefficient to apply to distance  
PARAMETER (earth_radius_meters = 6370000.)

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(nxa,nya) :: practical_mask  ! 1/0 mask for whether inside NDFD mask
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

REAL, DIMENSION(nxa,nya,npthresh) :: cdf ! 
REAL, DIMENSION(nxa,nya) :: cdf_mae !  CDF MAE from center point
REAL, DIMENSION(nxa,nya) :: cdf_diff !  CDF RMSE from center point
REAL :: cdfdiff, cdf_max

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
REAL, DIMENSION(npthresh) :: pthreshes ! precip thresholds CDF evaluated at
REAL, DIMENSION(nxa,nya) :: terrain ! terrain elevation
REAL, DIMENSION(nxa,nya) :: terrain_gradx ! terrain E-W gradient
REAL, DIMENSION(nxa,nya) :: terrain_grady ! terrain N-S gradient
REAL*8, DIMENSION(nxa,nya) :: ter_mean
REAL*8, DIMENSION(nxa,nya) :: ter_stddev

REAL*8 :: sum_fz2 
REAL*8 :: sum_fz
REAL*8 :: sum_cdf2
REAL*8 :: sum_cdf
REAL*8 :: sum_gradx 
REAL*8 :: sum_gradx2
REAL*8 :: sum_grady 
REAL*8 :: sum_grady2 
REAL*8 :: sum_terht 
REAL*8 :: sum_terht2

REAL*8 :: fz_here 
REAL*8, DIMENSION(npthresh) :: cdf_here
REAL*8 :: gradx_here
REAL*8 :: grady_here
REAL*8 :: ter_here

REAL*8 :: fz_there 
REAL*8, DIMENSION(npthresh) :: cdf_there
REAL*8 :: gradx_there
REAL*8 :: grady_there
REAL*8 :: ter_there

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
DO imonth = 1,12

    ! ---- read in the parameters of the Gamma distribution for 
    !      climatology, and fraction zero.

    infile = TRIM(data_directory) // &
        'climatology_gamma_parameters_ndfd2p5_' // &
        cleadb // '_to_' // cleade // '_' // &
        cmonths(imonth) //'.nc'
    PRINT *, TRIM(infile)
    CALL read_climatology_parameters_ndfd2p5(nxa, nya, npthresh, infile, &
        conusmask, practical_mask, lonsa, latsa, fraction_zero, cdf, &
        pthreshes)  
    
    ! ---- loop thru grid points and precompute statistics for mean, std dev
    !      of fraction zero, and terrain and RMSE of CDF diffs.  These will be used to
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
            cdf_here(:) = cdf(ixa,jya,:)
            sum_fz2 = 0.0
            sum_fz  = 0.0
            sum_cdf = 0.0
            sum_gradx = 0.0
            sum_gradx2 = 0.0
            sum_grady = 0.0
            sum_grady2 = 0.0
            sum_terht = 0.0
            sum_terht2 = 0.0
            nsamps  = 0
            IF (practical_mask(ixa,jya) .gt. 0) THEN ! grid point has valid data
                DO ix2 = ixmin, ixmax
                    DO jy2 = jymin, jymax
                        IF (conusmask(ix2,jy2).gt. 0) THEN

                            fz_there = fraction_zero(ix2,jy2)
                            cdf_there(:) = cdf(ix2,jy2,:)
                            gradx_there = terrain_gradx(ix2,jy2)
                            grady_there = terrain_grady(ix2,jy2)
                            ter_there = terrain(ix2,jy2)
					 
                            ! --- quantify the difference between forecast distributions at 
                            !     the two grid points
					 
                            sum_fz2 = sum_fz2 + fz_there**2
                            sum_fz  = sum_fz + fz_there
                            !sum_cdf2 = sum_cdf2 + SUM(ABS(cdf_here(:) - cdf_there(:))**2)
                            sum_cdf = sum_cdf + SUM(ABS(cdf_here(:) - cdf_there(:)))
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

                cdf_mae(ixa,jya) = sum_cdf/ REAL(nsamps-1)
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
                cdf_mae(ixa,jya) = -99.99
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

    DO ixa = 1, nxa
    !DO ixa = 2143, 2143
    !DO ixa = nxa/2, nxa/2
        CALL cpu_time(time2)
        IF (MOD(ixa-1,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa, &
	  	    ' time elapsed, delta: ',time2, time2-time1,'nproc = ', &
            nproc,' total pts = ',SUM(REAL(practical_mask))
        DO jya = 1, nya
        !DO jya = 903, 903
        !DO jya = nya/2, nya/2

            IF (practical_mask(ixa,jya) .gt. 0) THEN ! grid point inside area with CCPA data
                
                !PRINT *,'conus, practical mask = ', conusmask(ixa,jya), practical_mask(ixa,jya)
                !PRINT *,'processing ixa, jya = ', ixa, jya
                !PRINT *,'lon, lat = ', lonsa(ixa,jya), latsa(ixa,jya)
                
                ! ---- for points truly in the CONUS, with better data than offshore points,
                !      make sure that the we are considering for supplemental locations only
                !      other conus points.
                
                IF (conusmask(ixa,jya) .gt. 0) THEN
                    maskout(:,:) = 1 - conusmask(:,:) 
                ELSE 
                    maskout(:,:) = 1 - practical_mask(:,:)
                ENDIF            
                
                ! ---- rather than examining all points for their closeness to the 
                !      current grid point's cdf, fraction zero, and terrain
                !      values, let's thin down the number of grid points to examine
                !      to a smaller set that have similar values
            
                cdf_max = cdf_mae(ixa,jya) * 0.5
                fraction_zero_min = MAX(0.00001,fraction_zero(ixa,jya) - fz_stddev(ixa,jya))
                fraction_zero_max = MIN(fraction_zero(ixa,jya) + fz_stddev(ixa,jya), 1.0)
                !PRINT *,'cdf_max = ', cdf_max
                !PRINT *,'fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max =  ',&
                !    fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max
                
                nproc = nproc+1
                imin = max(1, ixa-maxseparation)
                imax = min(nxa, ixa+maxseparation)
                jmin = max(1, jya-maxseparation)
                jmax = min(nya, jya+maxseparation)
                !PRINT *,'imin, ixa, imax = ', imin, ixa, imax
                !PRINT *,'jmin, jya, jmax = ', jmin, jya, jmax
                
                fz_here = fraction_zero(ixa,jya)
                cdf_here = cdf(ixa,jya,:)
                gradx_here = terrain_gradx(ixa,jya)
                grady_here = terrain_grady(ixa,jya)
                ter_here = terrain(ixa,jya)

                terht_coeff2 = terht_coeff * SQRT(ter_stddev(ixa,jya)) / stdmax_sqrt
                gradx_coeff2 = gradx_coeff * SQRT(gradx_stddev(ixa,jya)) / gradxmax_sqrt
                grady_coeff2 = grady_coeff * SQRT(grady_stddev(ixa,jya)) / gradymax_sqrt

                nconsider = 0
                dist(:,:) = 999999999.
                cdf_diff(:,:) = 9999999.
                !CALL cpu_time(time2)
                !PRINT *,'before distance calculations ', time2
                DO ix2 = imin, imax
                    DO jy2 = jmin, jmax
                        dist(ix2,jy2) = MAX( ABS(ixa-ix2), ABS(jya-jy2))
                        cdf_there = cdf(ix2,jy2,:)
                        cdf_diff(ix2,jy2) = SUM(ABS( cdf_here(:)-cdf_there(:) ))
                        !CALL haversine(latsa(ixa,jya), lonsa(ixa,jya), &
                        !    latsa(ix2,jy2), lonsa(ix2,jy2), dist(ix2,jy2))
                    END DO
                END DO
                !CALL cpu_time(time2)
                !PRINT *,'after distance calculations ', time2
                        
                !PRINT *,'fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max =  ',&
                !    fraction_zero_min, fraction_zero(ixa,jya), fraction_zero_max
                
                ! ---- form a Boolean array as to whether or not to evaluate this point
                !      as a potential supplemental location.
                
                IF (conusmask(ixa,jya) .gt. 0) THEN
                    test_diffs_here = conusmask > 0 .and. dist < maxseparation .and. &
                    cdf_diff <= cdf_max .and. fraction_zero >= fraction_zero_min .and. &
                    fraction_zero <= fraction_zero_max
                ELSE 
                    test_diffs_here = practical_mask > 0 .and. dist < maxseparation .and. &
                    cdf_diff <= cdf_max .and. fraction_zero >= fraction_zero_min .and. &
                    fraction_zero <= fraction_zero_max
                ENDIF    
                 
                numtrue = count(test_diffs_here)
                !PRINT *,'numtrue = ',numtrue
                IF (numtrue .gt. 0) THEN 
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
                                cdf_there = cdf(ix2,jy2,:)
                                gradx_there = terrain_gradx(ix2,jy2)
                                grady_there = terrain_grady(ix2,jy2)
                                ter_there = terrain(ix2,jy2)

                                ! --- quantify the difference between fraction zero, cdf, terrain characteristics

                                fz_diff = ABS(fz_here - fz_there) / fz_stddev(ixa,jya) 
                                cdfdiff =  SUM(ABS(cdf_here(:)-cdf_there(:) )) / cdf_mae(ixa,jya)
                                gradx_diff = ABS(gradx_here - gradx_there) / gradx_stddev(ixa,jya)
                                grady_diff = ABS(grady_here - grady_there) / grady_stddev(ixa,jya)
                                terr_diff = ABS(ter_here - ter_there) / ter_stddev(ixa,jya)

                                ! --- calculate differences between this potential supplemental location
                                !     and the grid point of interest for each of the various factors 
                                !     we choose to weight 
                            
                                difference(ktrtrue) = cdf_coeff*cdfdiff + fz_coeff*fz_diff + &
                                    terht_coeff2*terr_diff + gradx_coeff2*gradx_diff + &
                                    grady_coeff2*grady_diff + dist(ix2,jy2)*distpenalty   
                                i_indices(ktrtrue) = ix2
                                j_indices(ktrtrue) = jy2                 
                            
                            ENDIF ! test_diffs_here
                        END DO ! jy2
                    END DO ! ix2
                    !print *,'ktrtrue = ', ktrtrue
                    !PRINT *,'difference(1:ktrtrue) = ',difference(1:ktrtrue)
                
                    CALL cpu_time(time2)
                    !PRINT *,'after calculating differences ', time2
                    !PRINT *,'considered ', nconsider,' out of ',SUM(INT(conusmask(imin:imax, jmin:jmax)))

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
                    
                        !PRINT *,'isupp, x,y supplemental = ',isupp, xlocation_supp(ixa,jya,isupp), &
                        !    ylocation_supp(ixa,jya,isupp)
                        !PRINT *,'    difference(1:ktrtrue) = ',difference(1:ktrtrue)
                    
                        CALL cpu_time(time2)
                        IF (MOD(ixa,100) .eq. 0 .and. jya .eq. nya/2) THEN
                        !IF (ixa .eq. nxa/2 .and. jya .eq. nya/2) THEN
                        !IF (ixa .eq. 1201 .and. jya .eq. 801) THEN
                            PRINT *,'ixa, jya, isupp, xsupp, ysupp = ',ixa, jya, isupp, &
                                xlocation_supp(ixa,jya,isupp), ylocation_supp(ixa,jya,isupp)
                        ENDIF

                    END DO ! isupp
                    3000 CONTINUE
                    !CALL cpu_time(time2)
                    !PRINT *,'after supplemental locations ', time2
                    DEALLOCATE(difference, i_indices, j_indices)
                ELSE ! numtrue = 0
                    xlocation_supp(ixa,jya,1) = ixa
                    ylocation_supp(ixa,jya,1) = jya
                    xlocation_supp(ixa,jya,2:nsupplemental) = -99
                    ylocation_supp(ixa,jya,2:nsupplemental) = -99
                ENDIF
            ENDIF   ! conusmask
            
        END DO  ! jya
    END DO ! ixa
    
    ! --- save supplemental locations and penalty informationfor this month

    PRINT *,'min, max lonsa = ', minval(lonsa), maxval(lonsa)
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