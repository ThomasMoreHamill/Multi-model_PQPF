PROGRAM supplemental_locations_ndfd2p5

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
INTEGER, PARAMETER :: nxa = 2345 ! NDFD 2.5 km grid x dimensions for expanded grid.
INTEGER, PARAMETER :: nya = 1597 ! NDFD 2.5 km grid y dimensions for expanded grid
PARAMETER (nmonths = 12) ! # months of the year
PARAMETER (nsupplemental = 50) ! max number of CCPA grid supplemental locations  
    ! to be stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 25) ! user-defined minimum separation between 
    ! analysis grid locations and supp location (in km)
PARAMETER (maxseparation = 1000) ! user-defined maximum separation between  
    ! analysis gridlocation and supplemental location (in grid pts)

!PARAMETER (alpha_coeff = 0.15) ! coefficient to apply to Gamma dist alpha parameter
!PARAMETER (beta_coeff = 0.15) ! coefficient to apply to Gamma dist beta parameter
!PARAMETER (fz_coeff = 0.15) ! coefficient to apply to fraction zero parameter
!PARAMETER (terht_coeff = 0.25) ! coefficient to apply to terrain height differences
!PARAMETER (gradx_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
!PARAMETER (grady_coeff = 0.15)  ! coefficient to apply to terrain x gradient differences
!PARAMETER (distpenalty = .0001) ! coefficient to apply to distance 
PARAMETER (ialpha_coeff = INT(0.15*10000.)) ! coefficient to apply to Gamma dist alpha parameter
PARAMETER (ibeta_coeff = INT(0.15*10000.)) ! coefficient to apply to Gamma dist beta parameter
PARAMETER (ifz_coeff = INT(0.15*10000.)) ! coefficient to apply to fraction zero parameter
PARAMETER (iterht_coeff = INT(0.25*10000.)) ! coefficient to apply to terrain height differences
PARAMETER (igradx_coeff = INT(0.15*10000.))  ! coefficient to apply to terrain x gradient differences
PARAMETER (igrady_coeff = INT(0.15*10000.))  ! coefficient to apply to terrain x gradient differences
PARAMETER (idistpenalty = INT(.0001*10000.)) ! coefficient to apply to distance  

PARAMETER (earth_radius_meters = 6370000.)

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(nxa,nya) :: maskout  ! array used to block out nearby points !
                           ! from further consideration
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: xlocation_supp ! for each forecast point,  
                           ! a list of which other forecast points
INTEGER, DIMENSION(nxa,nya,nsupplemental) :: ylocation_supp ! have the closest climatologies 
                           ! and forecast-obs relationship
INTEGER minx
INTEGER miny

REAL, DIMENSION(nxa,nya) :: alphahat ! estimated Gamma dist shape parameter read in from file
REAL, DIMENSION(nxa,nya) :: betahat ! estimated Gamma dist scale parameter read in from file

INTEGER, DIMENSION(nxa,nya) :: ialphahat ! estimated Gamma dist shape parameter read in from file
INTEGER, DIMENSION(nxa,nya) :: ibetahat ! estimated Gamma dist scale parameter read in from file

REAL, DIMENSION(nxa,nya) :: alphahat_mean  ! local spatial mean of alpha parameter
REAL, DIMENSION(nxa,nya) :: alphahat_stddev ! local spatial std dev of alpha parameter
REAL, DIMENSION(nxa,nya) :: betahat_mean ! local spatial mean of beta parameter
REAL, DIMENSION(nxa,nya) :: betahat_stddev ! local spatial std dev of beta parameter

!INTEGER*4, DIMENSION(nxa,nya) :: ialphahat_mean  ! local spatial mean of alpha parameter
INTEGER*4, DIMENSION(nxa,nya) :: ialphahat_stddev ! local spatial std dev of alpha parameter
!INTEGER*4, DIMENSION(nxa,nya) :: ibetahat_mean ! local spatial mean of beta parameter
INTEGER*4, DIMENSION(nxa,nya) :: ibetahat_stddev ! local spatial std dev of beta parameter

!REAL, DIMENSION(nxa,nya) :: difference ! total penalty difference between (ixa,jya) 
!  and potential supplemental (ixn,jyn)
!REAL, DIMENSION(nxa,nya) :: difference_alpha ! Gamma dist alpha parameter penalty
!REAL, DIMENSION(nxa,nya) :: difference_beta ! Gamma dist beta parameter penalty
!REAL, DIMENSION(nxa,nya) :: difference_fz ! fraction zero penalty
!REAL, DIMENSION(nxa,nya) :: difference_terht  ! terrain height penalty
!REAL, DIMENSION(nxa,nya) :: difference_gradx  ! E-W terrain height gradient penalty
!REAL, DIMENSION(nxa,nya) :: difference_grady  ! N-S terrain height gradient penalty
!REAL, DIMENSION(nxa,nya) :: difference_dist ! physical distance penalty

INTEGER*4, DIMENSION(nxa,nya) :: idifference ! total penalty difference between (ixa,jya) 
!  and potential supplemental (ixn,jyn)
INTEGER*4, DIMENSION(nxa,nya) :: idifference_alpha ! Gamma dist alpha parameter penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_beta ! Gamma dist beta parameter penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_fz ! fraction zero penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_terht  ! terrain height penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_gradx  ! E-W terrain height gradient penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_grady  ! N-S terrain height gradient penalty
INTEGER*4, DIMENSION(nxa,nya) :: idifference_dist ! physical distance penalty

INTEGER, DIMENSION(nxa,nya) :: ifraction_zero ! estimated fraction zero parameter from file
REAL, DIMENSION(nxa,nya) :: fraction_zero ! estimated fraction zero parameter from file

REAL, DIMENSION(nxa,nya) :: fz_mean ! local spatial mean of fraction zero parameter
REAL, DIMENSION(nxa,nya) :: fz_stddev ! local spatial std dev of fraction zero parameter
REAL, DIMENSION(nxa,nya) :: gradx_mean ! local spatial mean of E-W terrain gradient
REAL, DIMENSION(nxa,nya) :: gradx_stddev ! local spatial std dev of E-W terrain gradient
REAL, DIMENSION(nxa,nya) :: grady_mean ! local spatial mean of N-S terrain gradient
REAL, DIMENSION(nxa,nya) :: grady_stddev ! local spatial std dev of N-S terrain gradient
!INTEGER, DIMENSION(nxa,nya) :: ifz_mean ! local spatial mean of fraction zero parameter
INTEGER, DIMENSION(nxa,nya) :: ifz_stddev ! local spatial std dev of fraction zero parameter
!INTEGER, DIMENSION(nxa,nya) :: igradx_mean ! local spatial mean of E-W terrain gradient
INTEGER, DIMENSION(nxa,nya) :: igradx_stddev ! local spatial std dev of E-W terrain gradient
!INTEGER, DIMENSION(nxa,nya) :: igrady_mean ! local spatial mean of N-S terrain gradient
INTEGER, DIMENSION(nxa,nya) :: igrady_stddev ! local spatial std dev of N-S terrain gradient
REAL, DIMENSION(nxa,nya) :: lonsa ! analysis grid longitudes (negative is W lon)
REAL, DIMENSION(nxa,nya) :: latsa ! analysis grid latitudes in degrees
REAL, DIMENSION(nxa,nya,nsupplemental) :: penalty  ! penalty function 
    ! associated with the supp location
REAL, DIMENSION(nxa,nya) :: terrain ! terrain elevation
REAL, DIMENSION(nxa,nya) :: terrain_gradx ! terrain E-W gradient
REAL, DIMENSION(nxa,nya) :: terrain_grady ! terrain N-S gradient
INTEGER, DIMENSION(nxa,nya) :: iterrain ! terrain elevation
INTEGER, DIMENSION(nxa,nya) :: iterrain_gradx ! terrain E-W gradient
INTEGER, DIMENSION(nxa,nya) :: iterrain_grady ! terrain N-S gradient
REAL, DIMENSION(nxa,nya) :: ter_mean
REAL, DIMENSION(nxa,nya) :: ter_stddev
!INTEGER, DIMENSION(nxa,nya) :: iter_mean
INTEGER, DIMENSION(nxa,nya) :: iter_stddev

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

! ---- initialize penalties

penalty = 0.


! ---- read in terrain facet information for smoothing 
!      at short, intermediate, and larger scales

infile = TRIM(data_directory) // 'blend.precip_const.2p5.nc'
PRINT *,'calling read_terrain_heights_2p5'
PRINT *,TRIM(infile)
CALL read_terrain_heights_2p5(nxa, nya, infile, terrain, &
    lonsa, latsa, conusmask)

! ---- calculate terrain gradients

PRINT *, 'calling calculate_terrain_gradients'
CALL calculate_terrain_gradients(nxa, nya, earth_radius_meters, &
    terrain, latsa, terrain_gradx, terrain_grady)
    
PRINT *,'terrain(nxa/3, nya/3) = ', terrain(nxa/3, nya/3) 
PRINT *,'terrain_gradx(nxa/3, nya/2) = ', terrain_gradx(nxa/3, nya/2)
PRINT *,'terrain_grady(nxa/3, nya/2) = ', terrain_grady(nxa/3, nya/2)
    
iterrain = INT(terrain*10000.)    
iterrain_gradx = INT(terrain_gradx*10000.)
iterrain_grady = INT(terrain_grady*10000.)
     
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
    ifraction_zero = INT(fraction_zero*10000.)
    ialphahat = INT(alphahat*10000.)
    ibetahat = INT(betahat*10000.)
             
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
            !isum_fz2 = 0
            !isum_fz  = 0
            !isum_alpha2 = 0
            !isum_alpha  = 0
            !isum_beta2 = 0
            !isum_beta = 0
            !isum_gradx = 0
            !isum_gradx2 = 0
            !isum_grady = 0
            !isum_grady2 = 0
            !isum_terht = 0
            !isum_terht2 = 0
            nsamps  = 0
            IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point has valid data
                fz_here = fraction_zero(ixa,jya)
                alpha_here = alphahat(ixa,jya)
                beta_here = betahat(ixa,jya)
                gradx_here = terrain_gradx(ixa,jya)
                grady_here = terrain_grady(ixa,jya)
                ter_here = terrain(ixa,jya)
                !ifz_here = ifraction_zero(ixa,jya)
                !ialpha_here = ialphahat(ixa,jya)
                !ibeta_here = ibetahat(ixa,jya)
                !igradx_here = iterrain_gradx(ixa,jya)
                !igrady_here = iterrain_grady(ixa,jya)
                !iter_here = iterrain(ixa,jya)
                DO ix2 = ixmin, ixmax
                    DO jy2 = jymin, jymax
                        IF (conusmask(ix2,jy2).gt. 0) THEN
                            fz_there = fraction_zero(ix2,jy2)
                            alpha_there = alphahat(ix2,jy2)
                            beta_there = betahat(ix2,jy2)
                            gradx_there = terrain_gradx(ix2,jy2)
                            grady_there = terrain_grady(ix2,jy2)
                            ter_there = terrain(ix2,jy2)
                            !ifz_there = ifraction_zero(ix2,jy2)
                            !ialpha_there = ialphahat(ix2,jy2)
                            !ibeta_there = ibetahat(ix2,jy2)
                            !igradx_there = iterrain_gradx(ix2,jy2)
                            !igrady_there = iterrain_grady(ix2,jy2)
                            !iter_there = iterrain(ix2,jy2)
					 
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
                            !isum_fz2 = isum_fz2 + ifz_there**2
                            !isum_fz  = isum_fz + ifz_there
                            !isum_alpha2 = isum_alpha2 + ialpha_there**2
                            !isum_alpha  = isum_alpha + ialpha_there
                            !isum_beta2 = isum_beta2 + ibeta_there**2
                            !isum_beta = isum_beta + ibeta_there
                            !isum_gradx = isum_gradx + igradx_there
                            !isum_gradx2 = isum_gradx2 + igradx_there**2
                            !isum_grady = isum_grady + igrady_there
                            !isum_grady2 = isum_grady2 + igrady_there**2
                            !isum_terht = isum_terht + iter_there
                            !isum_terht2 = isum_terht2 + iter_there**2
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
                    
                !ialphahat_mean(ixa,jya) = isum_alpha / nsamps
                !ialphahat_stddev(ixa,jya) = INT(SQRT((isum_alpha2-REAL(nsamps)*ialphahat_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
                !ibetahat_mean(ixa,jya) = isum_beta / nsamps
                !ibetahat_stddev(ixa,jya) = INT(SQRT((isum_beta2-REAL(nsamps)*ibetahat_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
                !ifz_mean(ixa,jya) = isum_fz / nsamps
                !ifz_stddev(ixa,jya) = INT(SQRT((isum_fz2-REAL(nsamps)*ifz_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
                !igradx_mean(ixa,jya) = isum_gradx / nsamps
                !igradx_stddev(ixa,jya) = INT(SQRT((isum_gradx2-REAL(nsamps)*igradx_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
                !igrady_mean(ixa,jya) = isum_grady / nsamps
                !igrady_stddev(ixa,jya) = INT(SQRT((isum_grady2-REAL(nsamps)*igrady_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
                !iter_mean(ixa,jya) = isum_terht  / nsamps
                !iter_stddev(ixa,jya) = INT(SQRT((isum_terht2-REAL(nsamps)*iter_mean(ixa,jya)**2) / &
                !    REAL(nsamps-1))) / 1000
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
                !ialphahat_mean(ixa,jya) = -99
                !ialphahat_stddev(ixa,jya) = -99
                !ibetahat_mean(ixa,jya) = -99
                !ibetahat_stddev(ixa,jya) = -99
                !ifz_mean(ixa,jya) = -99
                !ifz_stddev(ixa,jya) = -99
                !igradx_mean(ixa,jya) = -99
                !igradx_stddev(ixa,jya) = -99
                !igrady_mean(ixa,jya) = -99
                !igrady_stddev(ixa,jya) = -99
                !iter_mean(ixa,jya) = -99
                !iter_stddev(ixa,jya) = -99
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
    
    ialphahat_stddev = INT(10000.*alphahat_stddev)
    ibetahat_stddev = INT(10000.*betahat_stddev)
    ifz_stddev = INT(10000.*fz_stddev)
    igradx_stddev = INT(10000.*gradx_stddev)
    igrady_stddev = INT(10000.*grady_stddev)
    iter_stddev = INT(10000.*ter_stddev)

    stdmax = MAXVAL(ter_stddev)
    stdmax_sqrt = INT(REAL(stdmax)**0.5)
    gradxmax = MAXVAL(gradx_stddev)
    gradxmax_sqrt = INT(SQRT(REAL(gradxmax)))
    gradymax = MAXVAL(grady_stddev)
    gradymax_sqrt = INT(SQRT(REAL(gradymax)))
    
    ! =======================================================================================
    ! ---- find locations of supplemental locations that are the best fit for each grid point
    ! =======================================================================================
    
    CALL cpu_time(time1)
    PRINT *, 'finding supplemental locations'
    nproc = 0
   
    !DO ixa = 1, nxa
    DO ixa = nxa/3, nxa/3
        
        CALL cpu_time(time2)
        IF (MOD(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa, &
	  	    ' time elapsed, delta: ',time2, time2-time1,'nproc = ', &
            nproc,' total pts = ',SUM(REAL(conusmask))
            
        !DO jya = 1, nya
        DO jya = nya/2, nya/2

            maskout(:,:) = 0
            !difference_alpha(:,:) = 99999999. 
            !difference_beta(:,:) = 99999999.
            !difference_fz(:,:) = 99999999.
            !difference_terht(:,:) = 99999999.
            !difference_gradx(:,:) = 99999999.
            !difference_grady(:,:) = 99999999.
            !difference_dist(:,:) = 99999999.
            idifference_alpha(:,:) = 999999999
            idifference_beta(:,:) = 999999999
            idifference_fz(:,:) = 999999999
            idifference_terht(:,:) = 999999999
            idifference_gradx(:,:) = 999999999
            idifference_grady(:,:) = 999999999
            idifference_dist(:,:) = 999999999

            IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point inside area with CCPA data

                nproc = nproc+1
                imin = max(1,ixa-maxseparation)
                imax = min(nxa,ixa+maxseparation)
                jmin = max(1,jya-maxseparation)
                jmax = min(nya,jya+maxseparation)
                
                !fz_here = fraction_zero(ixa,jya)
                !alpha_here = alphahat(ixa,jya)
                !beta_here = betahat(ixa,jya)
                !gradx_here = terrain_gradx(ixa,jya)
                !grady_here = terrain_grady(ixa,jya)
                !ter_here = terrain(ixa,jya)
                ifz_here = ifraction_zero(ixa,jya)
                ialpha_here = ialphahat(ixa,jya)
                ibeta_here = ibetahat(ixa,jya)
                igradx_here = iterrain_gradx(ixa,jya)
                igrady_here = iterrain_grady(ixa,jya)
                iter_here = iterrain(ixa,jya)
                PRINT *,'imin, imax, jmin, jmax ', imin, imax, jmin, jmax
                PRINT *,'ifz_here, ialpha_here, ibeta_here, igradx_here, igrady_here, iter_here = '
                PRINT *,ifz_here, ialpha_here, ibeta_here, igradx_here, igrady_here, iter_here 

                !terht_coeff2 = terht_coeff * SQRT(ter_stddev(ixa,jya)) / istdmax_sqrt
                !gradx_coeff2 = gradx_coeff * SQRT(gradx_stddev(ixa,jya)) / igradxmax_sqrt
                !grady_coeff2 = grady_coeff * SQRT(grady_stddev(ixa,jya)) / igradymax_sqrt
                iterht_coeff2 = INT(iterht_coeff * SQRT(REAL(ter_stddev(ixa,jya))) / stdmax_sqrt)
                igradx_coeff2 = INT(igradx_coeff * SQRT(REAL(gradx_stddev(ixa,jya))) / gradxmax_sqrt)
                igrady_coeff2 = INT(igrady_coeff * SQRT(REAL(grady_stddev(ixa,jya))) / gradymax_sqrt)

                !DO ix2 = imin, imax
                !    DO jy2 = jmin, jmax
                DO ix2 = imin, imax
                    DO jy2 = nya/2, nya/2

                        ! ---- only consider this grid point if inside area with data, 
                        !      and is less than max separation, computed in 1/8-degree 
                        !      grid points.  Since there is less trustworthy data outside 
                        !      CONUS (lots of problems over ocean, over Canada at this  
                        !      point, supplement data only with points from inside CONUS

                        CALL haversine(latsa(ixa,jya), lonsa(ixa,jya), &
                            latsa(ix2,jy2), lonsa(ix2,jy2), dist)
                        idist = INT(dist*100.)

                        
                        IF (conusmask(ix2,jy2) .gt. 0 .and. maskout(ix2,jy2) &
                        .eq. 0 .and. dist .le. maxseparation) THEN
                            
                            PRINT *,'ix2, jy2, dist = ',ix2,jy2,dist
                            !fz_there = fraction_zero(ix2,jy2)
                            !alpha_there = alphahat(ix2,jy2)
                            !beta_there = betahat(ix2,jy2)
                            !gradx_there = terrain_gradx(ix2,jy2)
                            !grady_there = terrain_grady(ix2,jy2)
                            !ter_there = terrain(ix2,jy2)
                            ifz_there = ifraction_zero(ix2,jy2)
                            ialpha_there = ialphahat(ix2,jy2)
                            ibeta_there = ibetahat(ix2,jy2)
                            igradx_there = iterrain_gradx(ix2,jy2)
                            igrady_there = iterrain_grady(ix2,jy2)
                            iter_there = iterrain(ix2,jy2)
                            
                            PRINT *, 'ix2, ifz_here, ialpha_here, ibeta_here, igradx_here, igrady_here, iter_here'
                            PRINT 345, ix2, ifz_here, ialpha_here, ibeta_here, igradx_here, igrady_here, iter_here
                            PRINT *, 'ix2, ifz_there, ialpha_there, ibeta_there, igradx_there, igrady_there, iter_there'
                            PRINT 345, ix2, ifz_there, ialpha_there, ibeta_there, igradx_there, igrady_there, iter_there
                            345 format(i4,2x,6(i8,2x))
                            PRINT *, 'ix2, ifz_sd, ialpha_sd, ibeta_sd, igradx_sd, igrady_sd, iter_sd'
                            PRINT 345, ix2, ifz_stddev(ixa,jya) , ialphahat_stddev(ixa,jya), ibetahat_stddev(ixa,jya), &
                                igradx_stddev(ixa,jya), igrady_stddev(ixa,jya), iter_stddev(ixa,jya)
                            PRINT 346, ix2, fz_stddev(ixa,jya) , alphahat_stddev(ixa,jya), betahat_stddev(ixa,jya), &
                                gradx_stddev(ixa,jya), grady_stddev(ixa,jya), ter_stddev(ixa,jya)
                            346 format(i4,2x,6(f9.4,1x))

                            ! --- quantify the difference between fraction zero, alpha, beta 

                            !fz_diff = ABS(fz_here - fz_there) / fz_stddev(ixa,jya) 
                            !alpha_diff = ABS(alpha_here - alpha_there) / alphahat_stddev(ixa,jya)
                            !beta_diff = ABS(beta_here - beta_there) / betahat_stddev(ixa,jya)
                            !gradx_diff = ABS(gradx_here - gradx_there) / gradx_stddev(ixa,jya)
                            !grady_diff = ABS(grady_here - grady_there) / grady_stddev(ixa,jya)
                            !terr_diff = ABS(ter_here - ter_there) / ter_stddev(ixa,jya)
                            ifz_diff = ABS(ifz_here - ifz_there) / ifz_stddev(ixa,jya) 
                            ialpha_diff = ABS(ialpha_here - ialpha_there) / ialphahat_stddev(ixa,jya)
                            ibeta_diff = ABS(ibeta_here - ibeta_there) / ibetahat_stddev(ixa,jya)
                            igradx_diff = ABS(igradx_here - igradx_there) / igradx_stddev(ixa,jya)
                            igrady_diff = ABS(igrady_here - igrady_there) / igrady_stddev(ixa,jya)
                            iterr_diff = ABS(iter_here - iter_there) / iter_stddev(ixa,jya)

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
                            idifference_alpha(ix2,jy2) = ialpha_coeff*ialpha_diff 
                            idifference_beta(ix2,jy2) = ibeta_coeff*ibeta_diff
                            idifference_fz(ix2,jy2) = ifz_coeff*ifz_diff
                            idifference_terht(ix2,jy2) = iterht_coeff2*iterr_diff
                            idifference_gradx(ix2,jy2) = igradx_coeff2*igradx_diff
                            idifference_grady(ix2,jy2) = igrady_coeff2*igrady_diff
                            idifference_dist(ix2,jy2) = idist*idistpenalty
                            idifference(ix2,jy2) = idifference_alpha(ix2,jy2) + &
                                idifference_beta(ix2,jy2) + idifference_fz(ix2,jy2) + &
                                idifference_terht(ix2,jy2) + idifference_gradx(ix2,jy2) + &
                                idifference_grady(ix2,jy2) + idifference_dist(ix2,jy2)                     
                     
                        ELSE
                      
                            !difference_alpha(ix2,jy2) = 99999999.
                            !difference_beta(ix2,jy2) = 99999999.
                            !difference_fz(ix2,jy2) = 99999999.
                            !difference_terht(ix2,jy2) = 99999999.
                            !difference_gradx(ix2,jy2) = 99999999.
                            !difference_grady(ix2,jy2) = 99999999.
                            !difference_dist(ix2,jy2) = 99999999.
                            !difference(ix2,jy2) = 99999999.
                            idifference_alpha(ix2,jy2) = 999999999
                            idifference_beta(ix2,jy2) = 999999999
                            idifference_fz(ix2,jy2) = 999999999
                            idifference_terht(ix2,jy2) = 999999999
                            idifference_gradx(ix2,jy2) = 999999999
                            idifference_grady(ix2,jy2) = 999999999
                            idifference_dist(ix2,jy2) = 999999999
                            idifference(ix2,jy2) = 999999999

                        ENDIF
                    END DO ! jy2
                END DO ! ix2
                
                PRINT *,' ix2     alpha   beta   fz    terht    gradx   grady   dist    total'
                DO ix2 = imin, imax
                    IF(conusmask(ix2,jya) .eq. 1) THEN
                        PRINT 234,ix2,idifference_alpha(ix2,jya), idifference_beta(ix2,jya),&
                        idifference_fz(ix2,jya), idifference_terht(ix2,jya), &
                        idifference_gradx(ix2,jya), idifference_grady(ix2,jya), &
                        idifference_dist(ix2,jya), idifference(ix2,jya)
                        234 format(i4,2x,8(i9,1x))
                    ENDIF
                END DO
                stop
                
                

                DO isupp = 1, nsupplemental  ! number of supplemental locations

                    IF (isupp .eq. 1) THEN ! the first supplemental location
                                     !  is simply that orginal grid point
                                     
                        xlocation_supp(ixa,jya,isupp) = ixa
                        ylocation_supp(ixa,jya,isupp) = jya

                        ! ---- don't consider points right around grid pt of interest,
                        !      too strong a correlation (want quasi-independent samples)

                        idiffmin = 0
                        CALL mask_around_thisgridpt(ixa, jya, nxa, nya, &
                            minseparation, maskout)

                    ELSE

				        ! ---- now find the analysis grid point with the next 
				        !      closest similarity. Set this as supplemental location 
				        !      point and then mask around this to eliminate nearby
				        !      points from any future consideration

                        minx = -99
                        miny = -99
                        !diffmin = 999999.
                        idiffmin = 999999999

                        DO ix2 = imin, imax   !1, nxa
                            DO jy2 = jmin, jmax  !1, nya

                                !IF (difference(ix2,jy2) .lt. diffmin .and. &
                                !maskout(ix2,jy2) .eq. 0 .and. &
                                !conusmask(ix2,jy2) .gt. 0) THEN
                                IF (idifference(ix2,jy2) .lt. idiffmin .and. &
                                maskout(ix2,jy2) .eq. 0 .and. &
                                conusmask(ix2,jy2) .gt. 0) THEN
                                    !diffmin = difference(ix2,jy2)
                                    idiffmin = idifference(ix2,jy2)
                                    minx = ix2
                                    miny = jy2
                                END IF
                            END DO ! jy2
                        END DO  ! ix2
               
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
!$OMP END PARALLEL DO
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

END PROGRAM supplemental_locations_ndfd2p5
