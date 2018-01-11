PROGRAM create_supp_locns_ndfd2p5

! purpose:  For every (i,j) in on the 2/.5-km NDFD grid with valid data, 
!    determine a set of supplemental data locations.  Here we will 
!    leverage supplemental locations that were previously generated for
!    the coarser 1/8-degree CONUS grid, forcing the use of similar locations
!    on the 2.5 km grid where the 1/8-degree data was within its CONUS
!    mask.  For 2.5-km NDFD points that are nearest to a 1/8-degree point
!    outside the CONUS mask, we fall back to an algorithm of searching more
!    broadly for supplemental locations.
!
!    Base the supplemental locations on the similarity of analyzed CDFs 
!    as well as the local terrain heights and facets.  Do this separately 
!    for every month of the year. We want grid points that will have 
!    somewhat independent data, so make sure that the set of supplemental 
!    locations are not too close together.
! 
USE netcdf

PARAMETER (nx_ccpa = 515) ! CCPA 1/8-degree grid x dimensions for expanded grid.
PARAMETER (ny_ccpa = 262) ! CCPA 1/8-degree grid y dimensions for expanded grid
PARAMETER (nx_2p5 = 2345) ! NDFD 2.5 km grid x dimensions for expanded grid.
PARAMETER (ny_2p5 = 1597) ! NDFD 2.5 km grid y dimensions for expanded grid
PARAMETER (nmonths = 12)   ! # months of the year
PARAMETER (nsupp = 50) ! max number of CCPA grid supplemental locations to be 
                           ! stored for a given forecast grid box ! prev 20
PARAMETER (nsupp_ccpa = 100) ! max number of CCPA grid supplemental locations to be 
                        ! stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 28) ! user-defined minimum separation between 
    ! analysis grid locations and supp location (in grid pts)
PARAMETER (maxseparation = 800) ! user-defined maximum separation between analysis grid 
    ! location and supplemental location (in grid pts)
PARAMETER (npct = 20)     ! number of quantiles of CDF to compute at
PARAMETER (pa_coeff = 0.05) ! coefficient to apply to analysis CDF differences
PARAMETER (facet_coeff = 0.6) ! coefficient to apply to terrain facet differences
PARAMETER (terht_coeff = 0.4)  ! coefficient to apply to terrain height differences
PARAMETER (distpenalty = .0001)  ! coefficient to apply to distance 

INTEGER(2), DIMENSION(nx_ccpa,ny_ccpa) :: conusmask_ccpa  ! 1/0 mask for whether in CONUS
INTEGER(2), DIMENSION(nx_2p5,ny_2p5) :: conusmask_2p5  ! 1/0 mask for whether in NDFD grid
INTEGER(2), DIMENSION(nx_2p5,ny_2p5) :: maskout_2p5  ! array used to block out nearby points 
    ! from further consideration

INTEGER, DIMENSION(nx_ccpa,ny_ccpa,nsupp_ccpa) :: xlocation_ccpa ! for each forecast point,  
    ! a list of which other forecast points
INTEGER, DIMENSION(nx_ccpa,ny_ccpa,nsupp_ccpa) :: ylocation_ccpa ! have the closest climatologies 
    ! and forecast-obs relationship
INTEGER, DIMENSION(nx_2p5,ny_2p5,nsupp) :: xlocation_2p5 ! for each forecast point,  
    ! a list of which other forecast points
INTEGER, DIMENSION(nx_2p5,ny_2p5,nsupp) :: ylocation_2p5 ! have the closest climatologies 
    ! and forecast-obs relationship
    
INTEGER, DIMENSION(nx_2p5,ny_2p5) :: imin_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt
INTEGER, DIMENSION(nx_2p5,ny_2p5) :: imax_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt

INTEGER(4), DIMENSION(nx_2p5,ny_2p5) :: ifacet_short_2p5  ! facet associated with fine smoothing
INTEGER(4), DIMENSION(nx_2p5,ny_2p5) :: ifacet_medium_2p5 ! ... medium ....
INTEGER(4), DIMENSION(nx_2p5,ny_2p5) :: ifacet_long_2p5 ! ... long ...
INTEGER, DIMENSION(nx_2p5,ny_2p5) :: inearest_ccpa ! nearest ccpa point to this 2.5 km point
INTEGER, DIMENSION(nx_2p5,ny_2p5) :: jnearest_ccpa
INTEGER, DIMENSION(nx_ccpa,ny_ccpa) :: inearest_2p5 ! nearest 2.5 km point to this ccpa point
INTEGER, DIMENSION(nx_ccpa,ny_ccpa) :: jnearest_2p5

REAL, ALLOCATABLE, DIMENSION(:,:)   :: difference

REAL, DIMENSION(nx_2p5,ny_2p5)   :: fraction_zero
REAL, DIMENSION(nx_ccpa,ny_ccpa) :: lons_ccpa
REAL, DIMENSION(nx_ccpa,ny_ccpa) :: lats_ccpa
REAL, DIMENSION(nx_2p5,ny_2p5)   :: lons_2p5
REAL, DIMENSION(nx_2p5,ny_2p5)   :: lats_2p5
REAL, DIMENSION(nx_2p5,ny_2p5,npct) :: panal_quantiles
REAL, DIMENSION(npct)            :: pctv
REAL, DIMENSION(npct)            :: pa_here
REAL, DIMENSION(npct)            :: pa_there
REAL, DIMENSION(nx_2p5,ny_2p5)   :: pa_stddev
REAL, DIMENSION(nx_2p5,ny_2p5)   :: pa_mean
!REAL, DIMENSION(nx_2p5,ny_2p5, nsupp) :: penalty
REAL, DIMENSION(nx_2p5,ny_2p5)   :: terrain_height
REAL, DIMENSION(nx_2p5,ny_2p5)   :: ter_mean
REAL, DIMENSION(nx_2p5,ny_2p5)   :: ter_stddev

REAL time1, time2

REAL diffmin

CHARACTER (len=120) :: outfile
CHARACTER (len=3) :: cmonths(12)
CHARACTER (len=2) :: cmonthsin(12)
CHARACTER (len=3) :: cmonth
CHARACTER (len=20) :: cfields(npct)
CHARACTER (len=2) :: cmonthin
CHARACTER (len=3) :: csupp

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun', &
             'Jul','Aug','Sep','Oct','Nov','Dec'/
DATA cmonthsin /'01','02','03','04','05','06',&
    '07','08','09','10','11','12'/
             
DATA pctv / 0.1,0.3,0.5,0.7,1., 2.,3.,4.,5.,7., &
    10.,15.,20.,25.,30., 40.,50.,75.,100.,150./
    
DATA cfields /'nonexceedance_0p1mm','nonexceedance_0p3mm',&
    'nonexceedance_0p5mm','nonexceedance_0p7mm','nonexceedance_1mm', &
    'nonexceedance_2mm','nonexceedance_3mm','nonexceedance_4mm',&
    'nonexceedance_5mm','nonexceedance_7mm','nonexceedance_10mm',&
    'nonexceedance_15mm','nonexceedance_20mm','nonexceedance_25mm',&
    'nonexceedance_30mm','nonexceedance_40mm','nonexceedance_50mm',&
    'nonexceedance_75mm','nonexceedance_100mm','nonexceedance_150mm'/         

CALL getarg(1, cmonthin)
READ (cmonthin,'(i2)') imonthin
iprecompute = 1 ! = 0 if terrain std dev and such have NOT been precomputed

IF (nsupp .ge. 100) THEN
   WRITE (csupp,'(i3)') nsupp
ELSE
   WRITE (csupp,'(i2)') nsupp
   csupp = '0' // csupp(1:2)
ENDIF
!print *,'csupp = ',csupp

xlocation_ccpa = -99  ! initialize all locations to missing value.
ylocation_ccpa = -99
xlocation_2p5 = -99
ylocation_2p5 = -99
    
! ---- for each month, load the quantile information and rank correlation previously
!      computed as well as terrain height and facet information for the 
!      2.5-km grid

PRINT *, 'calling read_terr_ht_and_facets_2p5'
CALL read_terr_ht_and_facets_2p5 (nx_2p5, ny_2p5, &
    terrain_height, ifacet_short_2p5, ifacet_medium_2p5, ifacet_long_2p5, &
    conusmask_2p5, lons_2p5, lats_2p5)

DO imonth = imonthin, imonthin
!DO imonth =  4, 12
    
    cmonth = cmonths(imonth)
    
    ! ---- read in the supplemental locations previously defined on the 
    !      coarser CCPA grid.  Also get lat, lon, and associated 
    !      CONUS mask

    CALL read_supp_locns_ccpa(cmonth, nx_ccpa, ny_ccpa, nsupp_ccpa, &
        xlocation_ccpa, ylocation_ccpa, lons_ccpa, lats_ccpa, &
        conusmask_ccpa)
    
    ! ---- Determine the nearest CCPA grid indices for every NDFD 2.5-km grid point.
    !      As a by-product, also determine the 2.5-km grid point that is nearest to
    !      each CCPA point.

    CALL determine_nearest_ccpa(nx_2p5, ny_2p5, nx_ccpa, ny_ccpa, &
        lons_2p5, lats_2p5, lons_ccpa, lats_ccpa, conusmask_ccpa, &
        inearest_ccpa, jnearest_ccpa, inearest_2p5, jnearest_2p5)
        
    ! ---- load the CDF data on the 2.5-km grid

    CALL load_precipquantiles_2p5 (cmonthsin(imonth), nx_2p5, ny_2p5, npct, &
        cfields, panal_quantiles, fraction_zero)
    
    ! ---- precompute the minimum and maximum indices for quantiles 
    !      where we will consider evaluating the CDF at
   
    DO jy = 1, ny_2p5
        DO ix = 1, nx_2p5
            pa_here(:) = panal_quantiles(ix, jy,:)
            CALL compute_imin_imax_of_CDF (npct, &
                pa_here, imin_pa(ix,jy), imax_pa(ix,jy))
        END DO
    END DO

    ! ---- loop thru grid points and precompute statistics for mean, std dev
    !      of analysis CDF and terrain

    IF (iprecompute .eq. 0) THEN
        PRINT *,'pre-processing to determine local standard devs of CDF and terrain elevation'
        DO jy = 1, ny_2p5
            IF (MOD(jy,50) .eq. 0) PRINT *,'processing  jy = ',jy,' of ',ny_2p5
            jymin = MAX(1,jy-30)
            jymax = MIN(ny_2p5,jy+30)
            DO ix = 1, nx_2p5

                ! --- only bother searching in a box +/- 10 grid points to limit computations.

                ixmin = MAX(1,ix-30)
                ixmax = MIN(nx_2p5,ix+30)
                sum_pa2 = 0.0
                sum_pa  = 0.0
                sum_t   = 0.0
                sum_t2  = 0.0
                nsamps  = 0
                ter_here = 0.
                ter_there = 0.
                IF (conusmask_2p5(ix,jy) .gt. 0) THEN ! grid point has valid data
                    pa_here(:) = panal_quantiles(ix, jy,:)
                    ter_here = terrain_height(ix,jy)
                    DO jy2 = jymin, jymax
                        DO ix2 = ixmin, ixmax
                            IF (conusmask_2p5(ix2,jy2).gt. 0) THEN
                                pa_there(:) = panal_quantiles(ix2, jy2,:)
                                ter_there = terrain_height(ix2,jy2)
					 
                                ! --- quantify the difference between forecast distributions at 
                                !     the two grid points
					 
                                padiff = SUM(ABS(pa_here(imin_pa(ix,jy):imax_pa(ix,jy)) - &
                                    pa_there(imin_pa(ix,jy):imax_pa(ix,jy)))) / &
                                    REAL(imax_pa(ix,jy)-imin_pa(ix,jy)+1)
                                sum_t2     = sum_t2     + ter_there**2
                                sum_t      = sum_t      + ter_there
                                sum_pa2    = sum_pa2    + padiff**2
                                sum_pa     = sum_pa     + padiff
                                nsamps     = nsamps + 1
                            ENDIF
                        END DO ! jy2
                    END DO ! ix2

                    ! --- now calculate std deviation by shortcut (Wilks 2006, eq. 3.20)

                    pa_mean(ix,jy) = sum_pa / REAL(nsamps)
                    pa_stddev(ix,jy) = &
                        SQRT((sum_pa2-REAL(nsamps)*pa_mean(ix,jy)**2) / &
                        REAL(nsamps-1))
                    ter_mean(ix,jy) = sum_t  / REAL(nsamps)
                    ter_stddev(ix,jy) = &
                        SQRT((sum_t2-REAL(nsamps)*ter_mean(ix,jy)**2) / &
                        REAL(nsamps-1))
                ELSE
                    pa_mean(ix,jy)       = -99.99
                    pa_stddev(ix,jy)     = -99.99
                    ter_mean(ix,jy)      = -99.99
                    ter_stddev(ix,jy)    = -99.99
                END IF
            END DO ! ix
        END DO    ! jy

        ! ---- save this data to file
        
        outfile = 'pa_terrain_mean_stddev.dat'
        OPEN (UNIT=1, FILE=outfile, STATUS='replace', FORM='unformatted')
        WRITE (1) pa_mean
        WRITE (1) pa_stddev
        WRITE (1) ter_mean
        WRITE (1) ter_stddev
        CLOSE (1) 

    ELSE  ! iprecompute != 0
        
        outfile = 'pa_terrain_mean_stddev.dat'
        PRINT *,'reading pre-computed terrain and Pa data from ', TRIM(outfile)
        OPEN (UNIT=1, FILE=outfile, STATUS='old', FORM='unformatted')
        READ (1) pa_mean
        READ (1) pa_stddev
        READ (1) ter_mean
        READ (1) ter_stddev
        CLOSE (1)

    ENDIF ! iprecompute
    
    PRINT *,' ------  pa_mean ------'
    DO j = ny_2p5/2+10, ny_2p5/2-10,-1
        PRINT 123, j,(pa_mean(i,j),i=380,400)
        123 format(i4,1x,21(f8.3,1x))
    END DO

    PRINT *,' ------  pa_stddev ------'
    DO j = ny_2p5/2+10, ny_2p5/2-10,-1
        PRINT 123, j,(pa_stddev(i,j),i=380,400)
    END DO

    PRINT *,' ------  ter_mean ------'
    DO j = ny_2p5/2+10, ny_2p5/2-10,-1
        PRINT 123, j,(ter_mean(i,j),i=380,400)
    END DO

    PRINT *,' ------  ter_stddev ------'
    DO j = ny_2p5/2+10, ny_2p5/2-10,-1
        PRINT 123, j,(ter_stddev(i,j),i=380,400)
    END DO

    PRINT *, 'These arrays span longitudes from ',lons_2p5(380,ny_2p5/2), &
         ' to ',lons_2p5(400, ny_2p5/2)
    
    stdmax = MAXVAL(ter_stddev)
    stdmax_0p5 = stdmax**0.5
    
    ! ---- find locations of supplemental locations 
    !      that are the best fit for each grid point

    CALL cpu_time(time1)
    PRINT *, '*************** finding supplemental locations ****************'
    
    nproc = 0
    DO jy = 1, ny_2p5
        CALL cpu_time(time2)
        IF (MOD(jy,20) .eq. 0) PRINT *,'jy = ',jy,' of ',ny_2p5,&
      	    ' time elapsed, delta: ',time2, time2-time1
        DO ix = 1, nx_2p5
            
            IF (conusmask_2p5(ix,jy) .gt. 0) THEN ! grid point inside area with valid data
                maskout_2p5(:,:) = 0
                IF (inearest_ccpa(ix,jy) .ne. -99) THEN ! there is a nearby ccpa point
                    
                    ixccpa = inearest_ccpa(ix,jy) ! nearest ccpa index to this 2.5-km point
                    jyccpa = jnearest_ccpa(ix,jy)
                    
                    ! ---- for this grid point we have already found supplemental locations
                    !      on the CCPA grid.  Let's leverage this calculation, forcing
                    !      each of the NDFD grid supplemental locations to be close to
                    !      the the CCPA location.  
                    
                    DO isupp = 1, nsupp 
                        
                        ! ---- only process if 1/8-degree CCPA data was valid for this
                        !      supplemental location number.
                                                    
                        IF (xlocation_ccpa(ixccpa,jyccpa,isupp) .gt. 0) THEN 

                            ! ---- determine the index of the 2.5-km grid point that is 
                            !      approximately nearest to this 1/8-degree CCPA supplemental
                            !      location
                        
                            IF (isupp .eq. 1) THEN 
                                ipmin = ix
                                jpmin = jy
                                xlocation_2p5(ix,jy,isupp) = ix
                                ylocation_2p5(ix,jy,isupp) = jy
                            ELSE
                                ixc = xlocation_ccpa(ixccpa,jyccpa,isupp)
                                jyc = ylocation_ccpa(ixccpa,jyccpa,isupp)
                                ixs_2p5_center = inearest_2p5(ixc, jyc)
                                jys_2p5_center = jnearest_2p5(ixc, jyc)

                                imin = MAX(1,ixs_2p5_center-3)
                                imax = MIN(nx_2p5,ixs_2p5_center+3)
                                jmin = MAX(1,jys_2p5_center-3)
                                jmax = MIN(ny_2p5,jys_2p5_center+3)
                                isumc = SUM(INT(conusmask_2p5(imin:imax,jmin:jmax)))
                                IF (isumc .eq. 0) THEN 
                                    ! We have a weird situation where there is no 2.5-km data near this CCPA point.
                                    ! This happens, for instance, in the N. Columbia River Basin where for some 
                                    ! reason the NDFD is trimmed.
                                    imin = MAX(1,ixs_2p5_center-maxseparation)
                                    imax = MIN(nx_2p5,ixs_2p5_center+maxseparation)
                                    jmin = MAX(1,jys_2p5_center-maxseparation)
                                    jmax = MIN(ny_2p5,jys_2p5_center+maxseparation)
                                ENDIF
                                    
                                nxdiff = imax - imin +1
                                nydiff = jmax - jmin +1
                                ALLOCATE (difference(nxdiff,nydiff))
                                difference(:,:) = 99999.
            
                                ! ---- for the set of locations (imin:imax, jmin:jmax), 
                                !      find the single 2.5-km grid point that has the smallest
                                !      penalty difference from (ix,jy)
                        
                                CALL determine_penalties(ix, jy, nx_2p5, ny_2p5, nxdiff,  &
                                    nydiff, npct, imin, imax, jmin, jmax, maxseparation, &
                                    ifacet_short_2p5, ifacet_medium_2p5, &
                                    ifacet_long_2p5, panal_quantiles, ter_stddev, pa_stddev, &
                                    distpenalty, facet_coeff, pa_coeff, terht_coeff,  &
                                    terrain_height, conusmask_2p5, imin_pa, imax_pa, &
                                    stdmax_0p5, lons_2p5, lats_2p5, difference, penalty_analdiff, &
                                    penalty_distance, penalty_facetdiff, penalty_height, &
                                    penalty_overall)
                            
                                ! ---- set the 2.5 km location to be the one that has the smallest penalty
                        
                                dmin = 9999999.
                                ipmin = -99
                                jpmin = -99
                                DO jymin = jmin, jmax
                                    jdiff = jymin - jmin + 1
                                    DO ixmin = imin, imax
                                        idiff = ixmin - imin + 1
                                        IF (difference(idiff,jdiff) .lt. dmin) THEN
                                            dmin = difference(idiff, jdiff)
                                            ipmin = ixmin
                                            jpmin = jymin
                                        ENDIF
                                    END DO
                                END DO
                                xlocation_2p5(ix,jy,isupp) = ipmin
                                ylocation_2p5(ix,jy,isupp) = jpmin
                                DEALLOCATE (difference)
                                
                            END IF 
                            !PRINT *,'isupp = ',isupp,'. Minimum found at ',ipmin, jpmin
                        ELSE
                            xlocation_2p5(ix,jy,isupp) = -99
                            ylocation_2p5(ix,jy,isupp) = -99
                        ENDIF
                        
                    END DO  ! nsupp
                    
                ELSE ! inearest_ccpa(ix,jy) = -99, so we have not yet defined supplemental 
                     ! locations on the CCPA grid near to this 2.5-km gpoint
                    
                    ! ---- we need to do a more long-winded computational procedure here to 
                    !      search for supplemental locations.  Compute the penalty difference
                    !      for every grid point less than maxseparation distance 
                    
                    nproc = nproc + 1
                    pa_here(:) = panal_quantiles(ix, jy,:)
                    imin = max(1, ix-maxseparation)
                    imax = min(nx_2p5, ix+maxseparation)
                    jmin = max(1, jy-maxseparation)
                    jmax = min(ny_2p5, jy+maxseparation)
                    nxdiff = imax - imin +1
                    nydiff = jmax - jmin +1
                    ALLOCATE (difference(nxdiff,nydiff))
                    difference(:,:) = 99999.

                    CALL determine_penalties(ix, jy, nx_2p5, ny_2p5, nxdiff, &
                        nydiff, npct, imin, imax, jmin, jmax, maxseparation, &
                        ifacet_short_2p5, ifacet_medium_2p5, &
                        ifacet_long_2p5, panal_quantiles, ter_stddev, pa_stddev, &
                        distpenalty, facet_coeff, pa_coeff, terht_coeff, &
                        terrain_height, conusmask_2p5, imin_pa, imax_pa, &
                        stdmax_0p5, lons_2p5, lats_2p5, difference, penalty_analdiff, &
                        penalty_distance, penalty_facetdiff, penalty_height, &
                        penalty_overall)

                    DO isupp = 1, nsupp  ! number of supplemental locations

                        IF (isupp .eq. 1) THEN ! the first supplemental location
                                     !  is simply that orginal grid point
                                     
                            xlocation_2p5(ix,jy,isupp) = ix
                            ylocation_2p5(ix,jy,isupp) = jy

                            ! ---- don't consider points right around grid pt of interest,
                            !      too strong a correlation (want quasi-independent samples)

                            diffmin = 0.
                            CALL mask_around_thisgridpt(ix, jy, nx_2p5, ny_2p5, &
                                minseparation, lons_2p5, lats_2p5, maskout_2p5)

                        ELSE

				            ! ---- now find the analysis grid point with the next 
				            !      closest similarity. Set this as supplemental location 
				            !      point and then mask around this to eliminate nearby
				            !      points from any future consideration

                            minx1 = -99
                            miny1 = -99
                            diffmin = 999999.
                            DO jy2 = jmin, jmax  !1, nya
                                jdiff = jy2 - jmin + 1
                                DO ix2 = imin, imax   !1, nxa
                                    idiff = ix2 - imin + 1
                                    IF (difference(idiff,jdiff) .lt. diffmin .and. &
                                    maskout_2p5(ix2,jy2) .eq. 0 .and. conusmask_2p5(ix2,jy2) .gt. 0) THEN
                                        diffmin = difference(idiff,jdiff)
                                        minx1 = ix2
                                        miny1 = jy2
                                    END IF
                                END DO ! ix2
                            END DO  ! jy2
               
                            ! ---- finally, define the supplemental 
				            !      location on the forecast grid to be the 
                            !      grid point that had this closest fit as defined above.
                  
                            xlocation_2p5(ix,jy,isupp) = minx1
                            ylocation_2p5(ix,jy,isupp) = miny1

                            ! ---- make sure no other grid points very close to  
                            !      this new supplemental location are considered 
				            !      as another potential supplemental location.

                            CALL mask_around_thisgridpt(minx1, miny1, nx_2p5, ny_2p5, &
                                minseparation, lons_2p5, lats_2p5, maskout_2p5)

                        ENDIF ! isupp > 1

                        !penalty(ix,jy,isupp) = diffmin

                    END DO ! isupp
                    DEALLOCATE(difference)
                ENDIF ! IF (inearest_ccpa(ix,jy) .ne. -99) THEN
            ENDIF  ! conusmask_2p5
         END DO  ! ix
    END DO ! jy

    ! --- save out data for this month

    outfile = '/Users/thamill/refcst2/test/supplemental_locations_NDFD2p5_'//&
        cmonth//'.dat'
    PRINT *, 'writing to ', TRIM(outfile)
    OPEN (UNIT=1, FILE=outfile, STATUS='replace', FORM='unformatted')
    WRITE (1) nx_2p5, ny_2p5, nsupp, nx_ccpa, ny_ccpa, nsupp_ccpa
    WRITE (1) xlocation_2p5
    WRITE (1) ylocation_2p5
    WRITE (1) conusmask_2p5
    WRITE (1) lons_2p5
    WRITE (1) lats_2p5
    WRITE (1) xlocation_ccpa
    WRITE (1) ylocation_ccpa
    WRITE (1) conusmask_ccpa
    WRITE (1) lons_ccpa
    WRITE (1) lats_ccpa
    !WRITE (1) penalty
    CLOSE (1)

END DO  ! imonth
PRINT *,'Done'

CONTAINS                        

! =============================================================================

SUBROUTINE mask_around_thisgridpt(minx, miny, nx, ny, minseparation, &
    lons_2p5, lats_2p5, maskout)
! set maskout array around grid point of interest to 1 
!    to indicate nonconsideration of this in the future.
INTEGER, INTENT(IN) :: minx, miny, nx, ny, minseparation
REAL, INTENT(IN), DIMENSION(nx,ny) :: lats_2p5, lons_2p5
INTEGER(2), INTENT(INOUT), DIMENSION(nx,ny) :: maskout

imin_m = MAX(1,minx-minseparation)
imax_m = MIN(nx,minx+minseparation)
jmin_m = MAX(1,miny-minseparation)
jmax_m = MIN(ny,miny+minseparation)

DO j = jmin_m, jmax_m
    DO i = imin_m, imax_m
        CALL haversine (lats_2p5(minx,miny), lons_2p5(minx,miny), &
            lats_2p5(i,j), lons_2p5(i,j), dist) ! in km
        IF (dist .le. REAL(minseparation)*2.5) maskout(i,j) = 1
    END DO
END DO
RETURN
END SUBROUTINE mask_around_thisgridpt

! =============================================================================

SUBROUTINE compute_imin_imax_of_CDF (npct, pa_here, imin_c, imax_c)

! ---- when comparing CDFs at different grid points, we don't want 
!      the CDF comparison to either focus on the zero values 
!      or the extreme values of the distribution, which can be sensitive 
!      to outliers.  Return imin and imax, the grid indices bounding 
!      the computation of the CDFs.

INTEGER, INTENT(IN) :: npct
REAL, INTENT(IN), DIMENSION(npct) :: pa_here
INTEGER, INTENT(OUT) :: imin_c, imax_c

imin_c = 1
imax_c = npct
DO i = 1,npct
   IF (pa_here(i) .gt. 0.) THEN
      imin_c = i
      GOTO 1231
   ENDIF
END DO

1231 CONTINUE
DO i = npct,1,-1
   IF (pa_here(i) .le. 0.95) THEN
      imax_c = i
      GOTO 1241
   ENDIF
END DO

1241 CONTINUE

IF (imax_c .le. imin_c) THEN
   imax_c = npct
   imin_c = 1
ENDIF

IF (imin_c .lt. 1 .or. imax_c .gt. npct) THEN
   PRINT *,'imax_c, imin_c = ',imax_c, imin_c
   PRINT *,'pa_here = ',pa_here
ENDIF
RETURN
END SUBROUTINE compute_imin_imax_of_CDF

! ==============================================================================

SUBROUTINE determine_penalties(ix, jy, nx_2p5, ny_2p5, nxdiff, &
    nydiff, npct, imin, imax, jmin, jmax, maxseparation, &
    ifacet_short_2p5, ifacet_medium_2p5, &
    ifacet_long_2p5, panal_quantiles, ter_stddev, pa_stddev, &
    distpenalty, facet_coeff, pa_coeff, terht_coeff, &
    terrain_height, conusmask_2p5, imin_pa, imax_pa, stdmax_0p5, &
    lons_2p5, lats_2p5, difference, penalty_analdiff, &
    penalty_distance, penalty_facetdiff, penalty_height, &
    penalty_overall)
        
INTEGER, INTENT(IN) :: ix, jy, nx_2p5, ny_2p5, npct, imin, imax
INTEGER, INTENT(IN) :: jmin, jmax, maxseparation

INTEGER, INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: ifacet_short_2p5  ! facet associated with fine smoothing
INTEGER, INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: ifacet_medium_2p5 ! ... medium ....
INTEGER, INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: ifacet_long_2p5 ! ... long ...
INTEGER(2), INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: conusmask_2p5
INTEGER, INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: imin_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt
INTEGER, INTENT(IN), DIMENSION(nx_2p5,ny_2p5) :: imax_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5,npct) :: panal_quantiles
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5)   :: pa_stddev
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5)   :: terrain_height
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5)   :: ter_stddev
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5)   :: lons_2p5
REAL, INTENT(IN), DIMENSION(nx_2p5,ny_2p5)   :: lats_2p5

REAL, INTENT(IN) :: distpenalty
REAL, INTENT(IN) :: facet_coeff
REAL, INTENT(IN) :: pa_coeff
REAL, INTENT(IN) :: terht_coeff
REAL, INTENT(IN) :: stdmax_0p5

REAL, INTENT(OUT), DIMENSION(nxdiff, nydiff) :: difference
REAL, INTENT(OUT) :: penalty_facetdiff
REAL, INTENT(OUT) :: penalty_analdiff
REAL, INTENT(OUT) :: penalty_distance
REAL, INTENT(OUT) :: penalty_height
REAL, INTENT(OUT) :: penalty_overall

REAL facet_avg
REAL height_diff
REAL height_diff_penalty
REAL d2
REAL dist
REAL padiff

REAL, DIMENSION(npct) :: pa_here, pa_there

REAL smoothness_ratio

INTEGER idiff_short, idiff_medium, idiff_long

pa_here(:) = panal_quantiles(ix, jy,:)

difference(:,:) = 9999999.

DO jy2 = jmin, jmax
    jdiff = jy2 - jmin + 1
    DO ix2 = imin, imax
        idiff = ix2 - imin +1

        ! ---- only consider this grid point if inside area with data, 
        !      and is less than max separation, computed in 1/8-degree 
        !      grid points.  Since there is less trustworthy data outside 
        !      CONUS (lots of problems over ocean, over Canada at this  
        !      point, supplement data only with points from inside CONUS

        CALL haversine (lats_2p5(ix,jy), lons_2p5(ix,jy), &
            lats_2p5(ix2,jy2), lons_2p5(ix2,jy2), dist)
        penalty_distance = dist*distpenalty 
        !PRINT *,'ix2, jy2, dist, maxsep = ', ix2, jy2, dist, REAL(maxseparation)*2.5

        IF (conusmask_2p5(ix2,jy2) .gt. 0 .and. dist .le. REAL(maxseparation)*2.5) THEN

            ! --- quantify the difference between analyzed cumulative  
            !     distributions at the two grid points 
            
            pa_there(:) = panal_quantiles(ix2, jy2,:)
            padiff = SUM(ABS(pa_here(imin_pa(ix,jy):imax_pa(ix,jy))- &
                pa_there(imin_pa(ix,jy):imax_pa(ix,jy)))) / &
                REAL(imax_pa(ix,jy)-imin_pa(ix,jy)+1)
            penalty_analdiff = pa_coeff*padiff / pa_stddev(ix,jy)

            ! ---- quantify the difference in facets at each of 
            !      the three scales and convert this to a penalty 
		    !      for average facet differences

            idiff_short = &
                MIN(IABS(ifacet_short_2p5(ix,jy) - ifacet_short_2p5(ix2,jy2)) , &
                IABS(ifacet_short_2p5(ix,jy) - (ifacet_short_2p5(ix2,jy2)-8) ), &
                IABS(ifacet_short_2p5(ix,jy) - (ifacet_short_2p5(ix2,jy2)+8)))
            idiff_medium = &
                MIN(IABS(ifacet_medium_2p5(ix,jy) - ifacet_medium_2p5(ix2,jy2)) , &
                IABS(ifacet_medium_2p5(ix,jy) - (ifacet_medium_2p5(ix2,jy2)-8) ), &
                IABS(ifacet_medium_2p5(ix,jy) - (ifacet_medium_2p5(ix2,jy2)+8)))
            idiff_long = &
                MIN(IABS(ifacet_long_2p5(ix,jy) - ifacet_long_2p5(ix2,jy2)) , &
                IABS(ifacet_long_2p5(ix,jy) - (ifacet_long_2p5(ix2,jy2)-8) ), &
                IABS(ifacet_long_2p5(ix,jy) - (ifacet_long_2p5(ix2,jy2)+8)))
            facet_avg = REAL(idiff_short + idiff_medium + idiff_long) / 3.

            ! ---- facets really shouldn't matter much if the terrain  
            !      is smooth. calculate a weight between 0 and 1. 

            smoothness_ratio = ter_stddev(ix,jy)**0.5 / stdmax_0p5 

            ! ---- calculate the final facet weight

            IF (facet_avg .le. 4.0) THEN
                penalty_facetdiff = &
                facet_coeff*smoothness_ratio*(facet_avg/4.)
            ELSE
                penalty_facetdiff = facet_coeff*smoothness_ratio
            ENDIF

            ! ---- quantify difference in terrain ht 
		    !      and turn into terrain penalty

            height_diff = ABS(terrain_height(ix,jy) - terrain_height(ix2,jy2))
            height_diff_penalty = 1.0 - exp(-height_diff/(50.**2))
            penalty_height = terht_coeff*height_diff_penalty

            ! --- using information on the difference between forecast CDFs
            !     (pfdiff), and F/O regression relationship information
            !     (b0_mean, b1_mean), determine the single number that 
            !     quantifies the overall strength of relationship between
            !     forecast characteristics at these two grid points.

            difference(idiff, jdiff) = penalty_analdiff + &
                penalty_distance + penalty_facetdiff + penalty_height
            
        ELSE
            difference(idiff, jdiff) = 9999999.
        ENDIF
        penalty_overall = difference(idiff, jdiff)
    END DO ! ix2
END DO ! jy2
RETURN
END SUBROUTINE determine_penalties

! ======================================================================

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
REAL, PARAMETER :: deg_to_rad = 3.1415926/180.
    ! exploit intrinsic atan to generate pi/180 runtime constant
REAL :: rad

rad = degree*deg_to_rad
END SUBROUTINE to_radian


END PROGRAM create_supp_locns_ndfd2p5
