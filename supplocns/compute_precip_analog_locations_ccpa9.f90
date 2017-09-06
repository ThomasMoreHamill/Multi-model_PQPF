PROGRAM compute_precip_analog_locations_ccpa9

! purpose:  For every (i,j) in on the 1/8-degree CCPA grid with valid data, 
!    determine a set of supplemental data locations.
!    Base the supplemental locations on the similarity of analyzed CDFs 
!    as well as the local terrain heights and facets.  Do this separately 
!    for every month of the year. We want grid points that will have 
!    somewhat independent data, so make sure that the set of supplemental 
!    locations are not too close together.
! 
USE netcdf

PARAMETER (nyears = 12)    ! 2002-2013 here for the CCPA data
PARAMETER (ndaysub = 91)   ! number of days of training data for each yr, 
                           ! 15th of month +/- 45 days
PARAMETER (nxa_small = 464)! CCPA 1/8-degree precip analysis x grid dimensions, 
                           ! surrounding CONUS
PARAMETER (nya_small = 224)  ! y grid dimension
PARAMETER (nxa = 515)      ! CCPA 1/8-degree grid x dimensions for expanded grid.
PARAMETER (nya = 262)      ! CCPA 1/8-degree grid y dimensions for expanded grid
PARAMETER (ioffset = 21)   ! shift in x-dir for CCPA small grid within larger grid
PARAMETER (joffset = 22)   ! shift in y-dir for CCPA small grid within larger grid
PARAMETER (nmonths = 12)   ! # months of the year
PARAMETER (nmaxlist = 16)  ! max number of CCPA grid pts assoc'd with 
	                       ! any forecast grid pt.
PARAMETER (nanalogs = 100) ! max number of CCPA grid supplemental locations to be 
                           ! stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 6) ! user-defined minimum separation between 
                           ! analysis grid locations and supp location (in grid pts)
PARAMETER (maxseparation = 250) ! user-defined maximum separation between analysis grid 
                           ! location and supplemental location (in grid pts)
PARAMETER (npct = 107)     ! number of quantiles of CDF to compute at
PARAMETER (pa_coeff = 0.1) ! coefficient to apply to analysis CDF differences
PARAMETER (facet_coeff = 0.6) ! coefficient to apply to terrain facet differences
PARAMETER (terht_coeff = 0.4)  ! coefficient to apply to terrain height differences
PARAMETER (distpenalty = .001)  ! coefficient to apply to distance 
PARAMETER (dx = 52083.)    ! forecast grid spacing at equator

INTEGER*2, DIMENSION(nxa,nya) :: conusmask  ! 1/0 mask for whether in CONUS
INTEGER*2, DIMENSION(nxa,nya) :: maskout  ! array used to block out nearby points !
                           ! from further consideration
INTEGER, DIMENSION(nxa,nya,nanalogs) :: xlocationa ! for each forecast point,  
                           ! a list of which other forecast points
INTEGER, DIMENSION(nxa,nya,nanalogs) :: ylocationa ! have the closest climatologies 
                           ! and forecast-obs relationship
INTEGER minx
INTEGER miny
INTEGER*2, DIMENSION(nxa,nya) :: imin_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt
INTEGER*2, DIMENSION(nxa,nya) :: imax_pa  ! the index of minimum for which to compute 
                                          ! differences in CDFs for a grid pt

INTEGER*2, DIMENSION(nxa,nya) :: ifacet_short  ! facet associated with fine smoothing
INTEGER*2, DIMENSION(nxa,nya) :: ifacet_medium ! ... medium ....
INTEGER*2, DIMENSION(nxa,nya) :: ifacet_long ! ... long ...
REAL, DIMENSION(nxa,nya)      :: difference
REAL, DIMENSION(nxa,nya)      :: difference_pa
REAL, DIMENSION(nxa,nya)      :: difference_facet
REAL, DIMENSION(nxa,nya)      :: difference_ter
REAL, DIMENSION(nxa,nya)      :: difference_dist
REAL, DIMENSION(nxa,nya)      :: lonsa
REAL, DIMENSION(nxa,nya)      :: latsa
REAL, DIMENSION(nxa,nya,npct) :: panal_quantiles
REAL, DIMENSION(npct)         :: pctv
REAL, DIMENSION(npct)         :: pa_here
REAL, DIMENSION(npct)         :: pa_there
REAL, DIMENSION(nxa,nya)      :: pa_stddev
REAL, DIMENSION(nxa,nya)      :: pa_mean
REAL, DIMENSION(nxa,nya,nanalogs) :: penalty  ! penalty function 
!                                               associated with the supp location
REAL, DIMENSION(nxa,nya,nanalogs) :: penalty_pa  ! penalty for CDF differences
REAL, DIMENSION(nxa,nya,nanalogs) :: penalty_ter  ! penalty for terrain ht difference
REAL, DIMENSION(nxa,nya,nanalogs) :: penalty_facet  ! penalty for facet differences
REAL, DIMENSION(nxa,nya,nanalogs) :: penalty_dist  ! penalty for distance differences
REAL, DIMENSION(nxa,nya)      :: terrain
REAL, DIMENSION(nxa,nya)      :: ter_mean
REAL, DIMENSION(nxa,nya)      :: ter_stddev

REAL time1, time2

REAL diffmin

CHARACTER*120 outfile
CHARACTER*3, DIMENSION(12) :: cmonth
CHARACTER*2 cmonthin
CHARACTER*3 csupp

DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun', &
             'Jul','Aug','Sep','Oct','Nov','Dec'/

CALL getarg(1, cmonthin)
READ (cmonthin,'(i2)') imonthin

IF (nanalogs .ge. 100) THEN
   WRITE (csupp,'(i3)') nanalogs
ELSE
   WRITE (csupp,'(i2)') nanalogs
   csupp = '0' // csupp(1:2)
ENDIF
print *,'csupp = ',csupp

! ---- for each month, load the quantile 
!      information and rank correlation previously
!      computed, and then determine analogs

xlocation = -99  ! initialize all locations to missing value.
ylocation = -99
xlocationa = -99
ylocationa = -99

penalty = 0.
penalty_pa = 0.
penalty_ter = 0.
penalty_facet = 0.

! ---- read in terrain facet information for smoothing 
!      at short, intermediate, and larger scales

CALL read_facets(nxa, nya, ifacet_short, & 
	ifacet_medium, ifacet_long) ! 0.5, 3.0, 10.0

!DO imonth = imonthin, imonthin
DO imonth = 1,12

   ! ---- load the CDF data 

   PRINT *,'calling  load_precipquantiles_ccpa_analonly_v9'
   CALL load_precipquantiles_ccpa_analonly_v9 &
   	  (imonth, nxa, nya, npct, nxa_small, nya_small, ioffset, joffset, &
      conusmask, lonsa, latsa, pctv, panal_quantiles, terrain)
   PRINT *, 'pctv = ',pctv
   
   print *,'min, max conusmask = ',minval(conusmask), maxval(conusmask)

   ! ---- precompute the indices of the minimum and 
   !      maximum indices for quantiles where we will
   !      consider evaluating the CDF at

   PRINT *,'calling compute_imin_imax_of_CDF'
   DO jya = 1, nya
      DO ixa = 1, nxa
         pa_here(:) = panal_quantiles(ixa, jya,:)
         CALL compute_imin_imax_of_CDF (npct, pctv, &
		 	pa_here, imin_pa(ixa,jya), imax_pa(ixa,jya))
      END DO
   END DO

   ! ---- loop thru grid points and precompute statistics for mean, std dev
   !      of analysis CDF and terrain

   PRINT *,'pre-processing to determine local standard devs of CDF and terrain'
   DO ixa = 1, nxa
      IF (MOD(ixa,10) .eq. 0) PRINT *,'processing  ixa = ',ixa,' of ',nxa

      ! --- only bother searching in a box +/- 10 grid points to limit computations.

      ixmin = MAX(1,ixa-10)
      ixmax = MIN(nxa,ixa+10)
      DO jya = 1, nya
         jymin = MAX(1,jya-10)
         jymax = MIN(nya,jya+10)
         sum_pa2 = 0.0
         sum_pa  = 0.0
         sum_t   = 0.0
         sum_t2  = 0.0
         nsamps  = 0
         ter_here = 0.
         ter_there = 0.
         IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point has valid data
            pa_here(:) = panal_quantiles(ixa, jya,:)
            ter_here = terrain(ixa,jya)
            DO ix2 = ixmin, ixmax
               DO jy2 = jymin, jymax
                  IF (conusmask(ix2,jy2).gt. 0) THEN
                     pa_there(:) = panal_quantiles(ix2, jy2,:)
                     ter_there = terrain(ix2,jy2)
					 
                     ! --- quantify the difference between forecast distributions at 
                     !     the two grid points
					 
                     padiff = SUM(ABS(pa_here(imin_pa(ixa,jya):imax_pa(ixa,jya)) - &
                        pa_there(imin_pa(ixa,jya):imax_pa(ixa,jya)))) / &
                        REAL(imax_pa(ixa,jya)-imin_pa(ixa,jya)+1)
                     sum_t2     = sum_t2     + ter_there**2
                     sum_t      = sum_t      + ter_there
                     sum_pa2    = sum_pa2    + padiff**2
                     sum_pa     = sum_pa     + padiff
                     nsamps  = nsamps  + 1
                  ENDIF
               END DO
            END DO

            ! --- now calculate std deviation by shortcut (Wilks 2006, eq. 3.20)

            pa_mean(ixa,jya) = sum_pa / REAL(nsamps)
            pa_stddev(ixa,jya) = &
				SQRT((sum_pa2-REAL(nsamps)*pa_mean(ixa,jya)**2) / &
				REAL(nsamps-1))
            ter_mean(ixa,jya) = sum_t  / REAL(nsamps)
            ter_stddev(ixa,jya) = &
				SQRT((sum_t2-REAL(nsamps)*ter_mean(ixa,jya)**2) / &
				REAL(nsamps-1))
         ELSE
            pa_mean(ixa,jya)       = -99.99
            pa_stddev(ixa,jya)     = -99.99
            ter_mean(ixa,jya)      = -99.99
            ter_stddev(ixa,jya)    = -99.99
         ENDIF

      END DO ! jya
   END DO    ! ixa

   PRINT *,' ------  pa_mean ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(pa_mean(i,j),i=25,45)
      123 format(i3,1x,21(f6.2,1x))
   END DO

   PRINT *,' ------  pa_stddev ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(pa_stddev(i,j),i=25,45)
   END DO

   PRINT *,' ------  ter_mean ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(ter_mean(i,j),i=25,45)
   END DO

   PRINT *,' ------  ter_stddev ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(ter_stddev(i,j),i=25,45)
   END DO

   stdmax = MAXVAL(ter_stddev)
   stdmax_0p5 = stdmax**0.5

   ! ---- find locations of supplemental locations 
   !      that are the best fit for each grid point

   CALL cpu_time(time1)
   PRINT *, 'finding supplemental locations'
   nproc = 0
   DO ixa = 1, nxa
   !DO ixa = 35, 235, 20 !35,nya/2 near SFO
      CALL cpu_time(time2)
      IF (MOD(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa,&
	  	' time elapsed, delta: ',time2, time2-time1,'nproc = ',&
		nproc,' total pts = ',SUM(REAL(conusmask))
      DO jya = 1, nya

         maskout(:,:) = 0
         difference(:,:) = 99999999.
         difference_pa(:,:) = 99999999.
         difference_ter(:,:) = 99999999.
         difference_facet(:,:) = 99999999.
         difference_dist(:,:) = 99999999.

         IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point inside area with CCPA data

            nproc = nproc+1
            pa_here(:) = panal_quantiles(ixa, jya,:)
            imin = max(1,ixa-maxseparation)
            imax = min(nxa,ixa+maxseparation)
            jmin = max(1,jya-maxseparation)
            jmax = min(nya,jya+maxseparation)

            DO ix2 = imin, imax
               DO jy2 = jmin, jmax

                  ! ---- only consider this grid point if inside area with data, 
                  !      and is less than max separation, computed in 1/8-degree 
                  !      grid points.  Since there is less trustworthy data outside 
                  !      CONUS (lots of problems over ocean, over Canada at this  
                  !      point, supplement data only with points from inside CONUS

                  d2 = REAL((ixa-ix2)**2 + (jya-jy2)**2)
                  dist = SQRT(d2)

                  IF (conusmask(ix2,jy2) .gt. 0 .and. maskout(ix2,jy2) &
					 .eq. 0 .and. dist .le. maxseparation) THEN

                     ! --- quantify the difference between analyzed cumulative  
                     !     distributions at the two grid points 

                     pa_there(:) = panal_quantiles(ix2, jy2,:)
                     padiff = SUM(ABS(pa_here(imin_pa(ixa,jya):imax_pa(ixa,jya))- &
                          pa_there(imin_pa(ixa,jya):imax_pa(ixa,jya)))) / &
                          REAL(imax_pa(ixa,jya)-imin_pa(ixa,jya)+1)

                     ! ---- quantify the difference in facets at each of 
                     !      the three scales and convert this to a penalty 
					 !      for average facet differences

                     idiff_short = &
                        MIN(ABS(ifacet_short(ixa,jya) - ifacet_short(ix2,jy2)) , &
                        ABS(ifacet_short(ixa,jya) - (ifacet_short(ix2,jy2)-8) ), &
                        ABS(ifacet_short(ixa,jya) - (ifacet_short(ix2,jy2)+8)))
                     idiff_medium = &
                        MIN(ABS(ifacet_medium(ixa,jya) - ifacet_medium(ix2,jy2)) , &
                        ABS(ifacet_medium(ixa,jya) - (ifacet_medium(ix2,jy2)-8) ), &
                        ABS(ifacet_medium(ixa,jya) - (ifacet_medium(ix2,jy2)+8)))
                     idiff_long = &
                        MIN(ABS(ifacet_long(ixa,jya) - ifacet_long(ix2,jy2)) , &
                        ABS(ifacet_long(ixa,jya) - (ifacet_long(ix2,jy2)-8) ), &
                        ABS(ifacet_long(ixa,jya) - (ifacet_long(ix2,jy2)+8)))
                     facet_avg = REAL(idiff_short + idiff_medium + idiff_long) / 3.

                     ! ---- facets really shouldn't matter much if the terrain  
                     !      is smooth. calculate a weight between 0 and 1. 

                     smoothness_ratio = ter_stddev(ixa,jya)**0.5 / stdmax_0p5 

                     ! ---- calculate the final facet weight

                     IF (facet_avg .le. 4.0) THEN
                        penalty_facetdiff = &
						   facet_coeff*smoothness_ratio*(facet_avg/4.)
                     ELSE
                        penalty_facetdiff = facet_coeff*smoothness_ratio
                     ENDIF

                     ! ---- quantify difference in terrain ht 
					 !      and turn into terrain penalty

                     height_diff = ABS(terrain(ixa,jya) - terrain(ix2,jy2))
                     height_diff_penalty = 1.0 - exp(-height_diff/(50.**2))

                     ! --- using information on the difference between forecast CDFs
                     !     (pfdiff), and F/O regression relationship information
                     !     (b0_mean, b1_mean), determine the single number that 
                     !     quantifies the overall strength of relationship between
                     !     forecast characteristics at these two grid points.

                     difference(ix2,jy2) = pa_coeff*padiff / pa_stddev(ixa,jya) + &
                          dist*distpenalty + penalty_facetdiff*facet_coeff + &
                          terht_coeff*height_diff_penalty
                     difference_pa(ix2,jy2) = pa_coeff*padiff / pa_stddev(ixa,jya) 
                     difference_ter(ix2,jy2) = terht_coeff*height_diff_penalty
                     difference_facet(ix2,jy2) =  penalty_facetdiff*facet_coeff 
                     difference_dist(ix2,jy2) = dist*distpenalty
                  ELSE
                     difference(ix2,jy2) = 9999999.
                     difference_pa(ix2,jy2) = 9999999.
                     difference_ter(ix2,jy2) = 9999999.
                     difference_facet(ix2,jy2) = 9999999.
                     difference_dist(ix2,jy2) = 9999999.

                  ENDIF
               END DO ! jy2
            END DO ! ix2

            DO iana = 1, nanalogs  ! number of supplemental locations

               IF (iana .eq. 1) THEN ! the first supplemental location
                                     !  is simply that orginal grid point
                  xlocationa(ixa,jya,iana) = ixa
                  ylocationa(ixa,jya,iana) = jya

                  ! ---- don't consider points right around grid pt of interest,
                  !      too strong a correlation (want quasi-independent samples)

                  diffmin = 0.
                  diffmin_pa = 0.
                  diffmin_ter = 0.
                  diffmin_facet = 0.
                  diffmin_dist = 0.
                  diffmin = 0.
                  CALL mask_around_thisgridpt(ixa, jya, nxa, nya, &
				  	  minseparation, maskout)

               ELSE

				  ! ---- now find the analysis grid point with the next 
				  !      closest similarity. Set this as supplemental location 
				  !      point and then mask around this to eliminate nearby
				  !      points from any future consideration

                  minx = -99
                  miny = -99
                  diffmin = 999999.
                  diffmin_pa = 999999.
                  diffmin_ter = 999999.
                  diffmin_facet = 999999.
                  diffmin_dist = 999999.
                  diffmin = 999999.
                  DO ix2 = imin, imax   !1, nxa
                     DO jy2 = jmin, jmax  !1, nya

                        IF (difference(ix2,jy2) .lt. diffmin .and. &
						maskout(ix2,jy2) .eq. 0 .and. &
						conusmask(ix2,jy2) .gt. 0) THEN
                           diffmin_pa = difference_pa(ix2,jy2)
                           diffmin_ter = difference_ter(ix2,jy2)
                           diffmin_facet = difference_facet(ix2,jy2)
                           diffmin_dist = difference_dist(ix2,jy2)
                           diffmin = difference(ix2,jy2)
                           minx = ix2
                           miny = jy2
                        END IF
                     END DO ! jy2
                  END DO  ! ix2
               
                  ! ---- finally, define the supplemental 
				  !      location on the forecast grid to be the 
                  !      grid point that had this closest fit as defined above.
                  
                  xlocationa(ixa,jya,iana) = minx
                  ylocationa(ixa,jya,iana) = miny

                  ! ---- make sure no other grid points very close to  
                  !      this new supplemental location are considered 
				  !      as another potential supplemental location.

                  CALL mask_around_thisgridpt(minx, miny, nxa, nya, &
				  	  minseparation, maskout)

               ENDIF ! iana > 1

               penalty(ixa,jya,iana) = diffmin
               penalty_pa(ixa,jya,iana) = diffmin_pa
               penalty_ter(ixa,jya,iana) = diffmin_ter
               penalty_facet(ixa,jya,iana) = diffmin_facet
               penalty_dist(ixa,jya,iana) = diffmin_dist

            END DO ! iana
         ENDIF   ! conusmask
      END DO  ! jya
   END DO ! ixa

   ! --- save out data for this month

   outfile = '/Users/thamill/refcst2/test/ppn_analog_locns_ccpa9_'//&
   	   cmonth(imonth)//'.dat'
   PRINT *, 'writing to ', TRIM(outfile)
   OPEN (UNIT=1, FILE=outfile, STATUS='replace', FORM='unformatted')
   WRITE (1) nxf, nyf, nxa, nya, nanalogs
   WRITE (1) xlocationa
   WRITE (1) ylocationa
   WRITE (1) conusmask
   WRITE (1) lonsa
   WRITE (1) latsa
   WRITE (1) penalty
   WRITE (1) penalty_pa
   WRITE (1) penalty_ter
   WRITE (1) penalty_facet
   WRITE (1) penalty_dist
   CLOSE (1)

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
DO i = imin_m, imax_m
   DO j = jmin_m, jmax_m
      dist2 = REAL((minx-i)**2 + (miny-j)**2)
      dist = SQRT(dist2)
      IF (dist .le. minseparation) maskout(i,j) = 1
   END DO
END DO
RETURN
END SUBROUTINE mask_around_thisgridpt

! =============================================================================

SUBROUTINE compute_imin_imax_of_CDF (npct, pctv, pa_here, imin_c, imax_c)

! ---- when comparing CDFs at different grid points, we don't want 
!      the CDF comparison to either focus on the zero values 
!      or the extreme values of the distribution, which can be sensitive 
!      to outliers.  Return imin and imax, the grid indices bounding 
!      the computation of the CDFs.

INTEGER, INTENT(IN) :: npct
REAL, INTENT(IN), DIMENSION(npct) :: pctv
REAL, INTENT(IN), DIMENSION(npct) :: pa_here
INTEGER*2, INTENT(OUT) :: imin_c, imax_c

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
   IF (pctv(i) .le. 0.95) THEN
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
   PRINT *,'pctv = ',pctv
ENDIF
RETURN
END SUBROUTINE compute_imin_imax_of_CDF

! ==============================================================================

END PROGRAM compute_precip_analog_locations_ccpa9
