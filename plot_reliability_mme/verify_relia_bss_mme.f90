!  f2py -c -m verify_relia_bss_mme verify_relia_bss_mme.f90 read_prob_forecasts_postproc.f90 ran3.f sort.f check.f90 -I/usr/local/gfortran/include -L/usr/local/gfortran/lib  -L/opt/local/lib -lnetcdff

SUBROUTINE verify_relia_bss_mme(cleade, nclasses, rthresh, &
    date_list_anal, apcp_anal_t, nxa, nya, ndates, &
    relia_NCEP, relia_NCEP_05, relia_NCEP_95, frequse_NCEP, bss_NCEP, &
    relia_CMC, relia_CMC_05, relia_CMC_95, frequse_CMC, bss_CMC, &
    relia_ECMWF, relia_ECMWF_05, relia_ECMWF_95, frequse_ECMWF, bss_ECMWF, &
    relia_MME, relia_MME_05, relia_MME_95, frequse_MME, bss_MME)

! purpose:  generate reliability information, frequency of usage, and Brier
!  skill scores for post-processed NCEP, CMC, ECMWF, and MME forecasts

PARAMETER (nresa = 1000)
PARAMETER (nthreshes = 7)
CHARACTER*3, INTENT(IN) :: cleade
INTEGER, INTENT(IN) :: nclasses, nxa, nya, ndates
REAL, INTENT(IN) :: rthresh ! precip threshold amount
REAL, DIMENSION(nxa,nya,ndates), INTENT(IN) :: apcp_anal_t
CHARACTER*10, DIMENSION(ndates), INTENT(IN) :: date_list_anal


REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_NCEP, relia_NCEP_05, relia_NCEP_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_CMC, relia_CMC_05, relia_CMC_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_ECMWF, relia_ECMWF_05, relia_ECMWF_95	
REAL, DIMENSION(nclasses), INTENT(OUT) :: relia_MME, relia_MME_05, relia_MME_95
REAL, DIMENSION(nclasses), INTENT(OUT) :: frequse_NCEP, frequse_CMC, &
    frequse_ECMWF, frequse_MME
    
REAL, INTENT(OUT) :: bss_NCEP, bss_CMC, bss_ECMWF, bss_MME
REAL*4, DIMENSION(ndates) :: bss_NCEP_daily
REAL*4, DIMENSION(ndates) :: bss_CMC_daily
REAL*4, DIMENSION(ndates) :: bss_ECMWF_daily
REAL*4, DIMENSION(ndates) :: bss_MME_daily

REAL, DIMENSION(nthreshes) :: pthreshes

! f2py intent(in) nxa, nya, ndates, cleade, nclasses 
! f2py intent(in) rthresh, apcp_anal_t, date_list_anal
! f2py depend(ndates) date_list_anal
! f2py depend(nxa,nya,ndates) apcp_anal_t
! f2py intent(out) relia_NCEP, relia_NCEP_05, relia_NCEP_95
! f2py intent(out) relia_CMC, relia_CMC_05, relia_CMC_95
! f2py intent(out) relia_ECMWF, relia_ECMWF_05, relia_ECMWF_95
! f2py intent(out) relia_MME, relia_MME_05, relia_MME_95
! f2py intent(out) frequse_NCEP, frequse_CMC, frequse_ECMWF, frequse_MME
! f2py intent(out) bss_NCEP, bss_CMC, bss_ECMWF, bss_MME
! f2py depend(nclasses) relia_NCEP, relia_NCEP_05, relia_NCEP_95
! f2py depend(nclasses) relia_CMC, relia_CMC_05, relia_CMC_95
! f2py depend(nclasses) relia_ECMWF, relia_ECMWF_05, relia_ECMWF_95
! f2py depend(nclasses) relia_MME, relia_MME_05, relia_MME_95
! f2py depend(nclasses) frequse_NCEP, frequse_CMC
! f2py depend(nclasses) frequse_ECMWF, frequse_MME
! f2py depend(nthreshes) bss_NCEP, bss_CMC, bss_ECMWF, bss_MME

! --- now local variables

INTEGER*2, DIMENSION(nxa,nya) :: conusmask

REAL, DIMENSION(nxa,nya) :: climo_prob
REAL, DIMENSION(nxa,nya) :: rlonsa
REAL, DIMENSION(nxa,nya) :: rlatsa
REAL, DIMENSION(nxa,nya) :: prob_forecast_NCEP
REAL, DIMENSION(nxa,nya) :: prob_forecast_ECMWF
REAL, DIMENSION(nxa,nya) :: prob_forecast_CMC
REAL, DIMENSION(nxa,nya) :: prob_forecast_MME
REAL, DIMENSION(nxa,nya) :: prob_forecast

REAL*8, DIMENSION(0:nclasses-1,2) :: contab
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_sum

REAL*8, DIMENSION(0:nclasses-1,2) :: contab_NCEP
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_CMC
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_ECMWF
REAL*8, DIMENSION(0:nclasses-1,2) :: contab_MME

REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_NCEP_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_CMC_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_ECMWF_daily
REAL*8, DIMENSION(ndates, 0:nclasses-1,2) :: contab_MME_daily

REAL*8, DIMENSION(0:nclasses-1) :: relia, relia_05, relia_95
REAL*8, DIMENSION(nresa, 0:nclasses-1) :: relia_resa
REAL*8, DIMENSION(0:nclasses-1) :: frequse

REAL*8 :: bs
REAL*8 :: bs_climo

REAL*8 :: bs_NCEP
REAL*8 :: bs_CMC
REAL*8 :: bs_ECMWF
REAL*8 :: bs_MME

REAL*8, DIMENSION(ndates) :: bs_daily
REAL*8, DIMENSION(ndates) :: bs_climo_daily
REAL*8, DIMENSION(ndates) :: bs_NCEP_daily
REAL*8, DIMENSION(ndates) :: bs_CMC_daily
REAL*8, DIMENSION(ndates) :: bs_ECMWF_daily
REAL*8, DIMENSION(ndates) :: bs_MME_daily

REAL, DIMENSION(nresa) :: rsamps

CHARACTER*120 infile
CHARACTER*10 cyyyymmddhh
CHARACTER*3 clead_use
CHARACTER*5 cmodel



contab_NCEP = 0.
contab_CMC = 0.
contab_ECMWF = 0.
contab_MME = 0.

contab_NCEP_daily = 0.
contab_CMC_daily = 0.
contab_ECMWF_daily = 0.
contab_NCEP_daily = 0.

bs_climo = 0.

bs_NCEP = 0.
bs_CMC = 0.
bs_ECMWF = 0.
bs_MME = 0.

bs_NCEP_daily(:) = 0.0
bs_CMC_daily(:) = 0.0
bs_ECMWF_daily(:) = 0.0
bs_MME_daily(:) = 0.0

pid180 = 3.1415926/180.
clead_use = cleade
IF (clead_use .eq. '012') clead_use = '12'
IF (clead_use .eq. '024') clead_use = '24'
IF (clead_use .eq. '036') clead_use = '36'
IF (clead_use .eq. '048') clead_use = '48'
IF (clead_use .eq. '060') clead_use = '60'
IF (clead_use .eq. '072') clead_use = '72'
IF (clead_use .eq. '084') clead_use = '84'
IF (clead_use .eq. '096') clead_use = '96'

PRINT *,'max(apcp_anal_t) = ', maxval(apcp_anal_t)

! ---- loop thru days and verify

DO idate = 1, ndates

    cyyyymmddhh = date_list_anal(idate) 

    ! ---- read in the forecast for this particular initial date/time and forecast lead and resolution

    cmodel = 'NCEP'
    CALL read_prob_forecasts_postproc(nxa, nya, cyyyymmddhh, cmodel, clead_use, rthresh,&
        nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
        prob_forecast_NCEP, climo_prob)

    cmodel = 'CMC'
    CALL read_prob_forecasts_postproc(nxa, nya, cyyyymmddhh, cmodel, clead_use, rthresh,&
        nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
        prob_forecast_CMC, climo_prob)
            
    cmodel = 'ECMWF'
    CALL read_prob_forecasts_postproc(nxa, nya, cyyyymmddhh, cmodel, clead_use, rthresh,&
        nthreshes, pthreshes, conusmask, rlonsa, rlatsa, &
        prob_forecast_ECMWF, climo_prob)

    pmax_NCEP = MAXVAL(prob_forecast_NCEP)
    pmax_CMC = MAXVAL(prob_forecast_CMC)
    pmax_ECMWF = MAXVAL(prob_forecast_ECMWF)
    
    ! ---- only verify if there is good data for all three models.
    
    IF (pmax_NCEP .ge. 0.0 .and. pmax_CMC .ge. 0.0 .and. pmax_ECMWF .ge. 0.0) THEN

	    DO itype = 1, 4
		    contab = 0.
		    bs = 0.
            bs_daily(idate) = 0.
		    IF (itype .eq. 1) THEN ! MME
			    prob_forecast(:,:) = prob_forecast_NCEP
		    ELSE IF (itype .eq. 2) THEN
                prob_forecast(:,:) = prob_forecast_CMC
		    ELSE IF (itype .eq. 3) THEN
                prob_forecast(:,:) = prob_forecast_ECMWF     
    		ELSE IF (itype .eq. 4) THEN
                prob_forecast(:,:) = 0.25*prob_forecast_NCEP + &
                    0.25*prob_forecast_CMC + 0.5*prob_forecast_ECMWF  
		    ENDIF  
	   
	   	    ! ---- now let's verify probability forecast

		    DO i = 1, nxa
          	    DO j = 1, nya
             	    IF (conusmask(i,j) .gt. 0 .and. prob_forecast(i,j) .GE. 0.0 .and. &
             	    prob_forecast(i,j) .le. 1.0 .and. apcp_anal_t(i,j,idate) .GE. 0.0) THEN
            	 	    cfac = cos(rlatsa(i,j)*pid180)  ! cos of latitude to acct for unequal grid box size
            	 	    pclimo = climo_prob(i,j)
            	 	    p      = prob_forecast(i,j)
            	 	    ipcat  = nint(p*20)
            	 	    v      = apcp_anal_t(i,j,idate)
				 
            	 	    IF (v .GE. rthresh) THEN
               	    	    contab(ipcat,2) = contab(ipcat,2) + cfac
                    	    bs = bs + cfac*(1.-p)**2
                            bs_daily(idate) = bs_daily(idate) + cfac*(1.-p)**2
                    	    IF (itype .eq. 1) THEN
                                bs_climo = bs_climo + cfac*(1.-pclimo)**2
                                bs_climo_daily(idate) = bs_climo_daily(idate) + cfac*(1.-pclimo)**2
                            ENDIF
            	 	    ELSE
                    	    contab(ipcat,1) = contab(ipcat,1) + cfac
                    	    bs = bs + cfac * p**2
                            bs_daily(idate) = bs_daily(idate) + cfac * p**2
                    	    IF (itype .eq. 1) THEN
                                bs_climo = bs_climo + cfac*pclimo**2
                                bs_climo_daily(idate) = bs_climo_daily(idate) + cfac * pclimo**2
                            ENDIF
                 	    ENDIF
             	    ENDIF ! conusmask
      	  	    END DO   ! j = 1, nya
   	   	    END DO      ! i = 1, nxa
	   
	   	    IF (itype .eq. 1) THEN
		   	    contab_NCEP = contab_NCEP + contab
			    contab_NCEP_daily(idate,:,:) = contab(:,:)
		   	    bs_NCEP = bs_NCEP + bs
                bs_NCEP_daily(idate) = bs_daily(idate)
	   	    ELSE IF (itype .eq. 2) THEN
		   	    contab_CMC = contab_CMC + contab
			    contab_CMC_daily(idate,:,:) = contab(:,:)
		   	    bs_CMC = bs_CMC + bs
                bs_CMC_daily(idate) = bs_daily(idate)
	        ELSE IF (itype .eq. 3) THEN
		   	    contab_ECMWF = contab_ECMWF + contab
			    contab_ECMWF_daily(idate,:,:) = contab(:,:)
		   	    bs_ECMWF = bs_ECMWF + bs
                bs_ECMWF_daily(idate) = bs_daily(idate)
    	    ELSE IF (itype .eq. 4) THEN
    		   	contab_MME = contab_MME + contab
    			contab_MME_daily(idate,:,:) = contab(:,:)
    		   	bs_MME = bs_MME + bs
                bs_MME_daily(idate) = bs_daily(idate)
	        ENDIF 
	   
   	    END DO  ! itype
        
    ELSE
        bs_climo_daily(idate) = 0.
        bs_NCEP_daily(idate) = 0.
        bs_CMC_daily(idate) = 0.
        bs_ECMWF_daily(idate) = 0.
        bs_MME_daily(idate) = 0. 
    ENDIF !(iex)
  
END DO ! idate

DO itype = 1, 4
	
	IF (itype .eq. 1) THEN
		PRINT *,'NCEP'
		contab = contab_NCEP
		contab_daily = contab_NCEP_daily
		bs = bs_NCEP
	ELSE IF (itype .eq. 2) THEN
		PRINT *, 'CMC'
		contab = contab_CMC
		contab_daily = contab_CMC_daily
		bs = bs_CMC
	ELSE IF (itype .eq. 3) THEN
		print *, 'ECMWF'
		contab = contab_ECMWF
		contab_daily = contab_ECMWF_daily
		bs = bs_ECMWF	
    ELSE IF (itype .eq. 4) THEN
    	print *, 'MME'
    	contab = contab_MME
    	contab_daily = contab_MME_daily
    	bs = bs_MME	
    ENDIF
	
	! ---- with tallied contingency tables, now set mean reliability and frequency of use for mean

	!print *,'bs, bs_climo = ',bs, bs_climo
	ctot = SUM(contab)
	bss = 1. - bs / bs_climo
	relia(:) = -99.9999
	PRINT *,'  bss = ',bss
	PRINT *,'  p   reliability   freq of usage'
	DO icat = 0,20
		frequse(icat) = (contab(icat,2)  + contab(icat,1)) / ctot
		IF ((contab(icat,1) + contab(icat,2)) .gt. 0) THEN
			relia(icat) = contab(icat,2) / (contab(icat,1) + contab(icat,2))
			PRINT 203,float(icat*5),relia(icat)*100.,frequse(icat)
		ELSE
			PRINT 203,float(icat*5),relia(icat),frequse(icat)
		ENDIF
  	    203 format(f5.1,3x,2(f8.3,3x))
	END DO  !icat
	
	! ---- perform a resampling to generate the confidence intervals for reliability
	!      sampling the days with replacement
	
	idum = -12345
	relia_resa(:,:) = -99.9999
	DO iresa = 1, nresa
        PRINT *,'iresa = ',iresa
		contab_sum = 0.
		DO idate = 1, ndates
			2345 cran = ran3(idum)
			idate2 = MIN(1+NINT(cran*REAL(ndates)),ndates)
            PRINT *,'cran,idate2 = ',cran,idate2
            IF (bs_MME_daily(idate2) .gt. 0.0) THEN
			    contab_sum(:,:) = contab_sum(:,:) + contab_daily(idate2,:,:)
            ELSE
                GOTO 2345
            ENDIF
		END DO
		DO icat = 0,20
			IF ((contab_sum(icat,1) + contab_sum(icat,2)) .gt. 0) THEN
				relia_resa(iresa,icat) = contab_sum(icat,2) / (contab_sum(icat,1) + contab_sum(icat,2))
			ENDIF
		END DO ! icat
	END DO  !iresa
	
    ! ---- now find the 5th and 95th percentiles of the resampled distribution, accounting
	!      for the possibility that at high probabilities there may be no samples in some 
	!      cases.

	DO i = 0, nclasses-1
		rsamps(:) = relia_resa(:,i)
		CALL sort(nresa, rsamps)
		ibegin = 1
		DO iresa = 1, nresa
			IF (rsamps(iresa) .gt. -99.) THEN
				ibegin = iresa
				goto 3212
			ENDIF
		END DO ! iresa
3212	isamp = MIN(ibegin + NINT(0.05* (nresa-ibegin+1) ),nresa)
		relia_05(i) = rsamps(isamp)
		isamp = MIN(ibegin + NINT(0.95* (nresa-ibegin+1) ),nresa)
		relia_95(i) = rsamps(isamp)
	END DO ! i 
	
	! ---- copy to output arrays
		
	IF (itype .eq. 1) THEN
		relia_NCEP = relia
		relia_NCEP_05 = relia_05
		relia_NCEP_95 = relia_95
		frequse_NCEP = frequse
		bss_NCEP = bss
	ELSE IF (itype .eq. 2) THEN
		relia_CMC = relia
		relia_CMC_05 = relia_05
		relia_CMC_95 = relia_95
		frequse_CMC = frequse
		bss_CMC = bss
	ELSE IF (itype .eq. 3) THEN
		relia_ECMWF = relia
		relia_ECMWF_05 = relia_05
		relia_ECMWF_95 = relia_95
		frequse_ECMWF = frequse
		bss_ECMWF = bss
    ELSE IF (itype .eq. 4) THEN
    	relia_MME = relia
    	relia_MME_05 = relia_05
    	relia_MME_95 = relia_95
    	frequse_MME = frequse
    	bss_MME = bss
	ENDIF
	
END DO  ! itype

! ---- get vectors of daily brier skill scores

PRINT *,'idate   BSS (NCEP   CMC   ECMWF    MME) '
DO idate = 1, ndates
    IF (bs_climo_daily(idate) .gt. 0.0) THEN
        bss_NCEP_daily(idate) = 1. - bs_NCEP_daily(idate)/bs_climo_daily(idate)
        bss_CMC_daily(idate) = 1. - bs_CMC_daily(idate)/bs_climo_daily(idate)
        bss_ECMWF_daily(idate) = 1. - bs_ECMWF_daily(idate)/bs_climo_daily(idate)
        bss_MME_daily(idate) = 1. - bs_MME_daily(idate)/bs_climo_daily(idate)
    ELSE
        bss_NCEP_daily(idate) = -99.99
        bss_CMC_daily(idate) = -99.99
        bss_ECMWF_daily(idate) = -99.99
        bss_MME_daily(idate) = -99.99
    ENDIF
    PRINT 302, idate, bss_NCEP_daily(idate), bss_CMC_daily(idate), bss_ECMWF_daily(idate), &
        bss_MME_daily(idate)
    302 FORMAT(i3,2x,4(f9.2,1x))
END DO

RETURN
END SUBROUTINE verify_relia_bss_mme

