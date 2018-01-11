! FILE: get_quantiles_linear.f90
!  f2py -c -m get_quantiles_linear get_quantiles_linear.f90

! precip_qret_a, istat = get_quantiles_linear(nthresh,npct,nya,nxa,pctv,thresh,CDFa)

! -------

SUBROUTINE get_quantiles_linear(nthresh, npct, nj, ni, pctv, thresh, apcp_CDF, apcp_quantiles, istat)

INTEGER, INTENT(IN) :: nthresh, npct, nj, ni
REAL, INTENT(IN), DIMENSION(npct) :: pctv   ! 0.01 to 0.99 by 0.01, the prob values assoc'd w quantiles
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh ! the log-scaling for thresholds to evaluate CDF at
REAL, INTENT(IN), DIMENSION(nthresh, nj, ni) :: apcp_CDF
REAL, INTENT(OUT), DIMENSION(npct, nj, ni) :: apcp_quantiles
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, npct, nj, ni, pctv, thresh, apcp_CDF
!f2py intent(out) apcp_quantiles, istat
!f2py depend(npct) pctv
!f2py depend(nthresh) thresh
!f2py depend(nthresh,nj,ni) apcp_CDF
!f2py depend(npct,nj,ni) apcp_quantiles

REAL, DIMENSION(nthresh) :: apcp_CDF_1d

DO jy = 1, nj
   DO ix = 1, ni
      apcp_CDF_1d(:) = apcp_CDF(:,jy,ix)
      DO iquant = 1, npct
         t = pctv(iquant)
         DO ithr = 1, nthresh-1
            IF (apcp_CDF_1d(ithr) .le. t .and. apcp_CDF_1d(ithr+1) .gt. t) THEN
               ! ---- do linear interpolation
               fac = (t-apcp_CDF_1d(ithr)) / (apcp_CDF_1d(ithr+1) - apcp_CDF_1d(ithr))
               value = (1.-fac)*thresh(ithr) + fac*thresh(ithr+1)
               apcp_quantiles(iquant,jy,ix) = value
            ENDIF
         END DO
      END DO
   END DO
END DO
istat = 1
RETURN
END SUBROUTINE get_quantiles_linear

 


