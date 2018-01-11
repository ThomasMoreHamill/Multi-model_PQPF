! FILE: get_quantiles_linear.f90
!  f2py -c -m get_quantiles_linear_cmc get_quantiles_linear_cmc.f90

! precip_qret_a, istat = get_quantiles_linear_CMC(nthresh,npct,nya,nxa,pctv,thresh,CDFa)

! -------

SUBROUTINE get_quantiles_linear_cmc(nthresh, npct, nmembers, nj, ni, pctv, &
    thresh, apcp_CDF, apcp_quantiles, istat)

INTEGER, INTENT(IN) :: nthresh, npct, nmembers, nj, ni
REAL, INTENT(IN), DIMENSION(npct) :: pctv   ! 0.01 to 0.99 by 0.01, the prob values assoc'd w quantiles
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh ! the log-scaling for thresholds to evaluate CDF at
REAL, INTENT(IN), DIMENSION(nthresh, nmembers, nj, ni) :: apcp_CDF
REAL, INTENT(OUT), DIMENSION(npct, nmembers, nj, ni) :: apcp_quantiles
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, npct, nj, ni, pctv, thresh, apcp_CDF
!f2py intent(out) apcp_quantiles, istat
!f2py depend(npct) pctv
!f2py depend(nthresh) thresh
!f2py depend(nthresh,nmembers,nj,ni) apcp_CDF
!f2py depend(npct,nmembers,nj,ni) apcp_quantiles

REAL, DIMENSION(nthresh) :: apcp_CDF_1d

DO jy = 1, nj
    DO ix = 1, ni
        DO imem = 1, nmembers
            apcp_CDF_1d(:) = apcp_CDF(:,imem,jy,ix)
            DO iquant = 1, npct
                t = pctv(iquant)
                DO ithr = 1, nthresh-1
                    IF (apcp_CDF_1d(ithr) .le. t .and. apcp_CDF_1d(ithr+1) .gt. t) THEN
                        ! ---- do linear interpolation
                        fac = (t-apcp_CDF_1d(ithr)) / (apcp_CDF_1d(ithr+1) - apcp_CDF_1d(ithr))
                        value = (1.-fac)*thresh(ithr) + fac*thresh(ithr+1)
                        apcp_quantiles(iquant,imem,jy,ix) = value
                    ENDIF ! apcp_CDF_1d(ithr) .le. t .and. apcp_CDF_1d(ithr+1) .gt. t
                END DO ! ithr = 1, nthresh-1
            END DO ! iquant = 1, npct
        END DO ! imem = 1, nmembers
    END DO ! ix = 1, ni
END DO ! jy = 1, nj
istat = 1

RETURN
END SUBROUTINE get_quantiles_linear_cmc

 


