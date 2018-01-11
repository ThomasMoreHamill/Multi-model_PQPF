! FILE: get_cdf_cmc.f90
!  f2py -c -m get_cdf_cmc get_cdf_cmc.f90
! ================================================================================

SUBROUTINE get_cdf_cmc(nthresh, nmembers, nya, nxa, nsuppmax, yoffset, xoffset, &
    thresh, apcp_fcst_ens, xlocations, ylocations, nsupplemental, conusmask, &
    CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nmembers, nya, nxa, nsuppmax, yoffset, xoffset
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nmembers, nya,nxa) :: apcp_fcst_ens
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: xlocations
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: ylocations
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: nsupplemental, conusmask

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nmembers,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nmembers, nya, nxa, nsuppmax, thresh, yoffset, xoffset
!f2py intent(in) apcp_fcst_ens, xlocations, ylocations, nsupplemental, conusmask
!f2py depend(nthresh) thresh
!f2py depend(nmembers,nya,nxa) apcp_fcst_ens
!f2py depend(nsuppmax,nya,nxa) xlocations, ylocations
!f2py depend(nya,nxa) nsupplemental, conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

PRINT *,'get_cdf_cmc: nthresh, nmembers, nya, nxa, nsuppmax = ',&
    nthresh, nmembers, nya, nxa, nsuppmax
CDFwork = 0.
icount(:,:) = 0
rminv = minval(apcp_fcst_ens)
rmaxv = maxval(apcp_fcst_ens)
PRINT *,'rminv, rmaxv = ', rminv, rmaxv

DO jya = 1, nya
    !PRINT *,'jya = ',jya,' of ',nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. rminv .ge. -98.) THEN
            icount(jya,ixa) = icount(jya,ixa) + nsupplemental(jya,ixa)
            DO imem = 1, nmembers
                DO isupp = 1, nsupplemental(jya,ixa)
                    jy2 = ylocations(isupp,jya,ixa) - yoffset
                    ix2 = xlocations(isupp,jya,ixa) - xoffset
                    DO ithr = 1, nthresh
                        IF (apcp_fcst_ens(imem,jy2,ix2) .LE. thresh(ithr)) &
                            CDFwork(ithr,imem,jya,ixa) = CDFwork(ithr,imem,jya,ixa) + 1.0
                    END DO ! ithr
                END DO ! isupp
            END DO ! imem
        END IF ! conusmask, good fcst data
    END DO ! ixa
END DO ! jya

istat = 1
IF (rminv .lt. -98) istat = -1
RETURN
END SUBROUTINE get_cdf_cmc



