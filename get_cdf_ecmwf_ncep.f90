! FILE: get_cdf_precip_fcst_ecmwf_ncep.f90
!  f2py -c -m get_cdf_ecmwf_ncep get_cdf_ecmwf_ncep.f90
! ================================================================================

SUBROUTINE get_cdf_ecmwf_ncep(nthresh, nmembers, nya, nxa, nsuppmax, yoffset, xoffset, &
    nmembers2process, thresh, apcp_fcst_ens, xlocations, ylocations, nsupplemental, &
    conusmask, CDFwork, icount, istat)
    
    
INTEGER, INTENT(IN) :: nthresh, nmembers, nya, nxa, nsuppmax, yoffset, xoffset
INTEGER, INTENT(IN) :: nmembers2process
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: xlocations
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: ylocations
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: nsupplemental, conusmask
REAL, INTENT(IN), DIMENSION(nmembers, nya,nxa) :: apcp_fcst_ens

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nmembers, nya, nxa, nsuppmax, yoffset, xoffset, thresh
!f2py intent(in) nmembers2process
!f2py intent(in) apcp_fcst_ens, conusmask, xlocations, ylocations, nsupplemental
!f2py depend(nthresh) thresh
!f2py depend(nmembers,nya,nxa) apcp_fcst_ens
!f2py depend(nya,nxa) conusmask, nsupplemental
!f2py depend(nsuppmax,nya,nxa) xlocations, ylocations
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

!PRINT *,'get_cdf_ecmwf_ncep:  nthresh, nmembers, nya, nxa, nsuppmax = ', &
!    nthresh, nmembers, nya, nxa, nsuppmax
CDFwork = 0.
icount(:,:) = 0
rminv = minval(apcp_fcst_ens)
rmaxv = maxval(apcp_fcst_ens)
!print *,'min,max value of forecast = ',rminv, rmaxv
!print *,'conusmask(nya/2,1:nxa:5) = ',conusmask(nya/2,1:nxa:5)
!print *,'nmembers, nsupplemental(nya/2, nxa/2) = ',nmembers, nsupplemental(nya/2, nxa/2)
!print *,'thresh(1:16) =',thresh(1:16) 
!print *,'apcp_fcst_ens(:,nya/2,nxa/2) = ',apcp_fcst_ens(:,nya/2,nxa/2)
DO jya = 1, nya
    !PRINT *,'jya = ',jya,' of ',nya
    DO ixa = 1, nxa
        IF (conusmask(jya,ixa) .eq. 1 .and. rminv .ge. -98.) THEN
            icount(jya,ixa) = icount(jya,ixa) + nsupplemental(jya,ixa)*nmembers2process
            DO imem = 1, nmembers2process
                !IF (jya .eq. nya/2 .and. ixa .eq. nxa/2) &
                !    print *,'imem, icount(nya/2,nxa/2) = ',imem, icount(nya/2,nxa/2) 
                DO isupp = 1, nsupplemental(jya,ixa)
                    jy2 = ylocations(isupp,jya,ixa) - yoffset
                    ix2 = xlocations(isupp,jya,ixa) - xoffset
                    DO ithr = 1, nthresh
                        IF (apcp_fcst_ens(imem,jy2,ix2) .LE. thresh(ithr)) &
                            CDFwork(ithr,jya,ixa) = CDFwork(ithr,jya,ixa) + 1.0
                    END DO ! ithr
                END DO ! isupp
            END DO ! imem
        END IF ! conusmask, good fcst data
    END DO ! ixa
END DO ! jya

IF (rminv .ge. -98) THEN
    istat = 1
ELSE
    istat = -1
ENDIF

RETURN
END SUBROUTINE get_cdf_ecmwf_ncep



