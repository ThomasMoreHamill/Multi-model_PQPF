! FILE: get_cdf_precip_anal.f90
!  f2py -c -m get_cdf_precip_anal get_cdf_precip_anal.f90
! ======================


!CDFworka, icounta, istata = get_cdf_precip_anal(nthresh, nya, nxa, nsuppmax, \
!   thresh, apcp_anal, xlocations, ylocations, nsupplemental, conusmask_in)

SUBROUTINE get_cdf_precip_anal(nthresh, nsuppmax, nya, nxa, yoffset, xoffset, &
    thresh, apcp_anal, xlocations, ylocations, nsupplemental, conusmask, &
    CDFwork, icount, istat)
    
INTEGER, INTENT(IN) :: nthresh, nya, nxa, nsuppmax, yoffset, xoffset
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
REAL, INTENT(IN), DIMENSION(nya,nxa) :: apcp_anal
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: xlocations
INTEGER, INTENT(IN), DIMENSION(nsuppmax,nya,nxa) :: ylocations
INTEGER*2, INTENT(IN), DIMENSION(nya,nxa) :: nsupplemental, conusmask

INTEGER, INTENT(OUT), DIMENSION(nya,nxa) :: icount
REAL*8, INTENT(OUT), DIMENSION(nthresh,nya,nxa) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nsuppmax, nya, nxa, thresh, yoffset, xoffset
!f2py intent(in) apcp_anal, xlocations, ylocations, nsupplemental, conusmask
!f2py depend(nthresh) thresh
!f2py depend(nya,nxa) apcp_anal
!f2py depend(nsuppmax,nya,nxa) xlocations, ylocations
!f2py depend(nya,nxa) nsupplemental, conusmask
!f2py intent(out) CDFwork, icount, istat
!f2py depend(nthresh, nya, nxa) CDFwork
!f2py depend(nya,nxa) icount

!print *,'starting get_cdf_precip_anal'
CDFwork = 0.
icount(:,:) = 0
rmaxprecip = maxval(apcp_anal)

!print *,'apcp_anal(nya/2,:)= ', apcp_anal(nya/2,:)
!print *,'conusmask(nya/2,:) = ', conusmask(nya/2,:)
!print *,'xlocations(:,nya/2, nxa/2) = ',xlocations(:,nya/2, nxa/2)
!print *,'ylocations(:,nya/2, nxa/2) = ',ylocations(:,nya/2, nxa/2)
!print *,'nsupplemental(nya/2, nxa/2) = ', nsupplemental(nya/2, nxa/2)
!print *,'thresh(1:30) = ', thresh(1:30)

DO jya = 1, nya
    DO ixa = 1, nxa
        !DO jya = 2,2
!    DO ixa = 355,355
        IF (conusmask(jya,ixa) .eq. 1 .and. rmaxprecip .ge. 0.0) THEN
            icount(jya,ixa) = icount(jya,ixa) + nsupplemental(jya,ixa)
            DO isupp = 1, nsupplemental(jya,ixa)
                jy2 = ylocations(isupp,jya,ixa) - yoffset
                ix2 = xlocations(isupp,jya,ixa) - xoffset
                DO ithr = 1, nthresh
                    IF (apcp_anal(jy2,ix2) .LE. thresh(ithr)) &
                        CDFwork(ithr,jya,ixa) = CDFwork(ithr,jya,ixa) + 1.0
                END DO ! ithr
            END DO ! isupp
        END IF ! conusmask
    END DO ! ixa
END DO ! jya

IF (rmaxprecip .lt. 0.0) THEN
   istat = -1
ELSE
   istat = 1
ENDIF

!PRINT *,'CDFwork(1:15,nya/2,nxa/2) = ', CDFwork(1:15,nya/2,nxa/2)

RETURN
END SUBROUTINE get_cdf_precip_anal



