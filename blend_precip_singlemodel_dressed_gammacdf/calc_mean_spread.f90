SUBROUTINE calc_mean_spread(n25, nxa, nya, nmembers, &
    ensemble_ccpa_x25, conusmask, ensmean, stddev, POP)
    
USE netcdf

INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
REAL, INTENT(IN), DIMENSION(n25, nxa, nya, nmembers) :: ensemble_ccpa_x25
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: ensmean, stddev
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: POP

REAL*8 :: sumxi, sumxi2, sn

ensmean(:,:) = 0.0
stddev(:,:) = 0.0
POP(:,:) = 0.0


!PRINT *,'n25, nxa, nya, nmembers = ', n25, nxa, nya, nmembers
!PRINT *,'ensemble_ccpa_x25(:,nxa/2, nya/2, :) = ', ensemble_ccpa_x25(:,nxa/2, nya/2, :) 
!PRINT *,'conusmask(1:nxa:10,nya/2) = ', conusmask(1:nxa:10,nya/2) 

DO jya = 1, nya
    !PRINT *,'jya = ',jya
    DO ixa = 1, nxa
        ntot = 0
        npos = 0
        sumxi = 0.0
        sumxi2 = 0.0
        sn = 0.0
        DO i25 = 1, n25
            DO imem = 1, nmembers
                e = ensemble_ccpa_x25(i25, ixa, jya, imem)
                sumxi = sumxi + e
                sumxi2 = sumxi2 + e**2
                sn = sn + 1.0
                IF (e .gt. 0.254) npos = npos+1
                ntot = ntot + 1
            END DO
        END DO
        ensmean(ixa,jya) = sumxi / sn
        stddev(ixa,jya) = (sumxi2 - sn*ensmean(ixa,jya)**2)/(sn-1.0)
        stddev(ixa,jya) = SQRT(stddev(ixa,jya))
        POP(ixa,jya) = REAL(npos) / REAL(ntot)
    END DO
END DO

RETURN 
END SUBROUTINE calc_mean_spread
        