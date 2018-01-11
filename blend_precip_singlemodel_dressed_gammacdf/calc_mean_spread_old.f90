SUBROUTINE calc_mean_spread(n25, nxa, nya, nmembers, &
    ensemble_ccpa_x25, conusmask, ensmean, stddev, &
    ensmean_pxform, stddev_pxform)
    
USE netcdf

INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
REAL, INTENT(IN), DIMENSION(n25, nxa, nya, nmembers) :: ensemble_ccpa_x25
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: ensmean, stddev
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: ensmean_pxform, stddev_pxform

REAL*8 :: sumxi, sumxi2, sn, sumxi_pxform, sumxi2_pxform

ensmean(:,:) = 0.0
stddev(:,:) = 0.0
DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(ixa, jya) .eq. 1) THEN
            sumxi = 0.0
            sumxi2 = 0.0
            !sumxi_pxform = 0.0
            !sumxi2_pxform = 0.0
            sn = 0.0
            DO i25 = 1, n25
                DO imem = 1, nmembers
                    e = ensemble_ccpa_x25(i25, ixa, jya, imem)
                    sumxi = sumxi + e
                    sumxi2 = sumxi2 + e**2
                    !e = ensemble_ccpa_x25(i25, ixa, jya, imem)**0.33
                    !sumxi_pxform = sumxi_pxform + e
                    !sumxi2_pxform = sumxi2_pxform + e**2
                    sn = sn + 1.0
                END DO
            END DO
            ensmean(ixa,jya) = sumxi / sn
            !ensmean_pxform(ixa,jya) = sumxi_pxform / sn
            ensmean_pxform(ixa,jya) = ensmean(ixa,jya)**0.33
            stddev(ixa,jya) = (sumxi2 - sn*ensmean(ixa,jya)**2)/(sn-1.0)
            stddev(ixa,jya) = SQRT(stddev(ixa,jya))
            !stddev_pxform(ixa,jya) = &
            !    (sumxi2_pxform - sn*ensmean_pxform(ixa,jya)**2)/(sn-1.0)
            !stddev_pxform(ixa,jya) = SQRT(stddev_pxform(ixa,jya))
            stddev_pxform(ixa,jya) = stddev(ixa,jya)**0.33
            
        END IF
    END DO
END DO

RETURN 
END SUBROUTINE calc_mean_spread
        