SUBROUTINE determine_nearest_ccpa(nx_2p5, ny_2p5, nx_ccpa, ny_ccpa, &
    lons_2p5, lats_2p5, lons_ccpa, lats_ccpa, conusmask_ccpa, &
    inearest_ccpa, jnearest_ccpa, inearest_2p5, jnearest_2p5)
    
INTEGER, INTENT(IN) :: nx_2p5, ny_2p5, nx_ccpa, ny_ccpa
REAL, INTENT(IN), DIMENSION(nx_2p5, ny_2p5) :: lons_2p5, lats_2p5
REAL, INTENT(IN), DIMENSION(nx_ccpa, ny_ccpa) :: lons_ccpa, lats_ccpa
INTEGER(2), INTENT(IN), DIMENSION(nx_ccpa, ny_ccpa) :: conusmask_ccpa

INTEGER, INTENT(OUT), DIMENSION(nx_2p5, ny_2p5) :: inearest_ccpa, jnearest_ccpa
INTEGER, INTENT(OUT), DIMENSION(nx_ccpa, ny_ccpa) :: inearest_2p5, jnearest_2p5


REAL, DIMENSION(nx_ccpa) :: lons_ccpa_1d
REAL, DIMENSION(ny_ccpa) :: lats_ccpa_1d
REAL, DIMENSION(nx_ccpa) :: lon_differences
REAL, DIMENSION(ny_ccpa) :: lat_differences
REAL, DIMENSION(nx_ccpa, ny_ccpa) :: sumi, sumj, rnsamps

inearest_2p5(:,:) = -99
jnearest_2p5(:,:) = -99
inearest_ccpa(:,:) = -99
jnearest_ccpa(:,:) = -99

! ---- CCPA grid is regular lat-lon, so effectively 1D

lons_ccpa_1d(:) = lons_ccpa(:,1)
lats_ccpa_1d(:) = lats_ccpa(1,:)
sumi(:,:) = 0.0
sumj(:,:) = 0.0
rnsamps(:,:) = 0.0

! ---- loop thru and find the nearest CCPA point.  If this point
!      is greater than 0.09 (~ 1/8th degree / 2 * sqrt(2)), 
!      then set to missing.

DO jy = 1, ny_2p5
    DO ix = 1, nx_2p5        
        rlon = lons_2p5(ix,jy)
        rlat = lats_2p5(ix,jy)
        iccpa = NINT(8.*(rlon-lons_ccpa(1,1)) + 1)  ! 8 b/c 1/8 degree 
        jccpa = NINT(8.* (rlat-lats_ccpa(1,1)) + 1)
        IF (iccpa .ge. 1 .and. iccpa .le. nx_ccpa .and. &
        jccpa .ge. 1 .and. jccpa .le. ny_ccpa) THEN
            IF (conusmask_ccpa(iccpa,jccpa) .eq. 1) THEN
                inearest_ccpa(ix,jy) = iccpa
                jnearest_ccpa(ix,jy) = jccpa  
                sumi(iccpa,jccpa) = sumi(iccpa,jccpa) + REAL(ix)  
                sumj(iccpa,jccpa) = sumj(iccpa,jccpa) + REAL(jy)  
                rnsamps(iccpa,jccpa) = rnsamps(iccpa,jccpa) + 1 
            ELSE
                inearest_ccpa(ix,jy) = -99
                jnearest_ccpa(ix,jy) = -99
            ENDIF            
        ELSE
            inearest_ccpa(ix,jy) = -99
            jnearest_ccpa(ix,jy) = -99
        ENDIF
    END DO
END DO

! ---- now let's determine the nearest 2.5 km grid point to each
!      CCPA point by averaging

DO jy = 1, ny_ccpa
    DO ix = 1, nx_ccpa
        IF (rnsamps(ix,jy) .gt. 0.0) THEN
            inearest_2p5(ix,jy) = NINT(sumi(ix,jy) / rnsamps(ix,jy))
            jnearest_2p5(ix,jy) = NINT(sumj(ix,jy) / rnsamps(ix,jy))
        ENDIF
    END DO
END DO

RETURN
END SUBROUTINE determine_nearest_ccpa
