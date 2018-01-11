SUBROUTINE calculate_terrain_gradients_r8(nxa, nya, earth_radius_meters, &
    terrain_height, latsa, terrain_gradx, terrain_grady)
    
INTEGER, INTENT(IN) :: nxa, nya
REAL, INTENT(IN) :: earth_radius_meters
REAL, INTENT(IN), DIMENSION (nxa,nya) :: terrain_height, latsa
REAL*8, INTENT(OUT), DIMENSION (nxa,nya) :: terrain_gradx, terrain_grady

REAL*8 dy, rpi, dx, dtdy, dtdx

rpi = 3.1415926
DO jya = 1, nya
    dy = (2.*rpi*earth_radius_meters) / (360.*8) ! diameter /  1/8 degree
    dx = dy*cos(latsa(1,jya)*rpi/180.)

    DO ixa = 1, nxa
        
        IF (jya .eq. 1) THEN
            dtdy = (terrain_height(ixa,2)-terrain_height(ixa,1))/dy
        ELSE IF (jya .eq. nya) THEN
            dtdy = (terrain_height(ixa,nya)-terrain_height(ixa,nya-1))/dy    
        ELSE 
            dtdy = (terrain_height(ixa,jya+1)-terrain_height(ixa,jya-1))/(2.*dy) 
        ENDIF
        terrain_grady(ixa,jya) = -1.*dtdy
        
        
        IF (ixa .eq. 1) THEN
            dtdx = (terrain_height(2,jya)-terrain_height(1,jya))/dx
        ELSE IF (ixa .eq. nxa) THEN
            dtdx = (terrain_height(nxa,jya)-terrain_height(nxa-1,jya))/dx    
        ELSE 
            dtdx = (terrain_height(ixa+1,jya)-terrain_height(ixa-1,jya))/(2.*dx) 
        ENDIF
        
        terrain_gradx(ixa,jya) = -1.*dtdx
    END DO ! ixa
END DO ! jya
RETURN
END SUBROUTINE calculate_terrain_gradients_r8