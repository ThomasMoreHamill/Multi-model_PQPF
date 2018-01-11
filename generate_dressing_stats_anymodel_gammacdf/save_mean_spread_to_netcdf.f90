SUBROUTINE save_mean_spread_to_netcdf(n25, nxa, nya, nmembers, &
    data_directory, cyyyymmddhh, cmodel, cleade, &
    ensemble_ccpa_x25, conusmask)
    
USE netcdf

INTEGER, INTENT(IN) :: n25, nxa, nya, nmembers
CHARACTER*(*), INTENT(IN) :: data_directory, cyyyymmddhh, cmodel, cleade
REAL, INTENT(IN), DIMENSION(n25, nxa, nya, nmembers) :: ensemble_ccpa_x25
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask

REAL*8 :: sumxi, sumxi2, sn, sumxi_pxform, sumxi2_pxform
REAL, DIMENSION(nxa, nya) :: ensmean, stddev, POP
REAL, DIMENSION(nxa, nya) :: ensmean_pxform, stddev_pxform
CHARACTER*256 outfile
INTEGER :: dimid_2d(2)

!PRINT *,'Subroutine save_mean_spread_to_netcdf'
!PRINT *,'maxval(ensemble_ccpa_x25*conusmask) = ', maxval(ensemble_ccpa_x25) 
!PRINT *,'sum conusmask = ', SUM(INT(conusmask))
!PRINT *,'nya, nxa, n25, nmembers = ', nya, nxa, n25, nmembers 
npoints = SUM(INT(conusmask))
ensmean(:,:) = 0.0
stddev(:,:) = 0.0
POP(:,:) = 0.0
DO jya = 1, nya
    DO ixa = 1, nxa
        nexceed = 0
        ntot = 0
        IF (conusmask(ixa, jya) .eq. 1) THEN
            sumxi = 0.0
            sumxi2 = 0.0
            sumxi_pxform = 0.0
            sumxi2_pxform = 0.0
            sn = 0.0
            
            DO i25 = 1, n25
                DO imem = 1, nmembers
                    e = ensemble_ccpa_x25(i25, ixa, jya, imem)
                    sumxi = sumxi + e
                    sumxi2 = sumxi2 + e**2
                    ntot = ntot+1
                    IF (e .gt. 0.254) nexceed = nexceed+1
                    e = ensemble_ccpa_x25(i25, ixa, jya, imem)**0.5
                    sumxi_pxform = sumxi_pxform + e
                    sumxi2_pxform = sumxi2_pxform + e**2
                    sn = sn + 1.0
                END DO
            END DO
            ensmean(ixa,jya) = sumxi / sn
            ensmean_pxform(ixa,jya) = sumxi_pxform / sn
            !ensmean_pxform(ixa,jya) = ensmean(ixa,jya)**0.33
            POP(ixa,jya) = REAL(nexceed) / REAL(ntot)
            !IF (jya .eq. nya/2) THEN
            !    PRINT *,'ixa, sumxi2, sn*ensmean(ixa,jya)**2, sn = ', ixa, sumxi2, sn*ensmean(ixa,jya)**2, sn 
            !ENDIF
            stddev(ixa,jya) = (sumxi2 - sn*ensmean(ixa,jya)**2)/(sn-1.0)
            stddev(ixa,jya) = SQRT(stddev(ixa,jya))
            !stddev_pxform(ixa,jya) = (sumxi2_pxform - sn*ensmean_pxform(ixa,jya)**2)/(sn-1.0)
            !stddev_pxform(ixa,jya) = SQRT(stddev_pxform(ixa,jya))
            stddev_pxform(ixa,jya) = stddev(ixa,jya)**0.33
            
        END IF
    END DO
END DO

outfile = TRIM(data_directory) // TRIM(cmodel) // '/' // &
    TRIM(cmodel) // '_' // TRIM(cyyyymmddhh) // '_mean_spread_' // TRIM(cleade) // 'h.nc'
PRINT *,'writing to ', TRIM(outfile)

! ---- Create the netCDF file.

CALL check( nf90_create(TRIM(outfile), NF90_CLOBBER, ncid) )

! ---- Define the array dimensions. NetCDF will hand back an ID for each.

CALL check( nf90_def_dim(ncid, "nxa", nxa, nxa_dimid) )
CALL check( nf90_def_dim(ncid, "nya", nya, nya_dimid) )
dimid_2d =  (/ nxa_dimid, nya_dimid /)

! ---- Define the variables and associated IDs

!PRINT *,'defining variables'
CALL check( nf90_def_var(ncid, "conusmask", &
    NF90_SHORT, dimid_2d, nconusmask_varid) )
CALL check( nf90_def_var(ncid, "ensemble_mean_qmapped_0p33", &
    NF90_FLOAT, dimid_2d, nmean_0p33_varid) )    
CALL check( nf90_def_var(ncid, "ensemble_stddev_qmapped_0p33", &
    NF90_FLOAT, dimid_2d, nstddev_0p33_varid) ) 
CALL check( nf90_def_var(ncid, "ensemble_mean_qmapped", &
    NF90_FLOAT, dimid_2d, nmean_varid) ) 
CALL check( nf90_def_var(ncid, "ensemble_mean_of_sqrt_qmapped", &
    NF90_FLOAT, dimid_2d, nmean_sqrt_varid) )    
       
CALL check( nf90_def_var(ncid, "ensemble_stddev_qmapped", &
    NF90_FLOAT, dimid_2d, nstddev_varid) )   
CALL check( nf90_def_var(ncid, "POP", &
    NF90_FLOAT, dimid_2d, npop_varid) )    
    
!PRINT *,'max ensemble_mean = ', &
!    maxval(ensmean)
!PRINT *,'max stddev = ', &
!    maxval(stddev)    
!PRINT *,'max ensemble_mean_qmapped_0p33 = ', &
!    maxval(ensmean_pxform)
!PRINT *,'max stddev_0p33 = ', &
!    maxval(stddev_pxform)

! --- End define mode. This tells netCDF we are done defining metadata.

CALL check( nf90_enddef(ncid) )

! ---- write the data.

CALL check( nf90_put_var(ncid, nconusmask_varid, conusmask))
CALL check( nf90_put_var(ncid, nmean_varid, ensmean))
CALL check( nf90_put_var(ncid, nstddev_varid, stddev))
CALL check( nf90_put_var(ncid, nmean_sqrt_varid, ensmean_pxform))
CALL check( nf90_put_var(ncid, npop_varid, POP))

! ---- Close the file. This frees up any internal netCDF resources
!      associated with the file, and flushes any buffers.

CALL check( nf90_close(ncid) )

RETURN 
END SUBROUTINE save_mean_spread_to_netcdf
        