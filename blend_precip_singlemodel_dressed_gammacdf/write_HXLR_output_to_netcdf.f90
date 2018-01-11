SUBROUTINE write_HXLR_output_to_netcdf(outfile_HXLR, nxa, nya, nthreshes, &
    pthreshes_in, rlonsa_in, rlatsa_in, conusmask_in, prob_forecast_in)
    
! tailored to write only the output of heteroskedastic ensemble logistic
! regression forecasts    
    
USE netcdf

CHARACTER*(*), INTENT(IN) :: outfile_HXLR
INTEGER, INTENT(IN) :: nxa, nya,  nthreshes
REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes_in
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlonsa_in, rlatsa_in 
INTEGER*2, INTENT(IN),  DIMENSION(nxa,nya) :: conusmask_in
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: prob_forecast_in
    
INTEGER :: dimid_3d(3), dimid_2d(2)
    
PRINT *,'writing to ',TRIM(outfile_HXLR)

CALL check( nf90_create(TRIM(outfile_HXLR), NF90_CLOBBER, ncid) )

! ---- Define the array dimensions. NetCDF will hand back an ID for each.

PRINT *,'array dimensions'
CALL check( nf90_def_dim(ncid, "nthreshes", nthreshes, nthreshes_dimid) )
CALL check( nf90_def_dim(ncid, "nxa", nxa, nxa_dimid) )
CALL check( nf90_def_dim(ncid, "nya", nya, nya_dimid) )
dimid_2d =  (/ nxa_dimid, nya_dimid /)
dimid_3d =  (/ nxa_dimid, nya_dimid, nthreshes_dimid /)

! ---- Define the variables and associated IDs

PRINT *,'defining variables'
CALL check( nf90_def_var(ncid, "pthreshes", &
    NF90_FLOAT, nthreshes_dimid, nthreshes_varid) )
CALL check( nf90_def_var(ncid, "conusmask", &
    NF90_SHORT, dimid_2d, nconusmask_varid) )    
CALL check( nf90_def_var(ncid, "rlonsa", &
    NF90_FLOAT, dimid_2d, nrlonsa_varid) )
CALL check( nf90_def_var(ncid, "rlatsa", &
    NF90_FLOAT, dimid_2d, nrlatsa_varid) )       
CALL check( nf90_def_var(ncid, "prob_forecast", &
    NF90_FLOAT, dimid_3d, npf_varid) )   

! --- End define mode. This tells netCDF we are done defining metadata.

CALL check( nf90_enddef(ncid) )

! ---- write the data.

PRINT *,'writing data'
CALL check( nf90_put_var(ncid, nthreshes_varid, pthreshes_in))
CALL check( nf90_put_var(ncid, nconusmask_varid, conusmask_in))
CALL check( nf90_put_var(ncid, nrlonsa_varid, rlonsa_in))
CALL check( nf90_put_var(ncid, nrlatsa_varid, rlatsa_in))
CALL check( nf90_put_var(ncid, npf_varid, prob_forecast_in))

! ---- Close the file. This frees up any internal netCDF resources
!      associated with the file, and flushes any buffers.

CALL check( nf90_close(ncid) )
PRINT *, "*** SUCCESS writing netcdf file "

RETURN
END SUBROUTINE write_HXLR_output_to_netcdf