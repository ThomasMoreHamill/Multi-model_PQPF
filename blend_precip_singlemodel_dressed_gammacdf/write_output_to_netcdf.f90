SUBROUTINE write_output_to_netcdf(outfile, nxa, nya, nens, nthreshes, &
    pthreshes_in, rlonsa_in, rlatsa_in, climo_prob_in, &
    conusmask_in, prob_forecast_raw_in, prob_forecast_qmapped_in, &
    prob_forecast_in, ensemble_mean_in)
    
USE netcdf

CHARACTER*(*), INTENT(IN) :: outfile
INTEGER, INTENT(IN) :: nxa, nya, nens, nthreshes
REAL, INTENT(IN), DIMENSION(nthreshes) :: pthreshes_in
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlonsa_in, rlatsa_in, ensemble_mean_in
INTEGER*2, INTENT(IN),  DIMENSION(nxa,nya) :: conusmask_in
REAL, INTENT(IN), DIMENSION(nxa,nya,nthreshes) :: prob_forecast_raw_in, &
    prob_forecast_qmapped_in, prob_forecast_in, climo_prob_in
    
INTEGER :: dimid_3d(3), dimid_2d(2)
    
PRINT *,'writing to ',TRIM(outfile)

CALL check( nf90_create(TRIM(outfile), NF90_CLOBBER, ncid) )

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
CALL check( nf90_def_var(ncid, "prob_forecast_raw", &
    NF90_FLOAT, dimid_3d, npfraw_varid) )
CALL check( nf90_def_var(ncid, "prob_forecast_qmapped", &
    NF90_FLOAT, dimid_3d, npfqmap_varid) )      
CALL check( nf90_def_var(ncid, "prob_forecast", &
    NF90_FLOAT, dimid_3d, npf_varid) ) 
CALL check( nf90_def_var(ncid, "climo_prob", &
    NF90_FLOAT, dimid_3d, nclimo_varid) ) 
CALL check( nf90_def_var(ncid, "ensemble_mean", &
    NF90_FLOAT, dimid_2d, nemean_varid) )    

! --- End define mode. This tells netCDF we are done defining metadata.

CALL check( nf90_enddef(ncid) )

! ---- write the data.

PRINT *,'writing data'
CALL check( nf90_put_var(ncid, nthreshes_varid, pthreshes_in))
CALL check( nf90_put_var(ncid, nconusmask_varid, conusmask_in))
CALL check( nf90_put_var(ncid, nrlonsa_varid, rlonsa_in))
CALL check( nf90_put_var(ncid, nrlatsa_varid, rlatsa_in))
CALL check( nf90_put_var(ncid, npfraw_varid, prob_forecast_raw_in))
CALL check( nf90_put_var(ncid, npfqmap_varid, prob_forecast_qmapped_in))
CALL check( nf90_put_var(ncid, npf_varid, prob_forecast_in))
CALL check( nf90_put_var(ncid, nclimo_varid, climo_prob_in))
CALL check( nf90_put_var(ncid, nemean_varid, ensemble_mean_in))

! ---- Close the file. This frees up any internal netCDF resources
!      associated with the file, and flushes any buffers.

CALL check( nf90_close(ncid) )
PRINT *, "*** SUCCESS writing netcdf file "

RETURN
END SUBROUTINE write_output_to_netcdf