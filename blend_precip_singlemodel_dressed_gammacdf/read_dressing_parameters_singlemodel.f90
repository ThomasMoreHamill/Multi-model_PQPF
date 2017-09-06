SUBROUTINE read_dressing_parameters_singlemodel (npvals, npcatvals, &
    nlo_int_hi_vals, nclim_p1_vals, nclim_vals, &
    infile_gamma_parameters, precip_values, gamma_shapes, &
    gamma_scales, fraction_zeros, gamma_shape_fclimpop, &
    gamma_scale_fclimpop, fraction_zeros_fclimpop, &
    climo_pop_thresholds)

USE netcdf
    
INTEGER, INTENT(IN) :: npvals, npcatvals, nlo_int_hi_vals, &
    nclim_p1_vals, nclim_vals
CHARACTER*(*), INTENT(IN) :: infile_gamma_parameters
    
REAL, INTENT(OUT), DIMENSION(npvals) :: precip_values
REAL, INTENT(OUT), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: gamma_shapes
REAL, INTENT(OUT), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: gamma_scales
REAL, INTENT(OUT), DIMENSION(npvals, npcatvals, nlo_int_hi_vals) :: fraction_zeros
REAL, INTENT(OUT), DIMENSION(nclim_p1_vals) :: gamma_shape_fclimpop
REAL, INTENT(OUT), DIMENSION(nclim_p1_vals) :: gamma_scale_fclimpop
REAL, INTENT(OUT), DIMENSION(nclim_p1_vals) :: fraction_zeros_fclimpop
REAL, INTENT(OUT), DIMENSION(nclim_vals) :: climo_pop_thresholds
    
PRINT *,TRIM(infile_gamma_parameters)

CALL check(nf90_open(TRIM(infile_gamma_parameters), NF90_NOWRITE, netid))

CALL check(nf90_inq_varid(netid, "precip_values", ivar))
CALL check(nf90_get_var(netid, ivar, precip_values, &
    start=(/1/), count=(/npvals/)))
        
CALL check(nf90_inq_varid(netid, "fraction_zeros", ivar))
CALL check(nf90_get_var(netid, ivar, fraction_zeros, &
    start=(/1,1,1/), count=(/npvals, npcatvals, nlo_int_hi_vals/)))
            
CALL check(nf90_inq_varid(netid, "gamma_shapes", ivar))
CALL check(nf90_get_var(netid, ivar, gamma_shapes, &
    start=(/1,1,1/), count=(/npvals, npcatvals, nlo_int_hi_vals/)))
    
CALL check(nf90_inq_varid(netid, "gamma_scales", ivar))
CALL check(nf90_get_var(netid, ivar, gamma_scales, &
    start=(/1,1,1/), count=(/npvals, npcatvals, nlo_int_hi_vals/)))
        
CALL check(nf90_inq_varid(netid, "gamma_shape_fclimpop", ivar))
CALL check(nf90_get_var(netid, ivar, gamma_shape_fclimpop, &
    start=(/1/), count=(/nclim_p1_vals/)))

CALL check(nf90_inq_varid(netid, "gamma_scale_fclimpop", ivar))
CALL check(nf90_get_var(netid, ivar, gamma_scale_fclimpop, &
    start=(/1/), count=(/nclim_p1_vals/)))
        
CALL check(nf90_inq_varid(netid, "fraction_zeros_fclimpop", ivar))
CALL check(nf90_get_var(netid, ivar, fraction_zeros_fclimpop, &
    start=(/1/), count=(/nclim_p1_vals/)))
    
CALL check(nf90_inq_varid(netid, "climo_pop_thresholds", ivar))
CALL check(nf90_get_var(netid, ivar, climo_pop_thresholds, &
    start=(/1/), count=(/nclim_vals/)))

CALL check(nf90_close(netid))
    
RETURN
END SUBROUTINE read_dressing_parameters_singlemodel
