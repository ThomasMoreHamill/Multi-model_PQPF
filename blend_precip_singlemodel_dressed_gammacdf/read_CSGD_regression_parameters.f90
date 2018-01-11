SUBROUTINE read_CSGD_regression_parameters(nxa, nya, infile_CSGD_regress_params, &
    ncsgd_params, csgd_parameters, rho)
    
    
USE netcdf
CHARACTER*(*), INTENT(IN) :: infile_CSGD_regress_params
INTEGER, INTENT(IN) :: ncsgd_params
REAL*8, INTENT(OUT), DIMENSION(ncsgd_params) :: csgd_parameters
REAL*4, DIMENSION(nxa,nya), INTENT(OUT) :: rho  !, shift

LOGICAL exist
CHARACTER*28 cfield
REAL, DIMENSION(ncsgd_params) :: csgd_parameters_rx4 

PRINT *,'in subroutine read_CSGD_regression_parameters'
INQUIRE (file=infile_CSGD_regress_params, exist=exist)
PRINT *, 'exist = ', exist
IF (exist) THEN

    CALL check (nf90_open(infile_CSGD_regress_params,NF90_NOWRITE,netid))

    cfield='rho_0p5'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,rho,&
        start=(/1,1/),count=(/nxa,nya/)))
    
    !cfield='shift'
    !CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    !CALL check(nf90_get_var(netid,ivar,shift,&
    !    start=(/1,1/),count=(/nxa,nya/)))
    
    cfield='alpha'
    CALL check(nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check(nf90_get_var(netid,ivar,csgd_parameters_rx4,&
        start=(/1/),count=(/ncsgd_params/)))
    CALL check(nf90_close(netid))
    csgd_parameters(:) = csgd_parameters_rx4(:)
    
   
ELSE

    PRINT *,'Unable to find file = ',TRIM(infile_CSGD_regress_params)
    PRINT *,'Stopping in subroutine rread_CSGD_regression_parameters'
    STOP

ENDIF

!csgd_parameters(1) = 0.04705411
!csgd_parameters(2) = 0.23096547
!csgd_parameters(3) = 0.3379612
!csgd_parameters(4) = 1.24825577
!csgd_parameters(5) = 0.35822056
!csgd_parameters(6) = 0.78102319  
    
RETURN
END SUBROUTINE read_CSGD_regression_parameters