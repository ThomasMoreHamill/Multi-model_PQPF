SUBROUTINE load_precipquantiles_2p5 (cmonth, nx_2p5, ny_2p5, npct, &
     cfields, panal_quantiles, fraction_zero)
!
! purpose: read in the data set for this month of the precipitation quantiles, 
!    both forecast and analyzed, as well as the rank correlation between 
!    forecast and analysis.
! 
USE netcdf

CHARACTER (LEN=2) :: cmonth
INTEGER, INTENT(IN) :: nx_2p5, ny_2p5, npct
CHARACTER (LEN=20), DIMENSION(npct) :: cfields

REAL, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5,npct) :: panal_quantiles
REAL, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: fraction_zero
REAL, DIMENSION(nx_2p5,ny_2p5) :: nonexceedance

CHARACTER (LEN=80) infile
CHARACTER (LEN=20) cfield
              
! ---- open the file, read in data

infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/' // &
    'gamma_climatology_parameters_' // cmonth // '.nc'

PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

DO ithresh = 1, npct
    cfield = cfields(ithresh)
    CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check (nf90_get_var(netid,ivar,nonexceedance, &
        start=(/1,1/),count=(/nx_2p5,ny_2p5/)))
    panal_quantiles(:,:,ithresh) = nonexceedance(:,:)
END DO

cfield = 'fraction_zero'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,fraction_zero, &
    start=(/1,1/),count=(/nx_2p5,ny_2p5/)))

! --- close file.

CALL check(nf90_close(netid))

RETURN
END subroutine load_precipquantiles_2p5
