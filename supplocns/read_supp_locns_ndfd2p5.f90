SUBROUTINE read_supp_locns_ndfd2p5(cmonth, nxa, nya, nsupp, &
    xlocation, ylocation, conusmask, lons, lats)
    
! compile with f2py -c -m read_supp_locns_ndfd2p5 read_supp_locns_ndfd2p5.f90 
    
CHARACTER*3, INTENT(IN) :: cmonth
INTEGER, INTENT(IN) :: nxa, nya, nsupp

INTEGER(2), INTENT(OUT), DIMENSION(nxa,nya) :: conusmask 
    ! 1/0 mask for whether in NDFD grid
INTEGER, INTENT(OUT), DIMENSION(nxa,nya,nsupp) :: xlocation
    ! for each forecast point,  a list of which other forecast points
INTEGER, INTENT(OUT), DIMENSION(nxa,nya,nsupp) :: ylocation
    ! have the closest climatologies and forecast-obs relationship
REAL, INTENT(OUT), DIMENSION(nxa,nya)   :: lons
REAL, INTENT(OUT), DIMENSION(nxa,nya)   :: lats

! f2py intent(in) cmonth, nxa, nya, nsupp
! f2py intent(out) conusmask, xlocation, ylocation, lons, lats
! f2py depend(nxa, nya) conusmask,  lons, lats
! f2py depend(nxa, nya, nsupp) xlocation, ylocation

CHARACTER (len=120) :: infile

! --- read data for this month

infile = '/Users/thamill/refcst2/test/supplemental_locations_NDFD2p5_'//cmonth//'.dat'
PRINT *, 'reading from ', TRIM(infile)
OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='unformatted')
READ (1) nx_2p5_in, ny_2p5_in, nsupp_in, nx_ccpa_in, ny_ccpa_in, nsupp_ccpa_in
PRINT *, 'nx_2p5_in, ny_2p5_in, nsupp_in, nx_ccpa_in, ny_ccpa_in, nsupp_ccpa_in'
PRINT *, nx_2p5_in, ny_2p5_in, nsupp_in, nx_ccpa_in, ny_ccpa_in, nsupp_ccpa_in
READ (1) xlocation
READ (1) ylocation
READ (1) conusmask
READ (1) lons
READ (1) lats
CLOSE (1)

RETURN
END SUBROUTINE read_supp_locns_ndfd2p5