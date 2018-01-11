SUBROUTINE read_terr_ht_and_facets_2p5 (nx_2p5, ny_2p5, &
    terrain_height, ifacet_short, ifacet_medium, ifacet_long,  &
    conusmask_2p5, lons_2p5, lats_2p5)
!
! purpose: read in previously calculated terrain height and facets..
! 
USE netcdf

INTEGER, INTENT(IN) :: nx_2p5, ny_2p5
INTEGER, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: ifacet_short ! terrain facet for small-scale smoothing
INTEGER, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: ifacet_medium ! terrain facet for med-scale smoothing
INTEGER, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: ifacet_long ! terrain facet for large-scale smoothing
REAL, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: terrain_height
REAL, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: lons_2p5
REAL, INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: lats_2p5
INTEGER(2), INTENT(OUT), DIMENSION(nx_2p5,ny_2p5) :: conusmask_2p5

INTEGER(2), DIMENSION(nx_2p5,ny_2p5) :: iwork
CHARACTER(len=120) :: infile
CHARACTER(len=20) :: cfield

! --- read in the terrain height, lon, lat from netCDF file 
!     prepared by Eric Engle.

infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/blend.precip_const.2p5.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

cfield = 'terrain'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,terrain_height,&
    start=(/1,1/),count=(/nx_2p5,ny_2p5/)))
    
cfield = 'longitude'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lons_2p5,&
    start=(/1,1/),count=(/nx_2p5,ny_2p5/)))
lons_2p5 = lons_2p5 - 360.    
    
cfield = 'latitude'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lats_2p5,&
    start=(/1,1/),count=(/nx_2p5,ny_2p5/)))        
    
CALL check(nf90_close(netid))

! --- read in the revised conus mask prepared by Tom Hamill 
!     (see revise_mask_ndfd2p5.py).

infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/blend.precip_const.2p5.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))
cfield = 'conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,conusmask_2p5,&
    start=(/1,1/),count=(/nx_2p5,ny_2p5/)))
CALL check(nf90_close(netid))

! ---- read in the terrain facets after the terrain heights
!      have been subjected to three different levels of smoothing

DO isize = 1, 3  ! short, medium, long smooth of terrain

    IF (isize .eq. 1) THEN
        infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/' // &
            'NDFD2p5_grid_conus_terrain_facets_smooth1.nc'
    ELSE IF (isize .eq. 2) THEN
        infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/' // &
            'NDFD2p5_grid_conus_terrain_facets_smooth10.nc'
    ELSE
        infile = '/data/thamill/Rf2_tests/ccpa_v1/ndfd2.5/' // &
            'NDFD2p5_grid_conus_terrain_facets_smooth20.nc'
    ENDIF

    PRINT *,'netid, reading from ',netid, TRIM(infile)
    CALL check (nf90_open(infile,NF90_NOWRITE,netid))
    cfield = 'terrain_facets'
    CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
    CALL check (nf90_get_var(netid,ivar,iwork,start=(/1,1/),&
        count=(/nx_2p5,ny_2p5/)))
    CALL check(nf90_close(netid))

    IF (isize .eq. 1) THEN
        ifacet_short = iwork
    ELSE IF (isize .eq. 2) THEN
        ifacet_medium = iwork
    ELSE
        ifacet_long = iwork
    ENDIF

END DO


RETURN
END subroutine read_terr_ht_and_facets_2p5
