SUBROUTINE read_facets (nxa, nya, ifacet_short, ifacet_medium, ifacet_long)
!
! purpose: use terrain facets created by terrain_facets_ccpa_extended.py
!    for purposes of classifying the direction the terrain faces.
!    Adapts algorithm presented in Chris Daly articles
! 
USE netcdf

INTEGER, INTENT(IN) :: nxa, nya
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: ifacet_short ! terrain facet for small-scale smoothing
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: ifacet_medium ! terrain facet for small-scale smoothing
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: ifacet_long ! terrain facet for small-scale smoothing

INTEGER*2, DIMENSION(nxa,nya) :: iwork
CHARACTER*80 infile
CHARACTER*20 cfield

DO isize = 1, 3  ! short, medium, long smooth of terrain

   IF (isize .eq. 1) THEN
      infile = 'CCPA_grid_conus_terrain_facets_smooth0.5.nc'
   ELSE IF (isize .eq. 2) THEN
      infile = 'CCPA_grid_conus_terrain_facets_smooth3.0.nc'
   ELSE
      infile = 'CCPA_grid_conus_terrain_facets_smooth10.0.nc'
   ENDIF

   ! ---- open the file
   PRINT *,'netid, reading from ',netid, TRIM(infile)
   CALL check (nf90_open(infile,NF90_NOWRITE,netid))

   cfield = 'terrain_facets'
   CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))

   CALL check (nf90_get_var(netid,ivar,iwork,start=(/1,1/),count=(/nxa,nya/)))

   IF (isize .eq. 1) THEN
      ifacet_short = iwork
   ELSE IF (isize .eq. 2) THEN
      ifacet_medium = iwork
   ELSE
      ifacet_long = iwork
   ENDIF
   
   print *, 'max/min iwork = ',maxval(iwork), minval(iwork)

   ! --- close file.

   CALL check(nf90_close(netid))

END DO


RETURN
END subroutine read_facets
