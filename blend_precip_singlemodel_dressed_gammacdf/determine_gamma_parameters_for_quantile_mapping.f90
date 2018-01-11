SUBROUTINE determine_gamma_parameters_for_quantile_mapping(nxa, nya, &
    iyyyymmddhh, imonth, nens_qmap, data_directory, cmodel, cleade, &
    cmonths, conusmask, gamma_shape_qmap_forecast, gamma_scale_qmap_forecast, &
    fraction_zero_qmap_forecast, gamma_shape_qmap_analysis, &
    gamma_scale_qmap_analysis, fraction_zero_qmap_analysis)
    
USE netcdf

! --- this subroutine will read in the previous 60 days of synthesized
!     information on Gamma distributions as well as supplemental location 
!     information, calculating for each grid point an estimated fraction
!     zero (fraction of the samples with zero precip amount) and Gamma 
!     distribution parameters shape (alpha) and scale (beta).  This is done
!     for both forecast and analyzed data.
    
INTEGER, INTENT(IN) :: nxa, nya, iyyyymmddhh, nens_qmap
CHARACTER*(*), INTENT(IN) :: data_directory
CHARACTER*(*), INTENT(IN) :: cmodel
CHARACTER*(*), INTENT(IN) :: cleade

CHARACTER*3, DIMENSION(12) , INTENT(IN) :: cmonths
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask

REAL, INTENT(OUT), DIMENSION(nxa, nya, nens_qmap) :: gamma_shape_qmap_forecast
REAL, INTENT(OUT), DIMENSION(nxa, nya, nens_qmap) :: gamma_scale_qmap_forecast
REAL, INTENT(OUT), DIMENSION(nxa, nya, nens_qmap) :: fraction_zero_qmap_forecast
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: gamma_shape_qmap_analysis
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: gamma_scale_qmap_analysis
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: fraction_zero_qmap_analysis

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: xlocations_largedomain, ylocations_largedomain
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: xlocations, ylocations
INTEGER*2, ALLOCATABLE, DIMENSION(:,:) :: nsupplemental_largedomain, nsupplemental

! --- local variables

CHARACTER*20 cvar
CHARACTER*120 infile
CHARACTER*10 cyyyymmddhh
INTEGER xoffset, yoffset

REAL, ALLOCATABLE, DIMENSION(:,:) :: sumx_analysis, sumx_analysis_multiday
REAL, ALLOCATABLE, DIMENSION(:,:) :: sumlnx_analysis, sumlnx_analysis_multiday
INTEGER*2, ALLOCATABLE, DIMENSION(:,:) :: npositive_analysis, nzeros_analysis
INTEGER*4, ALLOCATABLE, DIMENSION(:,:) :: npositive_analysis_multiday, nzeros_analysis_multiday
INTEGER*4, ALLOCATABLE, DIMENSION(:,:) :: npositive_analysis_final, nzeros_analysis_final

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: sumx_forecast, sumx_forecast_multiday
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: sumlnx_forecast, sumlnx_forecast_multiday
INTEGER*2, ALLOCATABLE, DIMENSION(:,:,:) :: npositive_forecast, nzeros_forecast
INTEGER*4, ALLOCATABLE, DIMENSION(:,:,:) :: npositive_forecast_multiday, nzeros_forecast_multiday
INTEGER*4, ALLOCATABLE, DIMENSION(:,:,:) :: npositive_forecast_final, nzeros_forecast_final

LOGICAL exist

! ---- read in the supplemental locations dimension information 

infile = TRIM(data_directory) // 'supplemental_locations_eighth_degree_'// &
    cmonths(imonth) // '_v9.nc'
PRINT *,'reading from ', TRIM(infile)

CALL check(nf90_open(TRIM(infile),NF90_NOWRITE,ncid))

cvar = "xa"
CALL check(nf90_inq_dimid(ncid,cvar,idimid))
CALL check(nf90_inquire_dimension(ncid,idimid,cvar,nxa_large))

cvar = "ya"
CALL check(nf90_inq_dimid(ncid,cvar,idimid))
CALL check(nf90_inquire_dimension(ncid,idimid,cvar,nya_large))

cvar = "xsupp"
CALL check(nf90_inq_dimid(ncid,cvar,idimid))
CALL check(nf90_inquire_dimension(ncid,idimid,cvar,nsupp))
!print *,'after xsupp'
    
PRINT *,'nxa_large, nya_large, nsupp = ', nxa_large, nya_large, nsupp
ALLOCATE(xlocations(nxa,nya,nsupp), xlocations_largedomain(nxa_large, nya_large, nsupp))    
ALLOCATE(ylocations(nxa,nya,nsupp), ylocations_largedomain(nxa_large, nya_large, nsupp))     
ALLOCATE(nsupplemental(nxa,nya), nsupplemental_largedomain(nxa_large, nya_large)) 

! ---- read in the supplemental locations.  NOTE: they are defined on the 
!      enlarged domain that extends further west and south, and so we need 
!      to thin down the size of the domain and shift the values
!      of grid indices.

cvar = "xlocations"
CALL check(nf90_inq_varid(ncid,cvar,ivarid))
CALL check(nf90_get_var(ncid,ivarid,xlocations_largedomain))

cvar = "ylocations"
CALL check(nf90_inq_varid(ncid,cvar,ivarid))
CALL check(nf90_get_var(ncid,ivarid,ylocations_largedomain))

cvar = "nsupplemental"
CALL check(nf90_inq_varid(ncid,TRIM(cvar),ivarid))
CALL check(nf90_get_var(ncid,ivarid,nsupplemental_largedomain,&
    start=(/1,1/),count=(/nxa_large, nya_large/)))

CALL check(nf90_close(ncid))

! ---- do the shifting, excising the necessary part from the smaller grid and shifing indices.
!      (sorry, I moved processing back to a smaller grid covering CONUS)

xoffset = 21
yoffset = 22
isbegin = xoffset+1
isend = xoffset+nxa
jsbegin = yoffset+1
jsend = yoffset+nya
xlocations(:,:,:) = xlocations_largedomain(isbegin:isend,jsbegin:jsend,:)
ylocations(:,:,:) = ylocations_largedomain(isbegin:isend,jsbegin:jsend,:)
nsupplemental(:,:) = nsupplemental_largedomain(isbegin:isend,jsbegin:jsend)

DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(ixa,jya) .eq. 1) THEN 
            IF (nsupplemental(ixa,jya) .eq. 0) THEN
                PRINT *,'invalid nsupplemental value at ixa,jya = ',ixa,jya
                PRINT *,'stopping'
                stop
            endif
            DO isupp = 1, nsupplemental(ixa,jya)
                xlocations(ixa,jya,isupp) = xlocations(ixa,jya,isupp) - xoffset 
                ylocations(ixa,jya,isupp) = ylocations(ixa,jya,isupp) - yoffset
            ENDDO
        ENDIF
    END DO
END DO
!print *,'after shifting'
                
! ---- now, for the input date and the forecast lead time, we are going to read
!      in the prior 60 days of information necessary to build the gamma CDFs.

ALLOCATE(sumx_analysis(nxa,nya), sumx_analysis_multiday(nxa,nya))
ALLOCATE(sumlnx_analysis(nxa,nya), sumlnx_analysis_multiday(nxa,nya))
ALLOCATE(npositive_analysis(nxa,nya), nzeros_analysis(nxa,nya))
ALLOCATE(npositive_analysis_multiday(nxa,nya), nzeros_analysis_multiday(nxa,nya))
ALLOCATE(npositive_analysis_final(nxa,nya), nzeros_analysis_final(nxa,nya))

ALLOCATE(sumx_forecast(nxa,nya,nens_qmap), &
    sumx_forecast_multiday(nxa,nya,nens_qmap))
ALLOCATE(sumlnx_forecast(nxa,nya,nens_qmap), &
    sumlnx_forecast_multiday(nxa,nya,nens_qmap))
ALLOCATE(npositive_forecast(nxa,nya,nens_qmap), &
    nzeros_forecast(nxa,nya,nens_qmap))
ALLOCATE(npositive_forecast_multiday(nxa,nya,nens_qmap), &
    nzeros_forecast_multiday(nxa,nya,nens_qmap))
ALLOCATE(npositive_forecast_final(nxa,nya,nens_qmap), &
    nzeros_forecast_final(nxa,nya,nens_qmap))

sumx_analysis_multiday(:,:) = 0.0
sumlnx_analysis_multiday(:,:) = 0.0
npositive_analysis_multiday(:,:) = 0
nzeros_analysis_multiday(:,:) = 0

sumx_forecast_multiday(:,:,:) = 0.0
sumlnx_forecast_multiday(:,:,:) = 0.0
npositive_forecast_multiday(:,:,:) = 0
nzeros_forecast_multiday(:,:,:) = 0

READ (cleade,'(i3)') ileade
ndayshift = -1 - NINT(REAL(ileade)/24.)  ! the most recent day with valid forecast/analysis data

ktr = 1
DO iday = ndayshift-60, ndayshift-1

    sumx_analysis(:,:) = 0.0
    sumlnx_analysis(:,:) = 0.0
    npositive_analysis(:,:) = 0
    nzeros_analysis(:,:) = 0

    sumx_forecast(:,:,:) = 0.0
    sumlnx_forecast(:,:,:) = 0.0
    npositive_forecast(:,:,:) = 0
    nzeros_forecast(:,:,:) = 0


    ishift_in_hours = iday*24
    !print *,'iday, ishift_in_hours,iyyyymmddhh = ', iday, ishift_in_hours,iyyyymmddhh
    CALL updat(iyyyymmddhh,ishift_in_hours,jyyyymmddhh)
    WRITE (cyyyymmddhh,'(i10)') jyyyymmddhh
    !PRINT *,'jyyyymmddhh, cyyyymmddhh = ', jyyyymmddhh, cyyyymmddhh
    
    ! ---- build filename, read in this day's information
    
    infile = TRIM(data_directory) // TRIM(cmodel) // '/' // TRIM(cmodel) // &
        '_D_statistics_data_IC' //cyyyymmddhh // '_' // TRIM(cleade) // 'h.nc'
    print *,'reading ', TRIM(infile)
    INQUIRE (file = infile, exist=exist)
    !PRINT *,'exist = ',exist
    IF (exist) THEN

        CALL check(nf90_open(TRIM(infile),NF90_NOWRITE,ncid))

        ! ---- read in the supplemental locations.  NOTE: they 
        !      are defined on the enlarged domain that extends 
        !      further west and south, and so we need to thin 
        !      down the size of the domain and shift the values
        !      of grid indices.

        cvar = "sumx_analysis"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,sumx_analysis))

        cvar = "sumlnx_analysis"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,sumlnx_analysis))

        cvar = "npositive_analysis"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,npositive_analysis))
    
        cvar = "nzeros_analysis"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,nzeros_analysis))
        
        cvar = "sumx_forecast"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,sumx_forecast))

        cvar = "sumlnx_forecast"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,sumlnx_forecast))

        cvar = "npositive_forecast"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,npositive_forecast))
    
        cvar = "nzeros_forecast"
        CALL check(nf90_inq_varid(ncid,cvar,ivarid))
        CALL check(nf90_get_var(ncid,ivarid,nzeros_forecast))    
    
        CALL check(nf90_close(ncid))
    
        !PRINT *,'npositive_analysis(1:nxa:5, nya/2) = ', &
        !    npositive_analysis(1:nxa:5, nya/2)
        !PRINT *,'nzeros_analysis(1:nxa:5, nya/2) = ', &
        !    nzeros_analysis(1:nxa:5, nya/2)
        !PRINT *,'sumx_analysis(1:nxa:5, nya/2) = ', &
        !    sumx_analysis(1:nxa:5, nya/2)        
        !PRINT *,'sumlnx_analysis(1:nxa:5, nya/2) = ', &
        !    sumlnx_analysis(1:nxa:5, nya/2)  
    
        sumx_analysis_multiday = sumx_analysis_multiday + sumx_analysis
        sumlnx_analysis_multiday = sumlnx_analysis_multiday + sumlnx_analysis
        npositive_analysis_multiday = npositive_analysis_multiday + INT(npositive_analysis)
        nzeros_analysis_multiday = nzeros_analysis_multiday + INT(nzeros_analysis)
    
        sumx_forecast_multiday = sumx_forecast_multiday + sumx_forecast
        sumlnx_forecast_multiday = sumlnx_forecast_multiday + sumlnx_forecast
        npositive_forecast_multiday = npositive_forecast_multiday + INT(npositive_forecast)
        nzeros_forecast_multiday = nzeros_forecast_multiday + INT(nzeros_forecast)
    ELSE
        PRINT *,'This file does not exist.'
    END IF 
END DO

! ---- now, using supplemental location information, build up the final arrays used
!      to calculate D statistic.

sumx_analysis(:,:) = 0.0
sumlnx_analysis(:,:) = 0.0
npositive_analysis_final(:,:) = 0
nzeros_analysis_final(:,:) = 0

sumx_forecast(:,:,:) = 0.0
sumlnx_forecast(:,:,:) = 0.0
npositive_forecast_final(:,:,:) = 0
nzeros_forecast_final(:,:,:) = 0

DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(ixa,jya) .eq. 1) THEN
            DO isupp = 1, nsupplemental(ixa,jya)
                ixn = xlocations(ixa,jya,isupp)
                jyn = ylocations(ixa,jya,isupp)
                IF (conusmask(ixn,jyn) .eq. 1) THEN
                    sumx_analysis(ixa,jya) = sumx_analysis(ixa,jya) + &
                        sumx_analysis_multiday(ixn,jyn)
                    sumlnx_analysis(ixa,jya) = sumlnx_analysis(ixa,jya) + &
                        sumlnx_analysis_multiday(ixn,jyn)
                    npositive_analysis_final(ixa,jya) = npositive_analysis_final(ixa,jya) + &
                        npositive_analysis_multiday(ixn,jyn)
                    nzeros_analysis_final(ixa,jya) = nzeros_analysis_final(ixa,jya) + &
                        nzeros_analysis_multiday(ixn,jyn)
                     
                    DO imem = 1, nens_qmap
                        sumx_forecast(ixa,jya,imem) = sumx_forecast(ixa,jya,imem) + &
                            sumx_forecast_multiday(ixn,jyn,imem)
                        sumlnx_forecast(ixa,jya,imem) = sumlnx_forecast(ixa,jya,imem) + &
                            sumlnx_forecast_multiday(ixn,jyn,imem)
                        npositive_forecast_final(ixa,jya,imem) = &
                            npositive_forecast_final(ixa,jya,imem) + &
                            npositive_forecast_multiday(ixn,jyn,imem)
                        nzeros_forecast_final(ixa,jya,imem) = nzeros_forecast_final(ixa,jya,imem) + &
                            nzeros_forecast_multiday(ixn,jyn,imem)
                    END DO
                ELSE
                    PRINT *,'when processing ixa, jya = ',ixa, jya,' id supplemental outside CONUS',ixn,jyn
                    STOP
                ENDIF
            END DO
        ENDIF
    END DO
END DO   

!PRINT *,'sumx_analysis(:,nya/2) = ',sumx_analysis(:,nya/2)
!PRINT *,'sumlnx_analysis(:,nya/2) = ',sumlnx_analysis(:,nya/2)
!PRINT *,'npositive_analysis_final(:,nya/2) = ', npositive_analysis_final(:,nya/2)
!PRINT *,'nzeros_analysis_final(:,nya/2) = ', nzeros_analysis_final(:,nya/2)

!PRINT *,'sumx_forecast(:,nya/2,1) = ',sumx_forecast(:,nya/2,1)
!PRINT *,'sumlnx_forecast(:,nya/2,1) = ',sumlnx_forecast(:,nya/2,1)
!PRINT *,'npositive_forecast_final(:,nya/2,1) = ', npositive_forecast_final(:,nya/2,1)
!PRINT *,'nzeros_forecast_final(:,nya/2,1) = ', nzeros_forecast_final(:,nya/2,1)
    
! ---- calculate D statistic and from that Gamma distribution parameters and fraction zero.

gamma_shape_qmap_analysis(:,:) = -99.99
gamma_scale_qmap_analysis(:,:) = -99.99
fraction_zero_qmap_analysis(:,:) = -99.99
gamma_shape_qmap_forecast(:,:,:) = -99.99
gamma_scale_qmap_forecast(:,:,:) = -99.99
fraction_zero_qmap_forecast(:,:,:) = -99.99

DO jya = 1, nya
    DO ixa = 1, nxa
        IF (conusmask(ixa,jya) .eq. 1) THEN
            xbar_anal = sumx_analysis(ixa,jya) / REAL(npositive_analysis_final(ixa,jya))
            D_anal = alog(xbar_anal) - sumlnx_analysis(ixa,jya) / &
                REAL(npositive_analysis_final(ixa,jya))    
            rnumer2 = 1. + 4.*D_anal/3.
            gamma_shape_qmap_analysis(ixa,jya) = (1. + SQRT(rnumer2)) / (4.*D_anal)
            gamma_scale_qmap_analysis(ixa,jya) = &
                xbar_anal / gamma_shape_qmap_analysis(ixa,jya)
            fraction_zero_qmap_analysis(ixa,jya) = REAL(nzeros_analysis_final(ixa,jya)) / &
                REAL(nzeros_analysis_final(ixa,jya) + npositive_analysis_final(ixa,jya))
            DO imem = 1, nens_qmap
                xbar_fcst = sumx_forecast(ixa,jya,imem) / REAL(npositive_forecast_final(ixa,jya,imem))
                D_fcst = alog(xbar_fcst) - sumlnx_forecast(ixa,jya,imem)/ &
                    REAL(npositive_forecast_final(ixa,jya,imem))
                rnumer2 = 1. + 4.*D_fcst/3.
                gamma_shape_qmap_forecast(ixa,jya,imem) = (1. + SQRT(rnumer2)) / (4.*D_fcst)
                gamma_scale_qmap_forecast(ixa,jya,imem) = &
                    xbar_fcst / gamma_shape_qmap_forecast(ixa,jya,imem)
                fraction_zero_qmap_forecast(ixa,jya,imem) = REAL(nzeros_forecast_final(ixa,jya,imem)) / &
                    REAL(nzeros_forecast_final(ixa,jya,imem) + npositive_forecast_final(ixa,jya,imem))
            END DO
        ELSE
            gamma_shape_qmap_analysis(ixa,jya) = -99.99
            gamma_scale_qmap_analysis(ixa,jya) = -99.99
            fraction_zero_qmap_analysis(ixa,jya) = -99.99
            gamma_shape_qmap_forecast(ixa,jya,:) = -99.99
            gamma_scale_qmap_forecast(ixa,jya,:) = -99.99
            fraction_zero_qmap_forecast(ixa,jya,:) = -99.99
        ENDIF
    END DO
END DO    

!PRINT *,'sample gamma shape statistic for analysis = ', &
!    gamma_shape_qmap_analysis(1:nxa:5,nya/2)
!PRINT *,'sample gamma scale statistic for analysis = ', &
!    gamma_scale_qmap_analysis(1:nxa:5,nya/2)
!PRINT *,'sample fraction zero statistic for analysis = ', &
!    fraction_zero_qmap_analysis(1:nxa:5,nya/2)
!PRINT *,'sample fraction zero statistic for forecast = ', &
!    fraction_zero_qmap_forecast(1:nxa:5,nya/2,1)

!PRINT *,' in subroutine determine_gamma_parameters_for_quantile_mapping'
!PRINT *,'fraction_zero_qmap_forecast(nxa/2,nya/2,:) = ',&
!    fraction_zero_qmap_forecast(nxa/2,nya/2,:)
    

! ---- deallocate

DEALLOCATE(xlocations, xlocations_largedomain)    
DEALLOCATE(ylocations, ylocations_largedomain)     
DEALLOCATE(nsupplemental, nsupplemental_largedomain) 
DEALLOCATE(sumx_analysis, sumx_analysis_multiday)
DEALLOCATE(sumlnx_analysis, sumlnx_analysis_multiday)
DEALLOCATE(npositive_analysis, nzeros_analysis)
DEALLOCATE(npositive_analysis_multiday, nzeros_analysis_multiday)
DEALLOCATE(npositive_analysis_final, nzeros_analysis_final)
DEALLOCATE(sumx_forecast, sumx_forecast_multiday)
DEALLOCATE(sumlnx_forecast, sumlnx_forecast_multiday)
DEALLOCATE(npositive_forecast, nzeros_forecast)
DEALLOCATE(npositive_forecast_multiday, nzeros_forecast_multiday)
DEALLOCATE(npositive_forecast_final, nzeros_forecast_final)

RETURN
END SUBROUTINE determine_gamma_parameters_for_quantile_mapping
