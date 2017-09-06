
SUBROUTINE read_prob_forecasts_weighted(cyyyymmddhh, &
	cleade, cthresh, nxa, nya, prob_forecast_raw, &
    prob_forecast_raw_NCEP, prob_forecast_raw_CMC, &
    prob_forecast_raw_ECMWF, prob_forecast_weighted, &
    prob_forecast_weighted_NCEP, prob_forecast_weighted_CMC, &
    prob_forecast_weighted_ECMWF, prob_forecast_weighted_avg, &
    climo_prob, rlonsa, rlatsa, conusmask)
	
! FILE: fortran_routine.f90
!  f2py -c -m read_prob_forecasts_weighted read_prob_forecasts_weighted.f90
!  this file will have a collection of fortran routines that can be interfaced to python
! ====================================================================================================

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*10, INTENT(IN) :: cyyyymmddhh
CHARACTER*3, INTENT(IN) :: cleade
CHARACTER*(*), INTENT(IN) :: cthresh
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: prob_forecast_raw, &
    prob_forecast_raw_NCEP, prob_forecast_raw_CMC,&
    prob_forecast_raw_ECMWF, prob_forecast_weighted, &
    prob_forecast_weighted_NCEP, prob_forecast_weighted_CMC, &
    prob_forecast_weighted_ECMWF, prob_forecast_weighted_avg, &
    climo_prob, rlonsa, rlatsa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
CHARACTER*120 infile, input_directory

REAL, DIMENSION(nxa,nya) :: prob_forecast_junk

! f2py intent(in) cyyyymmddhh, cleade, cthresh, nxa, nya
! f2py intent(out) conusmask, prob_forecast_raw, prob_forecast_raw_NCEP
! f2py intent(out) prob_forecast_raw_CMC, prob_forecast_raw_ECMWF
! f2py intent(out) prob_forecast_weighted, prob_forecast_weighted_NCEP
! f2py intent(out) prob_forecast_weighted_CMC, prob_forecast_weighted_ECMWF
! f2py intent(out) prob_forecast_weighted_avg
! f2py intent(out) climo_prob, rlonsa, rlatsa, prob_forecast_qmap_x9

! f2py depend(nxa,nya) conusmask, prob_forecast_raw, prob_forecast_raw_NCEP
! f2py depend(nxa,nya) prob_forecast_raw_CMC, prob_forecast_raw_ECMWF
! f2py depend(nxa,nya) prob_forecast_weighted, prob_forecast_weighted_NCEP 
! f2py depend(nxa,nya) prob_forecast_weighted_CMC
! f2py depend(nxa,nya) prob_forecast_weighted_ECMWF, prob_forecast_weighted_avg
! f2py depend(nxa,nya) climo_prob, rlonsa, rlatsa, prob_forecast_qmap_x9

! ---- read the PROB forecast with quantile mapping center point only

infile = '/Users/thamill/precip/ecmwf_data/ECMWF_NCEP_CMC_'//&
    TRIM(cleade)//'h_IC'//cyyyymmddhh//'_thresh'//TRIM(cthresh)//'.dat'

PRINT *,'reading data from ',TRIM(infile)
OPEN (unit=1,file=infile,status='old',form='unformatted')

READ (1)  nxain ,nyain
READ (1)  prob_forecast_raw
READ (1)  prob_forecast_raw_CMC
READ (1)  prob_forecast_raw_NCEP
READ (1)  prob_forecast_raw_ECMWF
READ (1)  prob_forecast_weighted
READ (1)  prob_forecast_weighted_CMC
READ (1)  prob_forecast_weighted_NCEP
READ (1)  prob_forecast_weighted_ECMWF
READ (1)  prob_forecast_weighted_avg
READ (1)  climo_prob
READ (1)  rlonsa
READ (1)  rlatsa
READ (1)  conusmask
CLOSE (1) 

RETURN
END SUBROUTINE read_prob_forecasts_weighted

