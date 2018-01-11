SUBROUTINE read_HXLR_coefficients(data_directory, cmodel, cyyyymmddhh, &
    cleade, b0_mean, b1_mean, b0_spread, b1_spread)
    
CHARACTER*(*), INTENT(IN) :: data_directory, cmodel, cyyyymmddhh, cleade
REAL, INTENT(OUT) :: b0_mean, b1_mean, b0_spread, b1_spread

CHARACTER*256 infile


! ----- set the regression coefficient file name, read in regression
!       coefficients.  For the time being, the regression coefficients
!       are hardwired.

!b0_mean = -0.495270 ! for 0.33 power transformed, transform before spread and mean calc
!b1_mean = 1.173933
!b0_spread = -0.264854
!b1_spread = 0.431978

b0_mean = -0.897879 ! for 0.33 power transformed, transform after spread and mean calc
b1_mean = 1.308752
b0_spread = -0.6557512
b1_spread = 0.0422231

!b0_mean = -0.614893 ! with no power transform
!b1_mean = 0.825213
!b0_spread = 0.4302449
!b1_spread = 0.5675933

b0_mean = -1.477 ! with no power transform
b1_mean = 0.821652
b0_spread = 0.8936721
b1_spread = 0.0

RETURN
END SUBROUTINE read_HXLR_coefficients