subroutine cumgam ( x, a, cum, ccum )

!*****************************************************************************80
!
!! CUMGAM evaluates the cumulative incomplete gamma distribution.
!
!  Discussion:
!
!    This routine computes the cumulative distribution function of the
!    incomplete gamma distribution, i.e., the integral from 0 to X of
!
!      (1/GAM(A))*EXP(-T)*T^(A-1) DT
!
!    where GAM(A) is the complete gamma function of A:
!
!      GAM(A) = integral from 0 to infinity of EXP(-T)*T^(A-1) DT
!
!  Author:
!
!    Barry Brown, James Lovato, Kathy Russell
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the upper limit of integration.
!
!    Input, real ( kind = 8 ) A, the shape parameter of the incomplete
!    Gamma distribution.
!
!    Output, real ( kind = 8 ) CUM, CCUM, the incomplete Gamma CDF and
!    complementary CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ccum
  real ( kind = 8 ) cum
  real ( kind = 8 ) x

  if ( x <= 0.0D+00 ) then

    cum = 0.0D+00
    ccum = 1.0D+00

  else

    call gamma_inc ( a, x, cum, ccum, 0 )

  end if

  return
end