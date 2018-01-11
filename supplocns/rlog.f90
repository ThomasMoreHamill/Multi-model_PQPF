function rlog ( x )

!*****************************************************************************80
!
!! RLOG computes X - 1 - LN(X).
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708:
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) RLOG, the value of the function.
!
  implicit none

  real ( kind = 8 ), parameter :: a  =  0.566749439387324D-01
  real ( kind = 8 ), parameter :: b  =  0.456512608815524D-01
  real ( kind = 8 ), parameter :: half = 0.5D+00
  real ( kind = 8 ), parameter :: p0 =  0.333333333333333D+00
  real ( kind = 8 ), parameter :: p1 = -0.224696413112536D+00
  real ( kind = 8 ), parameter :: p2 =  0.620886815375787D-02
  real ( kind = 8 ), parameter :: q1 = -0.127408923933623D+01
  real ( kind = 8 ), parameter :: q2 =  0.354508718369557D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) rlog
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: two =  2.0D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) x

  if ( x < 0.61D+00 ) then

    r = ( x - 0.5D+00 ) - 0.5D+00
    rlog = r - log ( x )

  else if ( x < 1.57D+00 ) then

    if ( x < 0.82D+00 ) then

      u = x - 0.7D+00
      u = u / 0.7D+00
      w1 = a - u * 0.3D+00

    else if ( x < 1.18D+00 ) then

      u = ( x - half ) - half
      w1 = 0.0D+00

    else if ( x < 1.57D+00 ) then

      u = 0.75D+00 * x - 1.0D+00
      w1 = b + u / 3.0D+00

    end if

    r = u / ( u + two )
    t = r * r
    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
    rlog = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

  else if ( 1.57D+00 <= x ) then

    r = ( x - half ) - half
    rlog = r - log ( x )

  end if

  return
end