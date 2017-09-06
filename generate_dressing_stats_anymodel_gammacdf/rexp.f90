function rexp ( x )

!*****************************************************************************80
!
!! REXP evaluates the function EXP(X) - 1.
!
!  Modified:
!
!    09 December 1999
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
!    Output, real ( kind = 8 ) REXP, the value of EXP(X)-1.
!
  implicit none

  real ( kind = 8 ), parameter :: p1 =  0.914041914819518D-09
  real ( kind = 8 ), parameter :: p2 =  0.238082361044469D-01
  real ( kind = 8 ), parameter :: q1 = -0.499999999085958D+00
  real ( kind = 8 ), parameter :: q2 =  0.107141568980644D+00
  real ( kind = 8 ), parameter :: q3 = -0.119041179760821D-01
  real ( kind = 8 ), parameter :: q4 =  0.595130811860248D-03
  real ( kind = 8 ) rexp
  real ( kind = 8 ) w
  real ( kind = 8 ) x

  if ( abs ( x ) <= 0.15D+00 ) then

    rexp = x * ( ( ( p2 * x + p1 ) * x + 1.0D+00 ) &
      / ( ( ( ( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00 ) )

  else

    w = exp ( x )

    if ( x <= 0.0D+00 ) then
      rexp = ( w - 0.5D+00 ) - 0.5D+00
    else
      rexp = w * ( 0.5D+00 + ( 0.5D+00 - 1.0D+00 / w ) )
    end if

  end if

  return
end