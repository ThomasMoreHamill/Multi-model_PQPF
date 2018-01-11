function exparg ( l )

!*****************************************************************************80
!
!! EXPARG returns the largest or smallest legal argument for EXP.
!
!  Discussion:
!
!    Only an approximate limit for the argument of EXP is desired.
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
!    Input, integer ( kind = 4 ) L, indicates which limit is desired.
!    If L = 0, then the largest positive argument for EXP is desired.
!    Otherwise, the largest negative argument for EXP for which the
!    result is nonzero is desired.
!
!    Output, real ( kind = 8 ) EXPARG, the desired value.
!
  implicit none

  integer ( kind = 4 ) b
  real ( kind = 8 ) exparg
  integer ( kind = 4 ) ipmpar
  integer ( kind = 4 ) l
  real ( kind = 8 ) lnb
  integer ( kind = 4 ) m
!
!  Get the arithmetic base.
!
  b = ipmpar(4)
!
!  Compute the logarithm of the arithmetic base.
!
  if ( b == 2 ) then
    lnb = 0.69314718055995D+00
  else if ( b == 8 ) then
    lnb = 2.0794415416798D+00
  else if ( b == 16 ) then
    lnb = 2.7725887222398D+00
  else
    lnb = log ( real ( b, kind = 8 ) )
  end if

  if ( l /= 0 ) then
    m = ipmpar(9) - 1
    exparg = 0.99999D+00 * ( m * lnb )
  else
    m = ipmpar(10)
    exparg = 0.99999D+00 * ( m * lnb )
  end if

  return
end