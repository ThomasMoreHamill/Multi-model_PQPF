function error_fc ( ind, x )

!*****************************************************************************80
!
!! ERROR_FC evaluates the complementary error function.
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
!    Input, integer ( kind = 4 ) IND, chooses the scaling.
!    If IND is nonzero, then the value returned has been multiplied by
!    EXP(X*X).
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) ERROR_FC, the value of the complementary
!    error function.
!
  implicit none

  real ( kind = 8 ), dimension ( 5 ) :: a = (/ &
     0.771058495001320D-04,  -0.133733772997339D-02, &
     0.323076579225834D-01,   0.479137145607681D-01, &
     0.128379167095513D+00 /)
  real ( kind = 8 ) ax
  real ( kind = 8 ), dimension(3) :: b = (/ &
    0.301048631703895D-02, &
    0.538971687740286D-01, &
    0.375795757275549D+00 /)
  real ( kind = 8 ) bot
  real ( kind = 8 ), parameter :: c = 0.564189583547756D+00
  real ( kind = 8 ) e
  real ( kind = 8 ) error_fc
  real ( kind = 8 ) exparg
  integer ( kind = 4 ) ind
  real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.36864857382717D-07, 5.64195517478974D-01, &
     7.21175825088309D+00, 4.31622272220567D+01, &
     1.52989285046940D+02, 3.39320816734344D+02, &
     4.51918953711873D+02, 3.00459261020162D+02 /)
  real ( kind = 8 ), dimension ( 8 ) :: q = (/  &
    1.00000000000000D+00, 1.27827273196294D+01, &
    7.70001529352295D+01, 2.77585444743988D+02, &
    6.38980264465631D+02, 9.31354094850610D+02, &
    7.90950925327898D+02, 3.00459260956983D+02 /)
  real ( kind = 8 ), dimension ( 5 ) :: r = (/ &
    2.10144126479064D+00, 2.62370141675169D+01, &
    2.13688200555087D+01, 4.65807828718470D+00, &
    2.82094791773523D-01 /)
  real ( kind = 8 ), dimension ( 4 ) :: s = (/ &
    9.41537750555460D+01, 1.87114811799590D+02, &
    9.90191814623914D+01, 1.80124575948747D+02 /)
  real ( kind = 8 ) t
  real ( kind = 8 ) top
  real ( kind = 8 ) w
  real ( kind = 8 ) x
!
!  ABS ( X ) <= 0.5
!
  ax = abs ( x )

  if ( ax <= 0.5D+00 ) then

    t = x * x

    top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) &
      + 1.0D+00

    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00

    error_fc = 0.5D+00 + ( 0.5D+00 &
      - x * ( top / bot ) )

    if ( ind /= 0 ) then
      error_fc = exp ( t ) * error_fc
    end if

    return

  end if
!
!  0.5 < abs ( X ) <= 4
!
  if ( ax <= 4.0D+00 ) then

    top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
      + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)

    bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
      + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)

    error_fc = top / bot
!
!  4 < ABS ( X )
!
  else

    if ( x <= -5.6D+00 ) then

      if ( ind == 0 ) then
        error_fc =  2.0D+00
      else
        error_fc =  2.0D+00  * exp ( x * x )
      end if

      return

    end if

    if ( ind == 0 ) then

      if ( 100.0D+00 < x ) then
        error_fc = 0.0D+00
        return
      end if

      if ( -exparg ( 1 ) < x * x ) then
        error_fc = 0.0D+00
        return
      end if

    end if

    t = ( 1.0D+00 / x )**2

    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
      + 1.0D+00

    error_fc = ( c - t * top / bot ) / ax

  end if
!
!  Final assembly.
!
  if ( ind /= 0 ) then

    if ( x < 0.0D+00 ) then
      error_fc =  2.0D+00  * exp ( x * x ) - error_fc
    end if

  else

    w = x * x
    t = w
    e = w - t
    error_fc = (( 0.5D+00 &
      + ( 0.5D+00 - e ) ) * exp ( - t ) ) * error_fc

    if ( x < 0.0D+00 ) then
      error_fc =  2.0D+00  - error_fc
    end if

  end if

  return
end