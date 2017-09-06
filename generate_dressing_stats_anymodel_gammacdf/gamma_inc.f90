subroutine gamma_inc ( a, x, ans, qans, ind )

!*****************************************************************************80
!
!! GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    Alfred Morris,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, X, the arguments of the incomplete
!    gamma ratio.  A and X must be nonnegative.  A and X cannot
!    both be zero.
!
!    Output, real ( kind = 8 ) ANS, QANS.  On normal output,
!    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
!    A or X is negative, or both are 0, or when the answer is
!    computationally indeterminate because A is extremely large
!    and X is very close to A.
!
!    Input, integer ( kind = 4 ) IND, indicates the accuracy request:
!    0, as much accuracy as possible.
!    1, to within 1 unit of the 6-th significant digit,
!    otherwise, to within 1 unit of the 3rd significant digit.
!
!  Local Parameters:
!
!     ALOG10 = LN(10)
!     RT2PIN = 1/SQRT(2*PI)
!     RTPI   = SQRT(PI)
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2n
  real ( kind = 8 ) a2nm1
  real ( kind = 8 ) acc
  real ( kind = 8 ), dimension ( 3 ) :: acc0 = (/ &
    5.0D-15, 5.0D-07, 5.0D-04 /)
  real ( kind = 8 ), parameter :: alog10 = 2.30258509299405D+00
  real ( kind = 8 ) am0
  real ( kind = 8 ) amn
  real ( kind = 8 ) an
  real ( kind = 8 ) an0
  real ( kind = 8 ) ans
  real ( kind = 8 ) apn
  real ( kind = 8 ) b2n
  real ( kind = 8 ) b2nm1
  real ( kind = 8 ) big(3)
  real ( kind = 8 ) c
  real ( kind = 8 ) c0
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) c4
  real ( kind = 8 ) c5
  real ( kind = 8 ) c6
  real ( kind = 8 ) cma
  real ( kind = 8 ) d0(13)
  real ( kind = 8 ) d1(12)
  real ( kind = 8 ) d2(10)
  real ( kind = 8 ) d3(8)
  real ( kind = 8 ) d4(6)
  real ( kind = 8 ) d5(4)
  real ( kind = 8 ) d6(2)
  real ( kind = 8 ) d10
  real ( kind = 8 ) d20
  real ( kind = 8 ) d30
  real ( kind = 8 ) d40
  real ( kind = 8 ) d50
  real ( kind = 8 ) d60
  real ( kind = 8 ) d70
  real ( kind = 8 ) e
  real ( kind = 8 ) e0
  real ( kind = 8 ) e00(3)
  real ( kind = 8 ) error_f
  real ( kind = 8 ) error_fc
  real ( kind = 8 ) g
  real ( kind = 8 ) gam1
  real ( kind = 8 ) gamma
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) iop
  real ( kind = 8 ) j
  real ( kind = 8 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  real ( kind = 8 ) qans
  real ( kind = 8 ) r
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rlog
  real ( kind = 8 ), parameter :: rt2pin = 0.398942280401433D+00
  real ( kind = 8 ) rta
  real ( kind = 8 ), parameter :: rtpi = 1.77245385090552D+00
  real ( kind = 8 ) rtx
  real ( kind = 8 ) s
  real ( kind = 8 ) sum1
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) tol
  real ( kind = 8 ) twoa
  real ( kind = 8 ) u
  real ( kind = 8 ) w
  real ( kind = 8 ) wk(20)
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x00(3)
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  data big(1)/20.0D+00/,big(2)/14.0D+00/,big(3)/10.0D+00/
  data e00(1)/0.25D-03/,e00(2)/0.25D-01/,e00(3)/0.14D+00/
  data x00(1)/31.0D+00/,x00(2)/17.0D+00/,x00(3)/9.7D+00/
  data d0(1)/0.833333333333333D-01/
  data d0(2)/-0.148148148148148D-01/
  data d0(3)/0.115740740740741D-02/,d0(4)/0.352733686067019D-03/
  data d0(5)/-0.178755144032922D-03/,d0(6)/0.391926317852244D-04/
  data d0(7)/-0.218544851067999D-05/,d0(8)/-0.185406221071516D-05/
  data d0(9)/0.829671134095309D-06/,d0(10)/-0.176659527368261D-06/
  data d0(11)/0.670785354340150D-08/,d0(12)/0.102618097842403D-07/
  data d0(13)/-0.438203601845335D-08/
  data d10/-0.185185185185185D-02/,d1(1)/-0.347222222222222D-02/
  data d1(2)/0.264550264550265D-02/,d1(3)/-0.990226337448560D-03/
  data d1(4)/0.205761316872428D-03/,d1(5)/-0.401877572016461D-06/
  data d1(6)/-0.180985503344900D-04/,d1(7)/0.764916091608111D-05/
  data d1(8)/-0.161209008945634D-05/,d1(9)/0.464712780280743D-08/
  data d1(10)/0.137863344691572D-06/,d1(11)/-0.575254560351770D-07/
  data d1(12)/0.119516285997781D-07/
  data d20/0.413359788359788D-02/,d2(1)/-0.268132716049383D-02/
  data d2(2)/0.771604938271605D-03/,d2(3)/0.200938786008230D-05/
  data d2(4)/-0.107366532263652D-03/,d2(5)/0.529234488291201D-04/
  data d2(6)/-0.127606351886187D-04/,d2(7)/0.342357873409614D-07/
  data d2(8)/0.137219573090629D-05/,d2(9)/-0.629899213838006D-06/
  data d2(10)/0.142806142060642D-06/
  data d30/0.649434156378601D-03/,d3(1)/0.229472093621399D-03/
  data d3(2)/-0.469189494395256D-03/,d3(3)/0.267720632062839D-03/
  data d3(4)/-0.756180167188398D-04/,d3(5)/-0.239650511386730D-06/
  data d3(6)/0.110826541153473D-04/,d3(7)/-0.567495282699160D-05/
  data d3(8)/0.142309007324359D-05/
  data d40/-0.861888290916712D-03/,d4(1)/0.784039221720067D-03/
  data d4(2)/-0.299072480303190D-03/,d4(3)/-0.146384525788434D-05/
  data d4(4)/0.664149821546512D-04/,d4(5)/-0.396836504717943D-04/
  data d4(6)/0.113757269706784D-04/
  data d50/-0.336798553366358D-03/,d5(1)/-0.697281375836586D-04/
  data d5(2)/0.277275324495939D-03/,d5(3)/-0.199325705161888D-03/
  data d5(4)/0.679778047793721D-04/
  data d60/0.531307936463992D-03/,d6(1)/-0.592166437353694D-03/
  data d6(2)/0.270878209671804D-03/
  data d70 / 0.344367606892378D-03/

  e = epsilon ( 1.0D+00 )

  if ( a < 0.0D+00 .or. x < 0.0D+00 ) then
    ans = 2.0D+00
    return
  end if

  if ( a == 0.0D+00 .and. x == 0.0D+00 ) then
    ans = 2.0D+00
    return
  end if

  if ( a * x == 0.0D+00 ) then
    if ( x <= a ) then
      ans = 0.0D+00
      qans = 1.0D+00
    else
      ans = 1.0D+00
      qans = 0.0D+00
    end if
    return
  end if

  iop = ind + 1
  if ( iop /= 1 .and. iop /= 2 ) iop = 3
  acc = max ( acc0(iop), e )
  e0 = e00(iop)
  x0 = x00(iop)
!
!  Select the appropriate algorithm.
!
  if ( 1.0D+00 <= a ) then
    go to 10
  end if

  if ( a == 0.5D+00 ) then
    go to 390
  end if

  if ( x < 1.1D+00 ) then
    go to 160
  end if

  t1 = a * log ( x ) - x
  u = a * exp ( t1 )

  if ( u == 0.0D+00 ) then
    ans = 1.0D+00
    qans = 0.0D+00
    return
  end if

  r = u * ( 1.0D+00 + gam1 ( a ) )
  go to 250

   10 continue

  if ( big(iop) <= a ) then
    go to 30
  end if

  if ( x < a .or. x0 <= x ) then
    go to 20
  end if

  twoa = a + a
  m = int ( twoa )

  if ( twoa == real ( m, kind = 8 ) ) then
    i = m / 2
    if ( a == real ( i, kind = 8 ) ) then
      go to 210
    end if
    go to 220
  end if

   20 continue

  t1 = a * log ( x ) - x
  r = exp ( t1 ) / gamma ( a )
  go to 40

   30 continue

  l = x / a

  if ( l == 0.0D+00 ) then
    ans = 0.0D+00
    qans = 1.0D+00
    return
  end if

  s = 0.5D+00 + ( 0.5D+00 - l )
  z = rlog ( l )
  if ( 700.0D+00 / a <= z ) then
    go to 410
  end if

  y = a * z
  rta = sqrt ( a )

  if ( abs ( s ) <= e0 / rta ) then
    go to 330
  end if

  if ( abs ( s ) <= 0.4D+00 ) then
    go to 270
  end if

  t = ( 1.0D+00 / a )**2
  t1 = ((( 0.75D+00 * t - 1.0D+00 ) * t + 3.5D+00 ) &
    * t - 105.0D+00 ) / ( a * 1260.0D+00 )
  t1 = t1 - y
  r = rt2pin * rta * exp ( t1 )

40    continue

  if ( r == 0.0D+00 ) then
    if ( x <= a ) then
      ans = 0.0D+00
      qans = 1.0D+00
    else
      ans = 1.0D+00
      qans = 0.0D+00
    end if
    return
  end if

  if ( x <= max ( a, alog10 ) ) then
    go to 50
  end if

  if ( x < x0 ) then
    go to 250
  end if

  go to 100
!
!  Taylor series for P/R.
!
50    continue

  apn = a + 1.0D+00
  t = x / apn
  wk(1) = t

  n = 20

  do i = 2, 20
    apn = apn + 1.0D+00
    t = t * ( x / apn )
    if ( t <= 1.0D-03 ) then
      n = i
      exit
    end if
    wk(i) = t
  end do

  sum1 = t

  tol = 0.5D+00 * acc

  do

    apn = apn + 1.0D+00
    t = t * ( x / apn )
    sum1 = sum1 + t

    if ( t <= tol ) then
      exit
    end if

  end do

  n_max = n - 1
  do m = 1, n_max
    n = n - 1
    sum1 = sum1 + wk(n)
  end do

  ans = ( r / a ) * ( 1.0D+00 + sum1 )
  qans = 0.5D+00 + ( 0.5D+00 - ans )
  return
!
!  Asymptotic expansion.
!
  100 continue

  amn = a - 1.0D+00
  t = amn / x
  wk(1) = t

  n = 20

  do i = 2, 20
    amn = amn - 1.0D+00
    t = t * ( amn / x )
    if ( abs ( t ) <= 1.0D-03 ) then
      n = i
      exit
    end if
    wk(i) = t
  end do

  sum1 = t

  do

    if ( abs ( t ) <= acc ) then
      exit
    end if

    amn = amn - 1.0D+00
    t = t * ( amn / x )
    sum1 = sum1 + t

  end do

  n_max = n - 1
  do m = 1, n_max
    n = n - 1
    sum1 = sum1 + wk(n)
  end do
  qans = ( r / x ) * ( 1.0D+00 + sum1 )
  ans = 0.5D+00 + ( 0.5D+00 - qans )
  return
!
!  Taylor series for P(A,X)/X**A
!
  160 continue

  an = 3.0D+00
  c = x
  sum1 = x / ( a + 3.0D+00 )
  tol = 3.0D+00 * acc / ( a + 1.0D+00 )

  do

    an = an + 1.0D+00
    c = -c * ( x / an )
    t = c / ( a + an )
    sum1 = sum1 + t

    if ( abs ( t ) <= tol ) then
      exit
    end if

  end do

  j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 / &
    ( a +  2.0D+00  ) ) * x + 1.0D+00 &
    / ( a + 1.0D+00 ) )

  z = a * log ( x )
  h = gam1 ( a )
  g = 1.0D+00 + h

  if ( x < 0.25D+00 ) then
    go to 180
  end if

  if ( a < x / 2.59D+00 ) then
    go to 200
  end if

  go to 190

  180 continue

  if ( -0.13394D+00 < z ) then
    go to 200
  end if

  190 continue

  w = exp ( z )
  ans = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
  qans = 0.5D+00 + ( 0.5D+00 - ans )
  return

200   continue

  l = rexp ( z )
  w = 0.5D+00 + ( 0.5D+00 + l )
  qans = ( w * j - l ) * g - h

  if ( qans < 0.0D+00 ) then
    ans = 1.0D+00
    qans = 0.0D+00
    return
  end if

  ans = 0.5D+00 + ( 0.5D+00 - qans )
  return
!
!  Finite sums for Q when 1 <= A and 2*A is an integer.
!
210   continue

  sum1 = exp ( - x )
  t = sum1
  n = 1
  c = 0.0D+00
  go to 230

220   continue

  rtx = sqrt ( x )
  sum1 = error_fc ( 0, rtx )
  t = exp ( -x ) / ( rtpi * rtx )
  n = 0
  c = -0.5D+00

  230 continue

  do while ( n /= i )
    n = n + 1
    c = c + 1.0D+00
    t = ( x * t ) / c
    sum1 = sum1 + t
  end do

  240 continue

  qans = sum1
  ans = 0.5D+00 + ( 0.5D+00 - qans )
  return
!
!  Continued fraction expansion.
!
250   continue

  tol = max ( 5.0D+00 * e, acc )
  a2nm1 = 1.0D+00
  a2n = 1.0D+00
  b2nm1 = x
  b2n = x + ( 1.0D+00 - a )
  c = 1.0D+00

  do

    a2nm1 = x * a2n + c * a2nm1
    b2nm1 = x * b2n + c * b2nm1
    am0 = a2nm1 / b2nm1
    c = c + 1.0D+00
    cma = c - a
    a2n = a2nm1 + cma * a2n
    b2n = b2nm1 + cma * b2n
    an0 = a2n / b2n

    if ( abs ( an0 - am0 ) < tol * an0 ) then
      exit
    end if

  end do

  qans = r * an0
  ans = 0.5D+00 + ( 0.5D+00 - qans )
  return
!
!  General Temme expansion.
!
270   continue

  if ( abs ( s ) <= 2.0D+00 * e .and. 3.28D-03 < a * e * e ) then
    ans =  2.0D+00
    return
  end if

  c = exp ( - y )
  w = 0.5D+00 * error_fc ( 1, sqrt ( y ) )
  u = 1.0D+00 / a
  z = sqrt ( z + z )

  if ( l < 1.0D+00 ) then
    z = -z
  end if

  if ( iop < 2 ) then

    if ( abs ( s ) <= 1.0D-03 ) then

      c0 = ((((((     &
              d0(7)   &
        * z + d0(6) ) &
        * z + d0(5) ) &
        * z + d0(4) ) &
        * z + d0(3) ) &
        * z + d0(2) ) &
        * z + d0(1) ) &
        * z - 1.0D+00 / 3.0D+00

      c1 = (((((      &
              d1(6)   &
        * z + d1(5) ) &
        * z + d1(4) ) &
        * z + d1(3) ) &
        * z + d1(2) ) &
        * z + d1(1) ) &
        * z + d10

      c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

      c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

      c4 = ( d4(2) * z + d4(1) ) * z + d40
      c5 = ( d5(2) * z + d5(1) ) * z + d50
      c6 = d6(1) * z + d60

      t = (((((( d70 &
        * u + c6 ) &
        * u + c5 ) &
        * u + c4 ) &
        * u + c3 ) &
        * u + c2 ) &
        * u + c1 ) &
        * u + c0

    else

      c0 = (((((((((((( &
              d0(13)   &
        * z + d0(12) ) &
        * z + d0(11) ) &
        * z + d0(10) ) &
        * z + d0(9)  ) &
        * z + d0(8)  ) &
        * z + d0(7)  ) &
        * z + d0(6)  ) &
        * z + d0(5)  ) &
        * z + d0(4)  ) &
        * z + d0(3)  ) &
        * z + d0(2)  ) &
        * z + d0(1)  ) &
        * z - 1.0D+00 / 3.0D+00

      c1 = ((((((((((( &
                d1(12) &
          * z + d1(11) &
        ) * z + d1(10) &
        ) * z + d1(9)  &
        ) * z + d1(8)  &
        ) * z + d1(7)  &
        ) * z + d1(6)  &
        ) * z + d1(5)  &
        ) * z + d1(4)  &
        ) * z + d1(3)  &
        ) * z + d1(2)  &
        ) * z + d1(1)  &
        ) * z + d10

      c2 = ((((((((( &
                d2(10) &
          * z + d2(9) &
        ) * z + d2(8) &
        ) * z + d2(7) &
        ) * z + d2(6) &
        ) * z + d2(5) &
        ) * z + d2(4) &
        ) * z + d2(3) &
        ) * z + d2(2) &
        ) * z + d2(1) &
        ) * z + d20

      c3 = ((((((( &
                d3(8) &
          * z + d3(7) &
        ) * z + d3(6) &
        ) * z + d3(5) &
        ) * z + d3(4) &
        ) * z + d3(3) &
        ) * z + d3(2) &
        ) * z + d3(1) &
        ) * z + d30

      c4 = ((((( d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z + d40

      c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50

      c6 = ( d6(2) * z + d6(1) ) * z + d60

      t = ((((((   &
            d70    &
        * u + c6 ) &
        * u + c5 ) &
        * u + c4 ) &
        * u + c3 ) &
        * u + c2 ) &
        * u + c1 ) &
        * u + c0

    end if

  else if ( iop == 2 ) then

    c0 = (((((      &
            d0(6)   &
      * z + d0(5) ) &
      * z + d0(4) ) &
      * z + d0(3) ) &
      * z + d0(2) ) &
      * z + d0(1) ) &
      * z - 1.0D+00 / 3.0D+00

    c1 = ((( d1(4) * z + d1(3) ) * z + d1(2) ) * z + d1(1) ) * z + d10
    c2 = d2(1) * z + d20
    t = ( c2 * u + c1 ) * u + c0

  else if ( 2 < iop ) then

    t = (( d0(3) * z + d0(2) ) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00

  end if

310   continue

  if ( 1.0D+00 <= l ) then
    qans = c * ( w + rt2pin * t / rta )
    ans = 0.5D+00 + ( 0.5D+00 - qans )
  else
    ans = c * ( w - rt2pin * t / rta )
    qans = 0.5D+00 + ( 0.5D+00 - ans )
  end if

  return
!
!  Temme expansion for L = 1
!
  330 continue

  if ( 3.28D-03 < a * e * e ) then
    ans =  2.0D+00
    return
  end if

  c = 0.5D+00 + ( 0.5D+00 - y )
  w = ( 0.5D+00 - sqrt ( y ) &
    * ( 0.5D+00 &
    + ( 0.5D+00 - y / 3.0D+00 ) ) / rtpi ) / c
  u = 1.0D+00 / a
  z = sqrt ( z + z )

  if ( l < 1.0D+00 ) then
    z = -z
  end if

  if ( iop < 2 ) then

    c0 = ((((((     &
            d0(7)   &
      * z + d0(6) ) &
      * z + d0(5) ) &
      * z + d0(4) ) &
      * z + d0(3) ) &
      * z + d0(2) ) &
      * z + d0(1) ) &
      * z - 1.0D+00 / 3.0D+00

    c1 = (((((      &
            d1(6)   &
      * z + d1(5) ) &
      * z + d1(4) ) &
      * z + d1(3) ) &
      * z + d1(2) ) &
      * z + d1(1) ) &
      * z + d10

    c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

    c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

    c4 = ( d4(2) * z + d4(1) ) * z + d40
    c5 = ( d5(2) * z + d5(1) ) * z + d50
    c6 = d6(1) * z + d60

    t = (((((( d70 &
      * u + c6 ) &
      * u + c5 ) &
      * u + c4 ) &
      * u + c3 ) &
      * u + c2 ) &
      * u + c1 ) &
      * u + c0

  else if ( iop == 2 ) then

    c0 = ( d0(2) * z + d0(1) ) * z - 1.0D+00 / 3.0D+00
    c1 = d1(1) * z + d10
    t = ( d20 * u + c1 ) * u + c0

  else if ( 2 < iop ) then

    t = d0(1) * z - 1.0D+00 / 3.0D+00

  end if

  go to 310
!
!  Special cases
!
  390 continue

  if ( x < 0.25D+00 ) then
    ans = error_f ( sqrt ( x ) )
    qans = 0.5D+00 + ( 0.5D+00 - ans )
  else
    qans = error_fc ( 0, sqrt ( x ) )
    ans = 0.5D+00 + ( 0.5D+00 - qans )
  end if

  return

  410 continue

  if ( abs ( s ) <= 2.0D+00 * e ) then
    ans =  2.0D+00
    return
  end if

  if ( x <= a ) then
    ans = 0.0D+00
    qans = 1.0D+00
  else
    ans = 1.0D+00
    qans = 0.0D+00
  end if

  return
end