!-----------------------------------------------------------------------------
!
! Evaluation of the gamma function
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000-2013 The R Core Team
!  Copyright (C) 2002-2004 The R Foundation
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!  SYNOPSIS
!
!    #include <Rmath.h>
!    double gammafn(double x);
!
!  DESCRIPTION
!
!    This function computes the value of the gamma function.
!
!  NOTES
!
!    This function is a translation into C of a Fortran subroutine
!    by W. Fullerton of Los Alamos Scientific Laboratory.
!    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
!
!    The accuracy of this routine compares (very) favourably
!    with those of the Sun Microsystems portable mathematical
!    library.
!
!    MM specialized the case of  n!  for n < 50 - for even better precision


FUNCTION gammafn (x)  RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  INTEGER, PARAMETER :: ngam = 22

  DOUBLE PRECISION, PARAMETER :: POS_INF = 9999.0d99
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -9999.0d99
  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION, PARAMETER :: M_PI = 3.141592653589793238462643383280d0
  DOUBLE PRECISION, PARAMETER :: M_LN_SQRT_2PI = 0.918938533204672741780329736406d0

  DOUBLE PRECISION, PARAMETER :: xmin = -170.5674972726612d0
  DOUBLE PRECISION, PARAMETER :: xmax = 171.61447887182298d0
  DOUBLE PRECISION, PARAMETER :: xsml = 2.2474362225598545d-308
  DOUBLE PRECISION, PARAMETER :: dxrel = 1.490116119384765696d-8

  DOUBLE PRECISION, PARAMETER :: gamcs(42) = (/ &
	+.8571195590989331421920062399942d-2, &
	+.4415381324841006757191315771652d-2, &
	+.5685043681599363378632664588789d-1, &
	-.4219835396418560501012500186624d-2, &
	+.1326808181212460220584006796352d-2, &
	-.1893024529798880432523947023886d-3, &
	+.3606925327441245256578082217225d-4, &
	-.6056761904460864218485548290365d-5, &
	+.1055829546302283344731823509093d-5, &
	-.1811967365542384048291855891166d-6, &
	+.3117724964715322277790254593169d-7, &
	-.5354219639019687140874081024347d-8, &
	+.9193275519859588946887786825940d-9, &
	-.1577941280288339761767423273953d-9, &
	+.2707980622934954543266540433089d-10, &
	-.4646818653825730144081661058933d-11, &
	+.7973350192007419656460767175359d-12, &
	-.1368078209830916025799499172309d-12, &
	+.2347319486563800657233471771688d-13, &
	-.4027432614949066932766570534699d-14, &
	+.6910051747372100912138336975257d-15, &
	-.1185584500221992907052387126192d-15, &
	+.2034148542496373955201026051932d-16, &
	-.3490054341717405849274012949108d-17, &
	+.5987993856485305567135051066026d-18, &
	-.1027378057872228074490069778431d-18, &
	+.1762702816060529824942759660748d-19, &
	-.3024320653735306260958772112042d-20, &
	+.5188914660218397839717833550506d-21, &
	-.8902770842456576692449251601066d-22, &
	+.1527474068493342602274596891306d-22, &
	-.2620731256187362900257328332799d-23, &
	+.4496464047830538670331046570666d-24, &
	-.7714712731336877911703901525333d-25, &
	+.1323635453126044036486572714666d-25, &
	-.2270999412942928816702313813333d-26, &
	+.3896418998003991449320816639999d-27, &
	-.6685198115125953327792127999999d-28, &
	+.1146998663140024384347613866666d-28, &
	-.1967938586345134677295103999999d-29, &
	+.3376448816585338090334890666666d-30, &
	-.5793070335782135784625493333333d-31 &
  /)

  INTEGER :: i, n
  DOUBLE PRECISION :: y, sinpiy, val, fn_val
  DOUBLE PRECISION :: stirlerr, chebyshev_eval, lgammacor

  ! If the argument is exactly zero or a negative integer then return NaN
  IF ( x.eq.0.0d0 ) THEN
     PRINT *, ' Warning! Gamma function evaluated at zero! '
     fn_val = NaN
     RETURN
  END IF
  IF ( x.lt.0.0d0 .and. x.eq.AINT(x) ) THEN
     PRINT *, ' Warning! Gamma function evaluated at negative integer! '
     fn_val = NaN
     RETURN
  END IF

  y = abs(x)

  IF ( y.le.10.0d0 ) THEN

     ! Compute gamma(x) for -10 <= x <= 10
     ! Reduce the interval and find gamma(1 + y) for 0 <= y < 1
     ! first of all.

     n = INT(x)
     IF (x .lt. 0.0d0) n = n-1

     y = x - DBLE(n)   ! n = floor(x) ==> y in [ 0, 1 )
     n = n-1

     val = chebyshev_eval (y * 2.0d0 - 1.0d0, gamcs, ngam) + 0.9375d0
     IF ( n.eq.0) THEN
        fn_val = val
        RETURN   ! x = 1.dddd = 1+y
     END IF

     IF (n .lt. 0) THEN
        ! compute gamma(x) for -10 <= x < 1

        IF (x .lt. -0.5d0 .and. abs((x-AINT(x-0.5d0)) / x) < dxrel) THEN
           PRINT *, ' Warning! Argument of gamma function is too near a negative integer, '
           PRINT *, '  precision is less than half of usual precision! '
        END IF

        ! The argument is so close to 0 that the result would overflow.
        IF (y .lt. xsml) THEN
           PRINT *, ' Warning! Argument of gamma function is too near a negative integer, '
           PRINT *, '  precision is less than half of usual precision! '
           IF (x .gt. 0.0d0) THEN
              fn_val = POS_INF
           ELSE
              fn_val = NEG_INF
           END IF
           RETURN
        END IF

        n = -n

        DO i = 0, n-1
           val = val / (x + DBLE(i))
        END DO
        fn_val = val
        RETURN

     ELSE
        ! gamma(x) for 2 <= x <= 10

        DO i = 1, n
           val = val * (y + DBLE(i))
        END DO
        fn_val = val
        RETURN
     END IF

  ELSE
     ! gamma(x) for y = |x| > 10

     IF (x .gt. xmax) THEN
        PRINT *, ' Warning! Overflow when evaluating the gamma function! '
        fn_val = POS_INF
        RETURN
     END IF

     IF (x .lt. xmin) THEN
        PRINT *, ' Warning! Underflow when evaluating the gamma function! '
        fn_val = 0.0d0
        RETURN
     END IF

     IF (y .lt. 50.0d0 .and. y .eq. AINT(y)) THEN  ! compute (n - 1)!
        val = 1.0d0
        DO i = 2, (INT(y)-1)
           val = val * DBLE(i)
        END DO
        fn_val = val
        RETURN

     ELSE  ! normal case
        IF (2.0d0*y .eq. AINT(2.0d0*y)) THEN
           val = (y-0.5d0) * log(y) - y + M_LN_SQRT_2PI + stirlerr(y)
        ELSE
           val = (y-0.5d0) * log(y) - y + M_LN_SQRT_2PI + lgammacor(y)
        END IF
        fn_val = exp(val)

	IF (x .gt. 0) RETURN

        IF (abs((x-AINT(x-0.5d0))/x) .lt. dxrel) THEN
           PRINT *, ' Warning! Argument of gamma function is too near a negative integer, '
           PRINT *, '  precision is less than half of usual precision! '
        END IF

!        sinpiy = sinpi(y)
        sinpiy = abs(sin(y*M_PI))

	IF (sinpiy .eq. 0.0d0) THEN  ! Negative integer arg - overflow
           PRINT *, ' Warning! Overflow when evaluating the gamma function! '
           fn_val = POS_INF
           RETURN
        END IF

        fn_val = -M_PI / (y * sinpiy * fn_val)
        RETURN

     END IF
  END IF

END FUNCTION gammafn





!  AUTHOR
!    Catherine Loader, catherine@research.bell-labs.com.
!    October 23, 2000.
!
!  Merge in to R:
!	Copyright (C) 2000, The R Core Team
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!
!  DESCRIPTION
!
!    Computes the log of the error term in Stirling's formula.
!      For n > 15, uses the series 1/12n - 1/360n^3 + ...
!      For n <=15, integers or half-integers, uses stored values.
!      For other n < 15, uses lgamma directly (don't use this to
!        write lgamma!)
!
! Merge in to R:
! Copyright (C) 2000, The R Core Team
! R has lgammafn, and lgamma is not part of ISO C
!
!
! stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
!             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
!             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
!
! see also lgammacor() in ./lgammacor.c  which computes almost the same!
!


FUNCTION stirlerr(n) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: n

  DOUBLE PRECISION, PARAMETER :: M_LN_SQRT_2PI = 0.918938533204672741780329736406d0

  DOUBLE PRECISION, PARAMETER :: S0 = 0.083333333333333333333d0         ! 1/12
  DOUBLE PRECISION, PARAMETER :: S1 = 0.00277777777777777777778d0       ! 1/360
  DOUBLE PRECISION, PARAMETER :: S2 = 0.00079365079365079365079365d0    ! 1/1260
  DOUBLE PRECISION, PARAMETER :: S3 = 0.000595238095238095238095238d0   ! 1/1680
  DOUBLE PRECISION, PARAMETER :: S4 = 0.0008417508417508417508417508d0  ! 1/1188

! error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.

  DOUBLE PRECISION, PARAMETER :: sferr_halves(31) = (/ &
	0.0d0,                           & ! n=0 - wrong, place holder only
	0.1534264097200273452913848d0,   & ! 0.5
	0.0810614667953272582196702d0,   & ! 1.0
	0.0548141210519176538961390d0,   & ! 1.5
	0.0413406959554092940938221d0,   & ! 2.0
	0.03316287351993628748511048d0,  & ! 2.5
	0.02767792568499833914878929d0,  & ! 3.0
	0.02374616365629749597132920d0,  & ! 3.5
	0.02079067210376509311152277d0,  & ! 4.0
	0.01848845053267318523077934d0,  & ! 4.5
	0.01664469118982119216319487d0,  & ! 5.0
	0.01513497322191737887351255d0,  & ! 5.5
	0.01387612882307074799874573d0,  & ! 6.0
	0.01281046524292022692424986d0,  & ! 6.5
	0.01189670994589177009505572d0,  & ! 7.0
	0.01110455975820691732662991d0,  & ! 7.5
	0.010411265261972096497478567d0, & ! 8.0
	0.009799416126158803298389475d0, & ! 8.5
	0.009255462182712732917728637d0, & ! 9.0
	0.008768700134139385462952823d0, & ! 9.5
	0.008330563433362871256469318d0, & ! 10.0
	0.007934114564314020547248100d0, & ! 10.5
	0.007573675487951840794972024d0, & ! 11.0
	0.007244554301320383179543912d0, & ! 11.5
	0.006942840107209529865664152d0, & ! 12.0
	0.006665247032707682442354394d0, & ! 12.5
	0.006408994188004207068439631d0, & ! 13.0
	0.006171712263039457647532867d0, & ! 13.5
	0.005951370112758847735624416d0, & ! 14.0
	0.005746216513010115682023589d0, & ! 14.5
	0.005554733551962801371038690d0  & ! 15.0
  /)

  DOUBLE PRECISION nn, fn_val
  DOUBLE PRECISION lgammafn

  IF (n .le. 15.0d0) THEN
     nn = n + n
     IF (nn .eq. AINT(nn)) THEN
        fn_val = sferr_halves(INT(nn)+1)
     ELSE
        fn_val = lgammafn(n + 1.0d0) - (n + 0.5d0)*log(n) + n - M_LN_SQRT_2PI
     END IF
     RETURN
  END IF

  nn = n*n
  
  IF (n .gt. 500.0d0) THEN
     fn_val = (S0-S1/nn) / n
     RETURN
  END IF

  IF (n .gt. 80.0d0) THEN
     fn_val = (S0-(S1-S2/nn)/nn) / n
     RETURN
  END IF

  IF (n .gt. 35.0d0) THEN
     fn_val = (S0-(S1-(S2-S3/nn)/nn)/nn) / n
     RETURN
  END IF

  ! else:  15 < n <= 35 :
  fn_val = (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn) / n
  RETURN

END FUNCTION stirlerr




!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!  SYNOPSIS
!
!    int chebyshev_init(double *dos, int nos, double eta)
!    double chebyshev_eval(double x, double *a, int n)
!
!  DESCRIPTION
!
!    "chebyshev_init" determines the number of terms for the
!    double precision orthogonal series "dos" needed to insure
!    the error is no larger than "eta".  Ordinarily eta will be
!    chosen to be one-tenth machine precision.
!
!    "chebyshev_eval" evaluates the n-term Chebyshev series
!    "a" at "x".
!
!  NOTES
!
!    These routines are translations into C of Fortran routines
!    by W. Fullerton of Los Alamos Scientific Laboratory.
!
!    Based on the Fortran routine dcsevl by W. Fullerton.
!    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).



FUNCTION chebyshev_init (dos, nos, eta) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nos

  DOUBLE PRECISION, INTENT(IN) :: eta, dos(nos)

  INTEGER i, ii, fn_val

  DOUBLE PRECISION err

  IF (nos .lt. 1) THEN
     fn_val = 0
     RETURN
  END IF

  err = 0.0d0
  i = 0			            ! just to avoid compiler warnings
  DO ii = 1, nos
     i = nos - ii
     err = err + abs(dos(i+1))
     IF (err .gt. eta) THEN
        fn_val = i
        EXIT
     END IF
  END DO

  fn_val = i
  RETURN

END FUNCTION chebyshev_init



FUNCTION chebyshev_eval (x, a, n) RESULT(fn_val)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n

  DOUBLE PRECISION, INTENT(IN) :: x, a(n)

  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION :: b0, b1, b2, twox, fn_val

  INTEGER :: i

  IF (n .lt. 1 .or. n .gt. 1000) THEN
     PRINT *, ' Error! Unable to evaluate the the n-term Chebyshev series!'
     fn_val = NaN
     RETURN
  END IF

  IF (x .lt. -1.1 .or. x .gt. 1.1) THEN
     PRINT *, ' Error! Unable to evaluate the the n-term Chebyshev series!'
     fn_val = NaN
     RETURN
  END IF

  twox = x * 2.0d0
  b2 = 0.0d0;  b1 = 0.0d0;  b0 = 0.0d0

  DO i = 1, n
     b2 = b1
     b1 = b0
     b0 = twox * b1 - b2 + a(n-i+1)
  END DO

  fn_val = (b0 - b2) * 0.5d0
  RETURN

END FUNCTION chebyshev_eval




!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000-2001 The R Core Team
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!  SYNOPSIS
!
!    #include <Rmath.h>
!    double lgammacor(double x);
!
!  DESCRIPTION
!
!    Compute the log gamma correction factor for x >= 10 so that
!
!    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
!
!    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
!
!  NOTES
!
!    This routine is a translation into C of a Fortran subroutine
!    written by W. Fullerton of Los Alamos Scientific Laboratory.
!
!  SEE ALSO
!
!    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
!    is faster and cleaner, but is only defined "fast" for half integers.


FUNCTION lgammacor (x) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  DOUBLE PRECISION, PARAMETER :: xbig = 94906265.62425156d0
  DOUBLE PRECISION, PARAMETER :: xmax = 3.745194030963158d306
  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION, PARAMETER :: algmcs(15) = (/ &
	+0.1666389480451863247205729650822d+0, &
	-0.1384948176067563840732986059135d-4, &
	+0.9810825646924729426157171547487d-8, &
	-0.1809129475572494194263306266719d-10, &
	+0.6221098041892605227126015543416d-13, &
	-0.3399615005417721944303330599666d-15, &
	+0.2683181998482698748957538846666d-17, &
	-0.2868042435334643284144622399999d-19, &
	+0.3962837061046434803679306666666d-21, &
	-0.6831888753985766870111999999999d-23, &
	+0.1429227355942498147573333333333d-24, &
	-0.3547598158101070547199999999999d-26, &
	+0.1025680058010470912000000000000d-27, &
	-0.3401102254316748799999999999999d-29, &
	+0.1276642195630062933333333333333d-30 &
  /)

  DOUBLE PRECISION tmp, chebyshev_eval, fn_val

  IF (x .lt. 10.0d0) THEN
     PRINT *, ' Error! Unable to evaluate the gamma correction term'
     fn_val = NaN
     RETURN
  ELSE IF (x .gt. xmax) THEN
     PRINT *, ' Warning! Underflow in "lgammacor"!'
  ELSE IF (x .lt. xbig) THEN
     tmp = 10 / x
     fn_val = chebyshev_eval (tmp * tmp * 2.0d0 - 1.0d0, algmcs, 5) / x
     RETURN
  END IF

  fn_val = 1.0d0 / (x * 12.0d0)
  RETURN

END FUNCTION lgammacor



!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000-2012 The R Core Team
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!  SYNOPSIS
!
!    #include <Rmath.h>
!    double lgammafn_sign(double x, int *sgn);
!    double lgammafn(double x);
!
!  DESCRIPTION
!
!    The function lgammafn computes log|gamma(x)|.  The function
!    lgammafn_sign in addition assigns the sign of the gamma function
!    to the address in the second argument if this is not NULL.
!
!  NOTES
!
!    This routine is a translation into C of a Fortran subroutine
!    by W. Fullerton of Los Alamos Scientific Laboratory.
!
!    The accuracy of this routine compares (very) favourably
!    with those of the Sun Microsystems portable mathematical
!    library.

FUNCTION lgammafn (x) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  DOUBLE PRECISION, PARAMETER :: xmax = 2.5327372760800758d+305
  DOUBLE PRECISION, PARAMETER :: dxrel = 1.490116119384765625d-8

  DOUBLE PRECISION, PARAMETER :: M_PI = 3.141592653589793238462643383280d0
  DOUBLE PRECISION, PARAMETER :: M_LN_SQRT_2PI = 0.918938533204672741780329736406d0
  DOUBLE PRECISION, PARAMETER :: M_LN_SQRT_PId2 = 0.225791352644727432363097614947d0

  DOUBLE PRECISION, PARAMETER :: POS_INF = 9999.0d99

  DOUBLE PRECISION :: ans, y, sinpiy, gammafn, lgammacor, fn_val


  IF (x .lt. 0.0d0 .and. x .eq. AINT(x)) THEN  ! Negative integer argument
     PRINT *, ' Warning! Gamma function evaluated at negative integer! '
     fn_val = POS_INF
     RETURN              ! +Inf, since lgamma(x) = log|gamma(x)|
  END IF

  y = abs(x)

  IF (y .lt. 1.0d-306) THEN
     fn_val = -log(y)         ! denormalized range, R change
     RETURN
  END IF

  IF (y .le. 10.0d0) THEN
     fn_val = log(abs(gammafn(x)))
     RETURN
  END IF

  ! ELSE  y = |x| > 10 ----------------------

  IF (y .gt. xmax) THEN
     PRINT *, ' Warning! Too large argument for gamma function! '
     fn_val = POS_INF
     RETURN
  END IF

  IF (x .gt. 0.0d0) THEN       ! i.e. y = x > 10
     IF (x .gt. 1.0d17) THEN
        fn_val = (x*(log(x) - 1.0d0))
     ELSE IF (x .gt. 4934720.0d0) THEN
        fn_val = M_LN_SQRT_2PI + (x - 0.5d0) * log(x) - x
     ELSE
        fn_val = M_LN_SQRT_2PI + (x - 0.5d0) * log(x) - x + lgammacor(x)
     END IF
     RETURN
  END IF

                                ! else: x < -10; y = -x
  sinpiy = abs(sin(y*M_PI))

  ans = M_LN_SQRT_PId2 + (x - 0.5d0) * log(y) - x - log(sinpiy) - lgammacor(y)

  IF (abs((x - AINT(x - 0.5d0)) * ans / x) < dxrel) THEN
     PRINT *, ' Warning! Argument of gamma function is too near a negative integer, '
     PRINT *, '  precision is less than half of usual precision! '
  END IF

  fn_val = ans
  RETURN

END FUNCTION lgammafn

