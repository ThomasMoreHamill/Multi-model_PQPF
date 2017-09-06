!-----------------------------------------------------------------------------
!
! Evaluation of the CDF of a gamma distribution
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!   Mathlib : A C Library of Special Functions
!   Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
!   Copyright (C) 2005-10 The R Foundation
!   Copyright (C) 2006-10 The R Core Team
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
!	#include <Rmath.h>
!
!	double pgamma (double x, double alph, double scale,
!		       int lower_tail, int log_p)
!
!	double log1pmx	(double x)
!	double lgamma1p (double a)
!
!	double logspace_add (double logx, double logy)
!	double logspace_sub (double logx, double logy)
!
!
!  DESCRIPTION
!
!	This function computes the distribution function for the
!	gamma distribution with shape parameter alph and scale parameter
!	scale.	This is also known as the incomplete gamma function.
!	See Abramowitz and Stegun (6.5.1) for example.
!
!  NOTES
!
!	Complete redesign by Morten Welinder, originally for Gnumeric.
!	Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
!



! Continued fraction for calculation of
!    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
!
! auxilary in log1pmx() and lgamma1p()

FUNCTION logcf (x, i, d, eps) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, i, d, eps 

  ! Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77
  DOUBLE PRECISION, PARAMETER :: scalefactor = 1.157921d+77

  DOUBLE PRECISION :: c1, c2, c3, c4, a1, a2, b1, b2, fn_val

  c1 = 2.0d0 * d
  c2 = i + d
  c4 = c2 + d
  a1 = c2
  b1 = i * (c2 - i * x)
  b2 = d * d * x
  a2 = c4 * c2 - b2

  IF (i.le.0.0d0 .or. d.lt.0.0d0) THEN
     PRINT *, 'Error in function logcf! Program terminated.'
     CALL EXIT
  END IF

  b2 = c4 * b1 - i * b2

  DO WHILE ( abs(a2 * b1 - a1 * b2) .gt. abs(eps * b1 * b2) )

     c3 = c2*c2*x
     c2 = c2 + d
     c4 = c4 + d
     a1 = c4 * a2 - c3 * a1
     b1 = c4 * b2 - c3 * b1

     c3 = c1 * c1 * x
     c1 = c1 + d
     c4 = c4 + d
     a2 = c4 * a1 - c3 * a2
     b2 = c4 * b1 - c3 * b2

     IF (abs(b2) .gt. scalefactor) THEN
        a1 = a1 / scalefactor
        b1 = b1 / scalefactor
        a2 = a2 / scalefactor
	b2 = b2 / scalefactor
     ELSE IF (abs(b2) .lt. 1.0d0/scalefactor) THEN
        a1 = a1 * scalefactor
        b1 = b1 * scalefactor
        a2 = a2 * scalefactor
        b2 = b2 * scalefactor
     END IF

  END DO

  fn_val = a2 / b2
  RETURN

END FUNCTION logcf



! Accurate calculation of log(1+x)-x, particularly for small x

FUNCTION log1pmx (x) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  DOUBLE PRECISION, PARAMETER :: minLog1Value = -0.79149064
  DOUBLE PRECISION, PARAMETER :: tol_logcf = 1.0d-14
  DOUBLE PRECISION, PARAMETER :: two = 2.0d0

  DOUBLE PRECISION :: r, y, logcf, log1p, fn_val

  IF (x.gt.1.0d0 .or. x.lt.minLog1Value) THEN
     fn_val = log1p(x) - x
  ELSE
     ! -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
     ! log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
     ! ---------------------------------------------
     ! S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)

     r = x / (2.0d0 + x)
     y = r * r
     IF (abs(x) .lt. 1.0d-2) THEN
        fn_val = r * ((((two / 9.0d0 * y + two / 7.0d0) * y + &
             two / 5.0d0) * y + two / 3.0d0) * y - x)
     ELSE
        fn_val = r * (2.0d0 * y * logcf (y, 3.0d0, 2.0d0, tol_logcf) - x)
     END IF
  END IF

  RETURN

END FUNCTION log1pmx



! Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5)

FUNCTION lgamma1p (a) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: a

  INTEGER, PARAMETER :: N = 40

  DOUBLE PRECISION, PARAMETER :: eulers_const = 0.5772156649015328606065120900824024d0

! coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 :

  DOUBLE PRECISION, PARAMETER :: coeffs(40) = (/ &
	0.3224670334241132182362075833230126d-0, &  ! = (zeta(2)-1)/2
	0.6735230105319809513324605383715000d-1, &  ! = (zeta(3)-1)/3
	0.2058080842778454787900092413529198d-1, &
	0.7385551028673985266273097291406834d-2, &
	0.2890510330741523285752988298486755d-2, &
	0.1192753911703260977113935692828109d-2, &
	0.5096695247430424223356548135815582d-3, &
	0.2231547584535793797614188036013401d-3, &
	0.9945751278180853371459589003190170d-4, &
	0.4492623673813314170020750240635786d-4, &
	0.2050721277567069155316650397830591d-4, &
	0.9439488275268395903987425104415055d-5, &
	0.4374866789907487804181793223952411d-5, &
	0.2039215753801366236781900709670839d-5, &
	0.9551412130407419832857179772951265d-6, &
	0.4492469198764566043294290331193655d-6, &
	0.2120718480555466586923135901077628d-6, &
	0.1004322482396809960872083050053344d-6, &
	0.4769810169363980565760193417246730d-7, &
	0.2271109460894316491031998116062124d-7, &
	0.1083865921489695409107491757968159d-7, &
	0.5183475041970046655121248647057669d-8, &
	0.2483674543802478317185008663991718d-8, &
	0.1192140140586091207442548202774640d-8, &
	0.5731367241678862013330194857961011d-9, &
	0.2759522885124233145178149692816341d-9, &
	0.1330476437424448948149715720858008d-9, &
	0.6422964563838100022082448087644648d-10, &
	0.3104424774732227276239215783404066d-10, &
	0.1502138408075414217093301048780668d-10, &
	0.7275974480239079662504549924814047d-11, &
	0.3527742476575915083615072228655483d-11, &
	0.1711991790559617908601084114443031d-11, &
	0.8315385841420284819798357793954418d-12, &
	0.4042200525289440065536008957032895d-12, &
	0.1966475631096616490411045679010286d-12, &
	0.9573630387838555763782200936508615d-13, &
	0.4664076026428374224576492565974577d-13, &
	0.2273736960065972320633279596737272d-13, &
	0.1109139947083452201658320007192334d-13 &  ! = (zeta(40+1)-1)/(40+1)
    /)

    DOUBLE PRECISION, PARAMETER :: c = 0.2273736845824652515226821577978691d-12  ! zeta(N+2)-1
    DOUBLE PRECISION, PARAMETER :: tol_logcf = 1d-14;
    DOUBLE PRECISION :: lgam, lgammafn, logcf, log1pmx, fn_val

    INTEGER :: i

    IF (abs(a) .ge. 0.5) THEN
       fn_val = lgammafn (a + 1.0d0)

    ! Abramowitz & Stegun 6.1.33 : for |x| < 2,
    !  <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
    !  where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
    !
    !  Here, another convergence acceleration trick is used to compute
    !  lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n

    ELSE
       lgam = c * logcf(-a/2.0d0, DBLE(N+2), 1.0d0, tol_logcf)
       DO i = N, 1, -1
          lgam = coeffs(i) - a * lgam
       END DO
       fn_val = (a * lgam - eulers_const) * a - log1pmx(a)
    END IF

    RETURN

END FUNCTION lgamma1p



! Compute the log of a sum from logs of terms, i.e.,
!
!     log (exp (logx) + exp (logy))
!
! without causing overflows and without throwing away large handfuls
! of accuracy.

FUNCTION logspace_add (logx, logy) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: logx, logy
  DOUBLE PRECISION :: log1p, fn_val

  fn_val = max(logx,logy) + log1p (exp(-abs(logx-logy)))
  RETURN

END FUNCTION logspace_add



! Compute the log of a difference from logs of terms, i.e.,
!
!     log (exp (logx) - exp (logy))
!
! without causing overflows and without throwing away large handfuls
! of accuracy.

FUNCTION logspace_sub (logx, logy) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: logx, logy

  DOUBLE PRECISION, PARAMETER :: M_LN2 = 0.693147180559945309417232121458d0

  DOUBLE PRECISION :: logy_minus_logx, log1p, expm1, fn_val

  logy_minus_logx = logy - logx

  IF (logy_minus_logx .gt. -M_LN2) THEN
     fn_val = logx + log(-expm1(logy_minus_logx))
  ELSE
     fn_val = logx + log1p(-exp(logy_minus_logx))
  END IF
  RETURN

END FUNCTION logspace_sub



! dpois_wrap (x_P_1,  lambda, g_log) ==
!   dpois (x_P_1 - 1, lambda, g_log) :=  exp(-L)  L^k / gamma(k+1) ,  k := x_P_1 - 1

FUNCTION dpois_wrap (x_plus_1, lambda, give_log) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x_plus_1, lambda

  LOGICAL, INTENT(IN) :: give_log

  DOUBLE PRECISION, PARAMETER :: M_cutoff = 3.196577d18

  DOUBLE PRECISION :: d, dpois_raw, lgammafn, fn_val

  IF (x_plus_1 .gt. 1.0d0) THEN
     fn_val = dpois_raw (x_plus_1-1.0d0, lambda, give_log)
     RETURN
  END IF

  IF (lambda .gt. abs(x_plus_1-1.0d0)*M_cutoff) THEN
     IF (give_log) THEN
        fn_val = -lambda - lgammafn(x_plus_1)
     ELSE
        fn_val = exp(-lambda - lgammafn(x_plus_1))
     END IF
     RETURN
  END IF

  d = dpois_raw (x_plus_1, lambda, give_log)

  IF (give_log) THEN
     fn_val = d + log(x_plus_1/lambda)
  ELSE
     fn_val = d * (x_plus_1/lambda)
  END IF

  RETURN

END FUNCTION dpois_wrap



! Abramowitz and Stegun 6.5.29 [right]

FUNCTION pgamma_smallx (x, alph, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, alph

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16
  DOUBLE PRECISION, PARAMETER :: M_LN2 = 0.693147180559945309417232121458d0

  DOUBLE PRECISION :: sum, c, n, term, f1, f2, lf2, f1m1, f2m1, arg, fn_val

  DOUBLE PRECISION :: dpois_raw, lgamma1p, log1p, expm1

  sum = 0.0d0
  c = alph
  n = 0.0d0
  term = 1.0d0

    ! Relative to 6.5.29 all terms have been multiplied by alph
    ! and the first, thus being 1, is omitted.

  DO WHILE ( abs(term) .gt. abs(sum)*DBL_EPSILON )
     n = n + 1.0d0
     c = c * (-x/n)
     term = c / (alph + n)
     sum = sum + term
  END DO

  IF (lower_tail) THEN

     IF (log_p) THEN
        f1 = log1p (sum)
     ELSE
        f1 = 1.0d0 + sum
     END IF

     IF (alph .gt. 1.0d0) THEN
        f2 = dpois_raw (alph, x, log_p)
        IF (log_p) THEN
           f2 = f2 + x
        ELSE
           f2 = f2 * exp(x)
        END IF
     ELSE IF (log_p) THEN
        f2 = alph * log(x) - lgamma1p(alph)
     ELSE
        f2 = (x**alph) / exp(lgamma1p(alph))
     END IF

     IF (log_p) THEN
        fn_val = f1 + f2
     ELSE
        fn_val = f1 * f2
     END IF

  ELSE
     lf2 = alph * log(x) - lgamma1p(alph)

     IF (log_p) THEN
        arg = log1p(sum) + lf2
        IF (arg .gt. -M_LN2) THEN
           fn_val = log(-expm1(arg))
        ELSE
           fn_val = log1p(-exp(arg))
        END IF
     ELSE
        f1m1 = sum
        f2m1 = expm1(lf2)
        fn_val = -(f1m1 + f2m1 + f1m1 * f2m1)
     END IF

  END IF
  RETURN

END FUNCTION pgamma_smallx



FUNCTION pd_upper_series (x, yy, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, yy

  LOGICAL, INTENT(IN) :: log_p

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16

  DOUBLE PRECISION :: y, term, sum, fn_val

  y = yy
  term = x / y
  sum = term

  DO
     y = y + 1.0d0
     term = term * (x/y)
     sum = sum + term
     IF (term .le. sum * DBL_EPSILON) EXIT
  END DO

    ! sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
    !	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
    !	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
    !	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}

  IF (log_p) THEN
     fn_val = log(sum)
  ELSE
     fn_val = sum
  END IF

  RETURN

END FUNCTION pd_upper_series



! Continued fraction for calculation of scaled upper-tail F_{gamma}
!  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]

FUNCTION pd_lower_cf (y, d) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: y, d

  ! Scalefactor:= (2^32)^8
  INTEGER, PARAMETER :: max_it = 200000
  DOUBLE PRECISION, PARAMETER :: scalefactor = 1.157921d+77
  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16

  INTEGER :: i
  DOUBLE PRECISION :: f, of, f0, c2, c3, c4,  a1, b1,  a2, b2, fn_val

  f = 0.0d0

  IF (y .eq. 0.0d0 ) THEN
     fn_val = 0.0d0
     RETURN
  END IF

  f0 = y/d
  ! Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE):
  IF ( abs(y-1.0d0) .lt. abs(d)*DBL_EPSILON) THEN
     fn_val = f0
     RETURN
  END IF

  IF (f0 .gt. 1.0d0)  f0 = 1.0d0
  c2 = y
  c4 = d ! original (y,d), *not* potentially scaled ones!
  
  a1 = 0.0d0;  b1 = 1.0d0
  a2 = y;      b2 = d
  
  DO WHILE (b2 .gt. scalefactor)
     a1 = a1 / scalefactor
     b1 = b1 / scalefactor
     a2 = a2 / scalefactor
     b2 = b2 / scalefactor
  END DO

  i = 0;  of = -1.0d0  ! far away

  DO WHILE (i .lt. max_it)

     i = i + 1;  c2 = c2 - 1.0d0;   c3 = DBLE(i) * c2;   c4 = c4 + 2.0d0
     ! c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd
     a1 = c4 * a2 + c3 * a1
     b1 = c4 * b2 + c3 * b1

     i = i + 1;   c2 = c2 - 1.0d0;   c3 = DBLE(i) * c2;   c4 = c4 + 2.0d0;
     ! c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even
     a2 = c4 * a1 + c3 * a2
     b2 = c4 * b1 + c3 * b2

     IF (b2 .gt. scalefactor) THEN
        a1 = a1 / scalefactor
        b1 = b1 / scalefactor
        a2 = a2 / scalefactor
        b2 = b2 / scalefactor
     END IF

     IF (b2 .ne. 0.0d0) THEN
        f = a2 / b2
        ! convergence check: relative; "absolute" for very small f :
        IF ( abs(f-of) .le. DBL_EPSILON * max(f0, abs(f))) THEN
           fn_val = f
           RETURN
        END IF
        of = f
     END IF

  END DO

  PRINT *, ' ** NON-convergence in pgamma: pd_lower_cf ** '
  fn_val = f  ! should not happen ...
  RETURN

END FUNCTION pd_lower_cf



FUNCTION pd_lower_series (lambda, yy) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) ::lambda, yy

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16

  DOUBLE PRECISION :: y, term, sum, f, pd_lower_cf, fn_val

  y = yy
  term = 1.0d0
  sum  = 0.0d0

  DO WHILE (y.ge.1.0d0 .and. term.gt.sum*DBL_EPSILON)
     term = term * (y/lambda)
     sum = sum + term
     y = y - 1.0d0
  END DO

    ! sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
    !	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
    !	   ~  y/lambda + o(y/lambda)
    
  IF (y .ne. floor(y)) THEN
     ! The series does not converge as the terms start getting
     ! bigger (besides flipping sign) for y < -lambda.

     ! FIXME: in quite few cases, adding  term*f  has no effect (f too small)
     !	  and is unnecessary e.g. for pgamma(4e12, 121.1)

     f = pd_lower_cf (y, lambda + 1.0d0 - y)
     sum = sum + term * f
  END IF

  fn_val = sum
  RETURN

END FUNCTION pd_lower_series



! Compute the following ratio with higher accuracy that would be had
! from doing it directly.
!
!		 dnorm (x, 0, 1, FALSE)
!	   ----------------------------------
!	   pnorm (x, 0, 1, lower_tail, FALSE)
!
! Abramowitz & Stegun 26.2.12

FUNCTION dpnorm (xx, l_tail, lp) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: xx, lp

  LOGICAL, INTENT(IN) :: l_tail

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16

  DOUBLE PRECISION :: x, term, sum, x2, i, d, dnorm, fn_val

  LOGICAL :: lower_tail

  ! So as not to repeat a pnorm call, we expect
  !
  !	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
  !
  ! but use it only in the non-critical case where either x is small
  ! or p==exp(lp) is close to 1.

  x = xx
  lower_tail = l_tail

  IF (x .lt. 0.0d0) THEN
     x = -x
     lower_tail = .not.lower_tail
  END IF

  IF (x.gt.10.0d0 .and. .not.lower_tail) THEN
     term = 1.0d0 / x
     sum = term
     x2 = x * x
     i = 1.0d0
     DO
        term = term * (-i/x2)
        sum = sum + term
        i = i + 2.0d0
        IF ( abs(term) .le. DBL_EPSILON*sum) EXIT
     END DO
     fn_val = 1.0d0 / sum
  ELSE
     d = dnorm (x, 0.0d0, 1.0d0, .false.)
     fn_val = d / exp(lp)
  END IF

  RETURN

END FUNCTION dpnorm



! Asymptotic expansion to calculate the probability that Poisson variate
! has value <= x.
! Various assertions about this are made (without proof) at
! http://members.aol.com/iandjmsmith/PoissonApprox.htm

FUNCTION ppois_asymp (x, lambda, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, lambda

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: coefs_a(7) = (/ &
	2.0d0/3.0d0, &
	-4.0d0/135.0d0, &
	8.0d0/2835.0d0, &
	16.0d0/8505.0d0, &
	-8992.0d0/12629925.0d0, &
	-334144.0d0/492567075.0d0, &
	698752.0d0/1477701225.0d0 &
        /)

  DOUBLE PRECISION, PARAMETER :: coefs_b(7) = (/ &
	1.0d0/12.0d0, &
	1.0d0/288.0d0, &
	-139.0d0/51840.0d0, &
	-571.0d0/2488320.0d0, &
	163879.0d0/209018880.0d0, &
	5246819.0d0/75246796800.0d0, &
	-534703531.0d0/902961561600.0d0 &
        /)

  DOUBLE PRECISION :: elfb, elfb_term
  DOUBLE PRECISION :: res12, res1_term, res1_ig, res2_term, res2_ig
  DOUBLE PRECISION :: dfm, pt_, s2pt, f, np, nd, n_d_over_p, fn_val

  DOUBLE PRECISION :: log1pmx, pnorm, dnorm, dpnorm, log1p

  INTEGER :: i

  dfm = lambda - x
  ! If lambda is large, the distribution is highly concentrated
  !  about lambda.  So representation error in x or lambda can lead
  !  to arbitrarily large values of pt_ and hence divergence of the
  !  coefficients of this approximation.

  pt_ = - log1pmx (dfm / x)
  s2pt = sqrt(2.0d0 * x * pt_)
  IF (dfm .lt. 0.0d0) s2pt = -s2pt

  res12 = 0.0d0
  res1_term = sqrt(x)
  res1_ig = res1_term
  res2_term = s2pt
  res2_ig = res2_term

  DO i = 1, 7
     res12 = res12 + res1_ig * coefs_a(i)
     res12 = res12 + res2_ig * coefs_b(i)
     res1_term = res1_term * (pt_ / i)
     res2_term = res2_term * ((2.0d0 * pt_) / DBLE(2*i + 1))
     res1_ig = (res1_ig / x) + res1_term
     res2_ig = (res2_ig / x) + res2_term
  END DO

  elfb = x
  elfb_term = 1.0d0

  DO i = 1, 7
     elfb = elfb + elfb_term * coefs_b(i)
     elfb_term = elfb_term / x
  END DO

  IF (.not.lower_tail) elfb = -elfb

  f = res12 / elfb

  np = pnorm(s2pt, 0.0d0, 1.0d0, .not.lower_tail, log_p)

  IF (log_p) THEN
     n_d_over_p = dpnorm(s2pt, .not.lower_tail, np)
     fn_val = np + log1p(f * n_d_over_p)
  ELSE
     nd = dnorm(s2pt, 0.0d0, 1.0d0, log_p)
     fn_val = np + f * nd
  END IF

  RETURN

END FUNCTION ppois_asymp




RECURSIVE FUNCTION pgamma_raw (x, alph, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, alph

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16
  DOUBLE PRECISION, PARAMETER :: DBL_MIN = 2.2250738585072014d-308
  DOUBLE PRECISION, PARAMETER :: M_LN2 = 0.693147180559945309417232121458d0
  DOUBLE PRECISION, PARAMETER :: POS_INF = 9999.0d99
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -9999.0d99

  DOUBLE PRECISION :: sum, d, f, fn_val

  DOUBLE PRECISION :: ppois_asymp, pgamma_smallx, pd_upper_series, &
       pd_lower_series, pd_lower_cf, dpois_wrap, expm1, log1p

  IF (x .le. 0.0d0) THEN
     IF(log_p .and. lower_tail)            fn_val = NEG_INF
     IF(log_p .and. .not.lower_tail)       fn_val = 0.0d0
     IF(.not.log_p .and. lower_tail)       fn_val = 0.0d0
     IF(.not.log_p .and. .not.lower_tail)  fn_val = 1.0d0
     RETURN
  END IF

  IF (x .ge. POS_INF) THEN
     IF(log_p .and. lower_tail)            fn_val = 0.0d0
     IF(log_p .and. .not.lower_tail)       fn_val = NEG_INF
     IF(.not.log_p .and. lower_tail)       fn_val = 1.0d0
     IF(.not.log_p .and. .not.lower_tail)  fn_val = 0.0d0
     RETURN
  END IF

  IF (x .lt. 1.0d0) THEN
     fn_val = pgamma_smallx (x, alph, lower_tail, log_p)

  ELSE IF (x .le. alph-1.0d0 .and. x .lt. 0.8d0*(alph+50.0d0)) THEN
     ! incl. large alph compared to x
     sum = pd_upper_series (x, alph, log_p)  ! = x/alph + o(x/alph)
     d = dpois_wrap (alph, x, log_p)

     IF (.not.lower_tail) THEN
        IF (log_p) THEN
           IF (d+sum .gt. -M_LN2) THEN
              fn_val = log(-expm1(d+sum))
           ELSE
              fn_val = log1p(-exp(d+sum))
           END IF
        ELSE
           fn_val = 1.0d0 - d * sum
        END IF
     ELSE
        IF (log_p) THEN
           fn_val = d + sum
        ELSE
           fn_val = d * sum
        END IF
     END IF

  ELSE IF (alph-1.0d0 .lt. x .and. alph .lt. 0.8d0*(x+50.0d0)) THEN
     ! incl. large x compared to alph
     d = dpois_wrap (alph, x, log_p)

     IF (alph .lt. 1.0d0) THEN
        IF (x * DBL_EPSILON .gt. 1.0d0 - alph) THEN
           IF (log_p) THEN
              sum = 0.0d0
           ELSE
              sum = 1.0d0
           END IF
        ELSE
           f = pd_lower_cf (alph, x - (alph - 1.0d0)) * x / alph
           ! = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1)
           IF (log_p) THEN
              sum = log(f)
           ELSE
              sum = f
           END IF
        END IF
     ELSE
        sum = pd_lower_series (x, alph - 1.0d0)  ! = (alph-1)/x + o((alph-1)/x)
        IF (log_p) THEN
           sum = log1p(sum)
        ELSE
           sum = 1.0d0 + sum
        END IF
     END IF

     IF (.not.lower_tail) THEN
        IF (log_p) THEN
           fn_val = d + sum
        ELSE
           fn_val = d * sum
        END IF
     ELSE
        IF (log_p) THEN
           IF (d+sum .gt. -M_LN2) THEN
              fn_val = log(-expm1(d+sum))
           ELSE
              fn_val = log1p(-exp(d+sum))
           END IF
        ELSE
           fn_val = 1 - d * sum
        END IF
     END IF

  ELSE  ! x >= 1 and x fairly near alph.

     fn_val = ppois_asymp (alph-1.0d0, x, (.not.lower_tail), log_p)

  END IF

  ! We lose a fair amount of accuracy to underflow in the cases
  ! where the final result is very close to DBL_MIN. In those
  ! cases, simply redo via log space.

  IF (.not.log_p .and. fn_val .lt. DBL_MIN / DBL_EPSILON) THEN
     fn_val = exp(pgamma_raw (x, alph, lower_tail, .true.))
  END IF

  RETURN

END FUNCTION pgamma_raw



FUNCTION pgamma (xx, alph, scale, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: xx, alph, scale

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: NEG_INF = -9999.0d99
  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION :: x, pgamma_raw, fn_val

  x = xx

  IF (alph .lt. 0.0d0 .or. scale .le. 0.0d0) THEN
     fn_val = NaN
     RETURN
  END IF

  x = x / scale

  IF (alph .eq. 0.0d0) THEN  ! limit case; useful e.g. in pnchisq()
     IF (x .le. 0.0d0) THEN
        IF(log_p .and. lower_tail)            fn_val = NEG_INF
        IF(log_p .and. .not.lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 1.0d0
     ELSE
        IF(log_p .and. lower_tail)            fn_val = 0.0d0
        IF(log_p .and. .not.lower_tail)       fn_val = NEG_INF
        IF(.not.log_p .and. lower_tail)       fn_val = 1.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 0.0d0
     END IF               ! <= assert  pgamma(0,0) ==> 0
     RETURN
  END IF

  fn_val = pgamma_raw (x, alph, lower_tail, log_p)
  RETURN

END FUNCTION pgamma



FUNCTION dpois_raw (x, lambda, give_log) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, lambda

  LOGICAL, INTENT(IN) :: give_log

  DOUBLE PRECISION, PARAMETER :: DBL_MIN = 2.2250738585072014d-308
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -9999.0d99
  DOUBLE PRECISION, PARAMETER :: M_2PI = 6.283185307179586476925286766559d0

  DOUBLE PRECISION :: lgammafn, stirlerr, bd0, fn_val

  ! x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
  ! lambda >= 0

  IF (lambda .eq. 0.0d0) THEN
     IF (x .eq. 0.0d0) THEN
        IF (give_log) THEN
           fn_val = 0.0d0
        ELSE
           fn_val = 1.0d0
        END IF
     ELSE
        IF (give_log) THEN
           fn_val = NEG_INF
        ELSE
           fn_val = 0.0d0
        END IF
     END IF
     RETURN
  END IF

  IF (x .lt. 0.0d0) THEN
     IF (give_log) THEN
        fn_val = NEG_INF
     ELSE
        fn_val = 0.0d0
     END IF
     RETURN
  END IF

  IF (x .le. lambda * DBL_MIN) THEN
     IF (give_log) THEN
        fn_val = -lambda
     ELSE
        fn_val = exp(-lambda)
     END IF
     RETURN
  END IF 

  IF (lambda .lt. x * DBL_MIN) THEN
     IF (give_log) THEN
        fn_val = -lambda + x*log(lambda) -lgammafn(x+1.0d0)
     ELSE
        fn_val = exp(-lambda + x*log(lambda) -lgammafn(x+1.0d0))
     END IF
     RETURN
  END IF 


  IF (give_log) THEN
     fn_val = -0.5*log(M_2PI*x) -stirlerr(x) - bd0(x,lambda)
  ELSE
     fn_val = exp( -stirlerr(x) - bd0(x,lambda) ) / sqrt(M_2PI*x)
  END IF
  RETURN

END FUNCTION dpois_raw



!  AUTHOR
!	Catherine Loader, catherine@research.bell-labs.com.
!	October 23, 2000.
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
!	Evaluates the "deviance part"
!	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
!		  =  x * log(x/M) + M - x
!	where M = E[X] = n*p (or = lambda), for	  x, M > 0
!
!	in a manner that should be stable (with small relative error)
!	for all x and M=np. In particular for x/np close to 1, direct
!	evaluation fails, and evaluation is based on the Taylor series
!	of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).


FUNCTION bd0 (x, np) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x, np

  DOUBLE PRECISION, PARAMETER :: DBL_MIN = 2.2250738585072014d-308
  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION :: ej, s, s1, v, fn_val

  INTEGER :: j

  IF (np .eq. 0.0d0) THEN  ! || !R_FINITE(x) || !R_FINITE(np)
     PRINT *, ' Error in bd0, returned NaN !'
     fn_val = NaN
     RETURN
  END IF

  IF (abs(x-np) .lt. 0.1d0*(x+np)) THEN
     v = (x-np)/(x+np)    ! might underflow to 0
     s = (x-np)*v         ! s using v -- change by MM

     IF(abs(s) .lt. DBL_MIN) THEN
        fn_val = s
        RETURN
     END IF

     ej = 2.0d0 * x * v
     v = v * v

     DO j = 1, 999          ! Taylor series; 1000: no infinite loop
                            ! as |v| < .1,  v^2000 is "zero"
        ej = ej * v                 ! = v^(2j+1)
        s1 = s + ej/DBLE(ISHFT(j,1)+1)
        IF (s1 .eq. s) THEN         ! last term was effectively 0
           fn_val = s1
           RETURN
        END IF
        s = s1
     END DO
  END IF
            ! else:  | x - np |  is not too small

  fn_val = x * log(x/np) + np - x
  RETURN

END FUNCTION bd0
