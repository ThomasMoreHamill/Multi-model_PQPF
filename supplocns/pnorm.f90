!-----------------------------------------------------------------------------
!
! Evaluation of the CDF of a Gaussian distribution
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998	    Ross Ihaka
!  Copyright (C) 2000-2013 The R Core Team
!  Copyright (C) 2003	    The R Foundation
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
!   #include <Rmath.h>
!
!   double pnorm5(double x, double mu, double sigma, int lower_tail,int log_p);
!	   {pnorm (..) is synonymous and preferred inside R}
!
!   void   pnorm_both(double x, double *cum, double *ccum,
!		       int i_tail, int log_p);
!
!  DESCRIPTION
!
!	The main computation evaluates near-minimax approximations derived
!	from those in "Rational Chebyshev approximations for the error
!	function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
!	transportable program uses rational functions that theoretically
!	approximate the normal distribution function to at least 18
!	significant decimal digits.  The accuracy achieved depends on the
!	arithmetic system, the compiler, the intrinsic functions, and
!	proper selection of the machine-dependent constants.
!
!  REFERENCE
!
!	Cody, W. D. (1993).
!	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
!	Special Function Routines and Test Drivers".
!	ACM Transactions on Mathematical Software. 19, 22-32.
!
!  EXTENSIONS
!
!  The "_both" , lower, upper, and log_p  variants were added by
!  Martin Maechler, Jan.2000;
!  as well as log1p() and similar improvements later on.
!
!  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
!  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
!  with the original Cody code.
!


FUNCTION pnorm (xx, mu, sigma, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: xx, mu, sigma

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: POS_INF = 1.0d99
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -1.0d99

  INTEGER :: i_tail
  DOUBLE PRECISION :: x, p, cp, fn_val

    ! Note: The structure of these checks has been carefully thought through.
    !  For example, if x == mu and sigma == 0, we get the correct answer 1.

  x = xx

  IF (sigma .lt. 0.0d0) THEN
     PRINT *, 'Error in function dnorm! Program terminated.'
     CALL EXIT
  END IF

  IF (sigma .eq. 0.0d0) THEN
     IF (x .lt. mu) THEN
        IF(log_p .and. lower_tail)            fn_val = NEG_INF
        IF(log_p .and. .not.lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 1.0d0
     ELSE
        IF(log_p .and. lower_tail)            fn_val = 0.0d0
        IF(log_p .and. .not.lower_tail)       fn_val = NEG_INF
        IF(.not.log_p .and. lower_tail)       fn_val = 1.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 0.0d0
     END IF
     RETURN
  END IF

  p = (x - mu) / sigma

  IF(p .le. NEG_INF .or. p .ge. POS_INF) THEN
     IF (x .lt. mu) THEN
        IF(log_p .and. lower_tail)            fn_val = NEG_INF
        IF(log_p .and. .not.lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. lower_tail)       fn_val = 0.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 1.0d0
     ELSE
        IF(log_p .and. lower_tail)            fn_val = 0.0d0
        IF(log_p .and. .not.lower_tail)       fn_val = NEG_INF
        IF(.not.log_p .and. lower_tail)       fn_val = 1.0d0
        IF(.not.log_p .and. .not.lower_tail)  fn_val = 0.0d0
     END IF
     RETURN
  END IF

  x = p

  IF (lower_tail) THEN
     i_tail = 0
  ELSE
     i_tail = 1
  END IF

  CALL pnorm_both(x, p, cp, i_tail, log_p)

  IF (lower_tail) THEN
     fn_val = p
  ELSE
     fn_val = cp
  END IF

  RETURN

END FUNCTION pnorm



SUBROUTINE pnorm_both (x, cum, ccum, i_tail, log_p)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x
  DOUBLE PRECISION, INTENT(INOUT) :: cum, ccum

  INTEGER, INTENT(IN) :: i_tail
  LOGICAL, INTENT(IN) :: log_p

  ! i_tail in {0,1,2} means: "lower", "upper", or "both" :
  !   if(lower) return  *cum := P[X <= x]
  !   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16
  DOUBLE PRECISION, PARAMETER :: DBL_MIN = 2.2250738585072014d-308
  DOUBLE PRECISION, PARAMETER :: SIXTEN	= 16.0d0                     ! Cutoff allowing exact "*" and "/"
  DOUBLE PRECISION, PARAMETER :: M_SQRT_32 = 5.6568542494923801952067548968380
  DOUBLE PRECISION, PARAMETER :: M_1_SQRT_2PI =	0.398942280401432677939946059934d0
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -1.0d99


  DOUBLE PRECISION, PARAMETER :: a(5) = (/ &
	2.2352520354606839287d0, &
	161.02823106855587881d0, &
	1067.6894854603709582d0, &
	18154.981253343561249d0, &
	0.065682337918207449113d0 &
  /)


  DOUBLE PRECISION, PARAMETER :: b(4) = (/ &
	47.20258190468824187d0, &
	976.09855173777669322d0, &
	10260.932208618978205d0, &
	45507.789335026729956d0 &
  /)

  DOUBLE PRECISION, PARAMETER :: c(9) = (/ &
	0.39894151208813466764d0, &
	8.8831497943883759412d0, &
	93.506656132177855979d0, &
	597.27027639480026226d0, &
	2494.5375852903726711d0, &
	6848.1904505362823326d0, &
	11602.651437647350124d0, &
	9842.7148383839780218d0, &
	1.0765576773720192317d-8 &
  /)

  DOUBLE PRECISION, PARAMETER :: d(8) = (/ &
	22.266688044328115691d0, &
	235.38790178262499861d0, &
	1519.377599407554805d0, &
	6485.558298266760755d0, &
	18615.571640885098091d0, &
	34900.952721145977266d0, &
	38912.003286093271411d0, &
	19685.429676859990727d0 &
  /)

  DOUBLE PRECISION, PARAMETER :: p(6) = (/ &
	0.21589853405795699d0, &
	0.1274011611602473639d0, &
	0.022235277870649807d0, &
	0.001421619193227893466d0, &
	2.9112874951168792d-5, &
	0.02307344176494017303d0 &
  /)

  DOUBLE PRECISION, PARAMETER :: q(5) = (/ &
	1.28426009614491121d0, &
	0.468238212480865118d0, &
	0.0659881378689285515d0, &
	0.00378239633202758244d0, &
	7.29751555083966205d-5 &
  /)

  INTEGER :: i
  LOGICAL :: lower, upper

  DOUBLE PRECISION :: xden, xnum, temp, del, eps, xsq, y, log1p


  ! Consider changing these :
  eps = DBL_EPSILON * 0.5d0

  ! i_tail in {0,1,2} =^= {lower, upper, both}
  lower = (i_tail .ne. 1)
  upper = (i_tail .eq. 0)

  y = abs(x)
  IF (y .le. 0.67448975d0) THEN   ! qnorm(3/4) = .6744.... -- earlier had 0.66291
     IF (y .gt. eps) THEN
        xsq = x * x
        xnum = a(5) * xsq
        xden = xsq
        DO i = 1, 3
           xnum = (xnum + a(i)) * xsq
           xden = (xden + b(i)) * xsq
        END DO
     ELSE
	xnum = 0.0d0;  xden = 0.0d0
     END IF

     temp = x * (xnum + a(4)) / (xden + b(4))
     IF (lower)  cum = 0.5d0 + temp
     IF (upper) ccum = 0.5d0 - temp
     IF (log_p) THEN
        IF (lower)  cum = log(cum)
        IF (upper) ccum = log(ccum)
     END IF
    
  ELSE IF (y .le. M_SQRT_32) THEN
     ! Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657

     xnum = c(9) * y
     xden = y
     DO i = 1, 7
        xnum = (xnum + c(i)) * y
        xden = (xden + d(i)) * y
     END DO
     temp = (xnum + c(8)) / (xden + d(8))

     xsq = AINT(y * SIXTEN) / SIXTEN
     del = (y - xsq) * (y + xsq)
     IF (log_p) THEN
        cum = (-xsq * xsq * 0.5d0) + (-del * 0.5d0) + log(temp)
        IF ((lower .and. x .gt. 0.0d0) .or. (upper .and. x .lt. 0.0d0)) THEN
           ccum = log1p(-exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp)
        END IF
     ELSE
        cum = exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp
        ccum = 1.0d0 - cum
     END IF

     IF (x .gt. 0.0d0) THEN ! swap  ccum <--> cum
        temp = cum
        IF (lower) cum = ccum
        ccum = temp
     END IF

     ! else |x| > sqrt(32) = 5.657 :
     ! the next two case differentiations were really for lower=T, log=F
     ! Particularly *not* for log_p !

     ! Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
     !
     ! Note that we do want symmetry(0), lower/upper -> hence use y

  ELSE IF ((log_p .and. y .lt. 1.0d170) &   ! avoid underflow below
     !  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
     ! Then, make use of  Abramowitz & Stegun, 26.2.13, something like
     !
     !  xsq = x*x;
     !
     !  if(xsq * DBL_EPSILON < 1.)
     !      del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
     !  else
     !      del = 0.;
     !      *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
     !      *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
     !
     !  swap_tail;
     !
     !  [Yes, but xsq might be infinite.]
     !
          .or. (lower .and. -37.5193d0 .lt. x  .and.  x .lt. 8.2924d0) &
          .or. (upper .and. -8.2924d0  .lt. x  .and.  x .lt. 37.5193d0) &
  ) THEN

  ! Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5)
     xsq = 1.0d0 / (x * x)  ! (1./x)*(1./x) might be better
     xnum = p(6) * xsq
     xden = xsq
     DO i = 1, 4
        xnum = (xnum + p(i)) * xsq
        xden = (xden + q(i)) * xsq
     END DO
     temp = xsq * (xnum + p(5)) / (xden + q(5))
     temp = (M_1_SQRT_2PI - temp) / y

     xsq = AINT(x * SIXTEN) / SIXTEN
     del = (x - xsq) * (x + xsq)
     IF (log_p) THEN
        cum = (-xsq * xsq * 0.5d0) + (-del * 0.5d0) + log(temp)
        IF ((lower .and. x .gt. 0.0d0) .or. (upper .and. x .le. 0.0d0)) THEN
           ccum = log1p(-exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp)
        END IF
     ELSE
        cum = exp(-xsq * xsq * 0.5d0) * exp(-del * 0.5d0) * temp
        ccum = 1.0d0 - cum
     END IF

     IF (x .gt. 0.0d0) THEN ! swap  ccum <--> cum
        temp = cum
        IF (lower) cum = ccum
        ccum = temp
     END IF

  ELSE  ! large x such that probs are 0 or 1
     IF (x .gt. 0.0d0) THEN
        IF (log_p) THEN
           cum = 0.0d0;  ccum = NEG_INF
        ELSE
           cum = 1.0d0;  ccum = 0.0d0
        END IF
     ELSE
        IF (log_p) THEN
           cum = NEG_INF;  ccum = 0.0d0
        ELSE
           cum = 0.0d0;  ccum = 1.0d0
        END IF
     END IF
  END IF


  ! do not return "denormalized" -- we do in R
  IF (log_p) THEN
     IF(cum .gt. -DBL_MIN)   cum = -0.0d0
     IF(ccum  .gt. -DBL_MIN) ccum = -0.0d0
  ELSE
     IF (cum .lt. DBL_MIN)   cum = 0.0d0
     IF (ccum .lt. DBL_MIN)  ccum = 0.0d0
  END IF

  RETURN

END SUBROUTINE pnorm_both






!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000, 2003, 2011 The R Core Team
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
!	double log1p(double x);
!
!  DESCRIPTION
!
!	Compute the relative error logarithm.
!
!			log(1 + x)
!
!  NOTES
!
!	This code is a translation of the Fortran subroutine `dlnrel'
!	written by W. Fullerton of Los Alamos Scientific Laboratory.
!
!
!  Every currently known platform has log1p (which is C99), 
!    but NetBSD/OpenBSD were at least at one time inaccurate */


FUNCTION log1p (x) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  ! series for log1p on the interval -.375 to .375
  !				     with weighted error   6.35e-32
  !				      log weighted error  31.20
  !			    significant figures required  30.93
  !				 decimal places required  32.01

  DOUBLE PRECISION, PARAMETER :: alnrcs(43) = (/ &
	+.10378693562743769800686267719098d+1, &
	-.13364301504908918098766041553133d+0, &
	+.19408249135520563357926199374750d-1, &
	-.30107551127535777690376537776592d-2, &
	+.48694614797154850090456366509137d-3, &
	-.81054881893175356066809943008622d-4, &
	+.13778847799559524782938251496059d-4, &
	-.23802210894358970251369992914935d-5, &
	+.41640416213865183476391859901989d-6, &
	-.73595828378075994984266837031998d-7, &
	+.13117611876241674949152294345011d-7, &
	-.23546709317742425136696092330175d-8, &
	+.42522773276034997775638052962567d-9, &
	-.77190894134840796826108107493300d-10, &
	+.14075746481359069909215356472191d-10, &
	-.25769072058024680627537078627584d-11, &
	+.47342406666294421849154395005938d-12, &
	-.87249012674742641745301263292675d-13, &
	+.16124614902740551465739833119115d-13, &
	-.29875652015665773006710792416815d-14, &
	+.55480701209082887983041321697279d-15, &
	-.10324619158271569595141333961932d-15, &
	+.19250239203049851177878503244868d-16, &
	-.35955073465265150011189707844266d-17, &
	+.67264542537876857892194574226773d-18, &
	-.12602624168735219252082425637546d-18, &
	+.23644884408606210044916158955519d-19, &
	-.44419377050807936898878389179733d-20, &
	+.83546594464034259016241293994666d-21, &
	-.15731559416479562574899253521066d-21, &
	+.29653128740247422686154369706666d-22, &
	-.55949583481815947292156013226666d-23, &
	+.10566354268835681048187284138666d-23, &
	-.19972483680670204548314999466666d-24, &
	+.37782977818839361421049855999999d-25, &
	-.71531586889081740345038165333333d-26, &
	+.13552488463674213646502024533333d-26, &
	-.25694673048487567430079829333333d-27, &
	+.48747756066216949076459519999999d-28, &
	-.92542112530849715321132373333333d-29, &
	+.17578597841760239233269760000000d-29, &
	-.33410026677731010351377066666666d-30, &
	+.63533936180236187354180266666666d-31 &
  /)

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -1.0d99
  DOUBLE PRECISION, PARAMETER :: NaN = -99.99d0

  DOUBLE PRECISION, PARAMETER :: xmin = -0.999999985d0

  DOUBLE PRECISION :: chebyshev_eval, fn_val

  IF (x .eq. 0.0d0) THEN
     fn_val = 0.0d0        ! speed
     RETURN
  END IF

  IF (x .eq. -1.0d0) THEN
     fn_val = NEG_INF
     RETURN
  END IF

  IF (x .lt. -1.0d0) THEN
     fn_val = NAN
     RETURN
  END IF

  IF (abs(x) .le. 0.375d0) THEN
     ! Improve on speed (only); again give result accurate to IEEE double precision:
     IF (abs(x) .lt. 0.5d0 * DBL_EPSILON) THEN
        fn_val = x
        RETURN
     END IF

     IF ( (0.0d0 .lt. x .and. x .lt. 1.0d-8) .or. (-1.0d-9 .lt. x .and. x .lt. 0.0d0)) THEN
        fn_val = x * (1.0d0 - 0.5d0 * x)
        RETURN
     END IF

     fn_val = x * (1.0d0 - x * chebyshev_eval (x/0.375d0, alnrcs, 22))
     RETURN

  END IF
                           ! else
  IF (x .lt. xmin) THEN
     PRINT *, ' Warning! Argument of log1p is too near -1, '
     PRINT *, '  precision is less than half of usual precision! '
  END IF

  fn_val = log(1.0d0 + x)
  RETURN

END FUNCTION log1p




!  Mathlib : A C Library of Special Functions
!  Copyright (C) 2002 The R Core Team
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
!	double expm1(double x);
!
!  DESCRIPTION
!
!	Compute the Exponential minus 1
!
!			exp(x) - 1
!
!      accurately also when x is close to zero, i.e. |x| << 1
!
!  NOTES
!
!	As log1p(), this is a standard function in some C libraries,
!	particularly GNU and BSD (but is neither ISO/ANSI C nor POSIX).
!
!  We supply a substitute for the case when there is no system one.


FUNCTION expm1 (x) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: x

  DOUBLE PRECISION, PARAMETER :: DBL_EPSILON = 2.22044604925031d-16

  DOUBLE PRECISION :: y, a, log1p, fn_val

  a = abs(x)

  IF (a .lt. DBL_EPSILON) THEN
     fn_val = x
     RETURN
  END IF

  IF (a .gt. 0.697d0) THEN
     fn_val = exp(x) - 1.0d0  ! negligible cancellation
     RETURN
  END IF

  IF (a .gt. 1.0d-8) THEN
     y = exp(x) - 1.0d0
  ELSE                    ! Taylor expansion, more accurate in this range
     y = (x / 2.0d0 + 1.0d0) * x
  END IF

  ! Newton step for solving   log(1 + y) = x   for y :
  ! WARNING: does not work for y ~ -1: bug in 1.5.0

  y = y - (1.0d0 + y) * (log1p(y) - x)

  fn_val = y
  RETURN

END FUNCTION expm1


