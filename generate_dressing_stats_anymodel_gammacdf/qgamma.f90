!-----------------------------------------------------------------------------
!
! Evaluation of the quantile function of a gamma distribution
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000--2011 The R Core Team
!  Copyright (C) 2004--2009 The R Foundation
!  based on AS 91 (C) 1979 Royal Statistical Society
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, a copy is available at
!  http://www.r-project.org/Licenses/
!
!  DESCRIPTION
!
!	Compute the quantile function of the gamma distribution.
!
!  NOTES
!
!	This function is based on the Applied Statistics
!	Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
!
!	R core improvements:
!	o  lower_tail, log_p
!       o  non-trivial result for p outside [0.000002, 0.999998]
! 	o  p ~ 1 no longer gives +Inf; final Newton step(s)
!
!  REFERENCES
!
!	Best, D. J. and D. E. Roberts (1975).
!	Percentage Points of the Chi-Squared Distribution.
!	Applied Statistics 24, page 385.


FUNCTION qchisq_appr (p, nu, g, lower_tail, log_p, tol) RESULT(fn_val)

     ! g = log Gamma(nu/2), tol = EPS1

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: p, nu, g, tol

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0

  DOUBLE PRECISION, PARAMETER :: C7  = 4.67d0
  DOUBLE PRECISION, PARAMETER :: C8  = 6.66d0
  DOUBLE PRECISION, PARAMETER :: C9  = 6.73d0
  DOUBLE PRECISION, PARAMETER :: C10 = 13.32d0

  DOUBLE PRECISION, PARAMETER :: M_LN2 = 0.693147180559945309417232121458d0

  DOUBLE PRECISION :: alpha, a, c, ch, p1, p2, q, t, x, Clog, lgam1pa, fn_val
  DOUBLE PRECISION :: log1p, expm1, lgamma1p, qnorm

  ! test arguments and initialise

  ! #ifdef IEEE_754
  !    if (ISNAN(p) || ISNAN(nu))
  !        return p + nu;
  ! #endif

  IF ((log_p .and. p.gt.0.0d0) .or. (.not.log_p .and. (p.lt.0.0d0 .or. p.gt.1.0d0))) THEN
     PRINT *, 'Error in "qchisq_appr", probability not inside [0,1] !'
     fn_val = NaN
     RETURN
  END IF

  IF (nu .le. 0.0d0) THEN
     PRINT *, 'Error in "qchisq_appr", nu <= 0 !'
     fn_val = NaN
     RETURN
  END IF

  alpha = 0.5d0 * nu    ! = [pq]gamma() shape
  c = alpha - 1.0d0

  IF (lower_tail) THEN
     IF (log_p) THEN
        p1 = p
     ELSE
        p1 = log(p)
     END IF
  ELSE
     IF (log_p) THEN
        IF (p .gt. -M_LN2) THEN
           p1 = log(-expm1(p))
        ELSE
           p1 = log1p(-exp(p))
        END IF
     ELSE
        p1 = log1p(-p)
     END IF
  END IF

  IF (nu .lt. -1.24d0*p1) THEN
	! for small chi-squared
	! log(alpha) + g = log(alpha) + log(gamma(alpha)) =
	!        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
	!  catastrophic cancellation when alpha << 1

     IF (alpha .lt. 0.5d0) THEN
	lgam1pa = lgamma1p(alpha)
     ELSE
        lgam1pa = log(alpha) + g
     END IF

     ch = exp((lgam1pa + p1)/alpha + M_LN2)

  ELSE IF (nu .gt. 0.32d0) THEN  !  using Wilson and Hilferty estimate
     x = qnorm(p, 0.0d0, 1.0d0, lower_tail, log_p)
     p1 = 2.0d0 / (9.0d0*nu)
     ch = nu * (x*sqrt(p1) + 1.0d0-p1)**3

     ! approximation for p tending to 1:
     IF (ch .gt. 2.2d0*nu + 6.0d0) THEN
        IF (lower_tail) THEN
           IF (log_p) THEN
              IF (p .gt. -M_LN2) THEN
                 Clog = log(-expm1(p))
              ELSE
                 Clog = log1p(-exp(p))
              END IF
           ELSE
              Clog = log1p(-p)
           END IF
        ELSE
           IF (log_p) THEN
              Clog = p
           ELSE
              Clog = log(p)
           END IF
        END IF
        ch = -2.0d0*(Clog - c*log(0.5d0*ch) + g)
     END IF

  ELSE   ! "small nu" : 1.24*(-log(p)) <= nu <= 0.32
     ch = 0.4d0
     IF (lower_tail) THEN
        IF (log_p) THEN
           IF (p .gt. -M_LN2) THEN
              Clog = log(-expm1(p))
           ELSE
              Clog = log1p(-exp(p))
           END IF
        ELSE
           Clog = log1p(-p)
        END IF
     ELSE
        IF (log_p) THEN
           Clog = p
        ELSE
           Clog = log(p)
        END IF
     END IF
     a = Clog + g + c*M_LN2

     DO
        q = ch
        p1 = 1.0d0 / (1.0d0+ch*(C7+ch))
        p2 = ch*(C9+ch*(C8+ch))
        t = -0.5d0 +(C7+2.0d0*ch)*p1 - (C9+ch*(C10+3.0d0*ch))/p2
        ch = ch - (1.0d0- exp(a+0.5d0*ch)*p2*p1)/t
	IF (abs(q-ch) .le. tol * abs(ch)) EXIT
     END DO

  END IF

  fn_val = ch
  RETURN

END FUNCTION qchisq_appr




FUNCTION qgamma (pp, alpha, scale, lower_tail, give_log) RESULT(fn_val)
!		   shape = alpha
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: pp, alpha, scale

  LOGICAL, INTENT(IN) :: lower_tail, give_log

  DOUBLE PRECISION, PARAMETER :: POS_INF = 9999.0d99
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -9999.0d99
  DOUBLE PRECISION, PARAMETER :: NaN = 9999.0d0
  DOUBLE PRECISION, PARAMETER :: DBL_MIN = 2.2250738585072014d-308

  DOUBLE PRECISION, PARAMETER ::  i420 = 1.0d0 / 420.0d0
  DOUBLE PRECISION, PARAMETER :: i2520 = 1.0d0 / 2520.0d0
  DOUBLE PRECISION, PARAMETER :: i5040 = 1.0d0 / 5040.0d0

  DOUBLE PRECISION, PARAMETER :: d1_p = 1.0d0 + 1.0d-7
  DOUBLE PRECISION, PARAMETER :: d1_m = 1.0d0 - 1.0d-7

  DOUBLE PRECISION, PARAMETER :: EPS1 = 1.0d-2
  DOUBLE PRECISION, PARAMETER :: EPS2 = 5.0d-7                               ! final precision of AS 91
  DOUBLE PRECISION, PARAMETER :: EPS_N = 1.0d-15                             ! precision of Newton step / iterations
  DOUBLE PRECISION, PARAMETER :: LN_EPS = -36.043653389117156d0              ! = log(.Machine$double.eps) iff IEEE_754
  DOUBLE PRECISION, PARAMETER :: M_LN2 = 0.693147180559945309417232121458d0

  DOUBLE PRECISION, PARAMETER :: pMIN = 1.0d-100         ! was 0.000002 = 2e-6
  DOUBLE PRECISION, PARAMETER :: pMAX = (1.0d0-1.0d-14)  ! was (1-1e-12) and 0.999998 = 1 - 2e-6

  INTEGER, PARAMETER :: MAXIT = 1000  ! was 20

  DOUBLE PRECISION :: p, p_, a, b, c, g, ch, ch0, p1, p2, q, t, x
  DOUBLE PRECISION :: s1, s2, s3, s4, s5, s6
  DOUBLE PRECISION :: lgammafn, qnorm, pgamma, pgamma_raw, dgam, qchisq_appr
  DOUBLE PRECISION :: expm1, fn_val

  INTEGER i, max_it_Newton

  LOGICAL :: log_p

  log_p = give_log
  p = pp

  max_it_Newton = 1

  ! test arguments and initialise

  ! #ifdef IEEE_754
  !    if (ISNAN(p) || ISNAN(alpha) || ISNAN(scale))
  !	return p + alpha + scale;
  ! #endif

    IF (log_p) THEN
       IF (p .gt. 0.0d0) THEN
          PRINT *, ' Error! "qgamma" is avaluated with probability outside [0,1] <log> !', p
          fn_val = NaN
          RETURN
       END IF
       IF (p .eq. 0.0d0) THEN  ! upper bound
          IF (lower_tail) THEN
             fn_val = POS_INF
          ELSE
             fn_val = 0.0d0
          END IF
          RETURN
       END IF
       IF (p .eq. NEG_INF) THEN
          IF (lower_tail) THEN
             fn_val = 0.0d0
          ELSE
             fn_val = POS_INF
          END IF
          RETURN
       END IF

    ELSE   ! !log_p
       IF (p .lt. 0.0d0 .or. p .gt. 1.0d0) THEN
          PRINT *, ' Error! "qgamma" is avaluated with probability outside [0,1] !', p
          fn_val = NaN
          RETURN
       END IF
       IF (p .eq. 0.0) THEN
          IF (lower_tail) THEN
             fn_val = 0.0d0
          ELSE
             fn_val = POS_INF
          END IF
          RETURN
       END IF
       IF (p .eq. 1.0d0) THEN
          IF (lower_tail) THEN
             fn_val = POS_INF
          ELSE
             fn_val = 0.0d0
          END IF
          RETURN
       END IF
    END IF

    IF ( alpha.lt.0.0d0 .or. scale.le.0.0 ) THEN
       PRINT *, ' Error! Invalid scale or shape parameter in "qgamma" !'
       fn_val = NaN
       RETURN
    END IF

    IF ( alpha.eq.0.0d0 ) THEN  ! all mass at 0
       fn_val = 0.0d0
       RETURN
    END IF

    IF ( alpha.lt.1.0d-10 ) THEN   ! Warning seems unnecessary now:
       ! #ifdef _DO_WARN_qgamma_
       PRINT *,  'Warning! Value of shape parameter in "qgamma" is extremely small !'
       ! #endif
       max_it_Newton = 7   ! may still be increased below
    END IF

    IF (log_p) THEN            ! lower_tail prob (in any case)
       IF (lower_tail) THEN
          p_ = exp(p)
       ELSE
          p_ = - expm1(p)
       END IF
    ELSE
       IF (lower_tail) THEN
          p_ = p
       ELSE
          p_ = 0.5d0 - p + 0.5d0
       END IF
    END IF

    g = lgammafn(alpha)    ! log Gamma(v/2)


    ! ----- Phase I : Starting Approximation
    ch = qchisq_appr(p, 2.0d0*alpha, g, lower_tail, log_p, EPS1)

    IF ( ch.eq.POS_INF ) THEN    !    if(!R_FINITE(ch))
       ! forget about all iterations!
       max_it_Newton = 0
       GO TO 100
    END IF
  
    IF (ch .lt. EPS2)  THEN   ! Corrected according to AS 91; MM, May 25, 1999
       max_it_Newton = 20
       GO TO 100              ! and do Newton steps
    END IF

    ! FIXME: This (cutoff to {0, +Inf}) is far from optimal
    ! -----  when log_p or !lower_tail, but NOT doing it can be even worse

    IF ( p_.gt.pMAX .or. p_.lt.pMIN ) THEN
       ! did return ML_POSINF or 0.;  much better:
       max_it_Newton = 20
       GO TO 100              ! and do Newton steps
    END IF


    ! ----- Phase II: Iteration
    ! Call pgamma() [AS 239] and calculate seven term taylor series

    c = alpha - 1.0d0
    s6 = (120.0d0 + c*(346.0d0+127.0d0*c)) * i5040  ! used below, is "const"

    ch0 = ch            ! save initial approx
    DO i = 1, MAXIT
       q = ch
       p1 = 0.5d0 * ch
       p2 = p_ - pgamma_raw (p1, alpha, .true., .false.)
       ! #ifdef IEEE_754
       !     if(!R_FINITE(p2) || ch <= 0)
       ! #else
       !     if(errno != 0 || ch <= 0)
       ! #endif
       IF (p2 .eq. POS_INF .or. ch .le. 0.0d0) THEN
          ch = ch0
          max_it_Newton = 27
          GO TO 100
       END IF         ! was return ML_NAN

       t = p2 * exp(alpha*M_LN2+g+p1-c*log(ch))
       b = t/ch
       a = 0.5d0*t - b*c
       s1 = (210.0d0+ a*(140.0d0+a*(105.0d0+a*(84.0d0+a*(70.0d0+60.0d0*a))))) * i420
       s2 = (420.0d0+ a*(735.0d0+a*(966.0d0+a*(1141.0d0+1278.0d0*a)))) * i2520
       s3 = (210.0d0+ a*(462.0d0+a*(707.0d0+932.0d0*a))) * i2520
       s4 = (252.0d0+ a*(672.0d0+1182.0d0*a) + c*(294.0d0+a*(889.0d0+1740.0d0*a))) * i5040
       s5 = (84.0d0+2264.0d0*a + c*(1175.0d0+606.0d0*a)) * i2520

       ch = ch + t*(1.0d0+0.5d0*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))))
       IF (abs(q - ch) .lt. EPS2*ch)  GO TO 100
       IF (abs(q - ch) .gt. 0.1d0*ch) THEN  ! diverging? -- also forces ch > 0
          IF (ch .lt. q) THEN
             ch = 0.9d0 * q
          ELSE
             ch = 1.1d0 * q
          END IF
       END IF
    END DO


    ! no convergence in MAXIT iterations -- but we add Newton now...

100 CONTINUE

!     PR# 2214 : From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
!   --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision

!    With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
!   
!    Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
!   		    - also for lower_tail = FALSE	 or log_p = TRUE
!    		    - optionally *iterate* Newton

    x = 0.5d0 * scale * ch

    IF ( max_it_Newton.gt.0 ) THEN
       IF (.not.log_p) THEN   ! always use log scale
          p = log(p)
          log_p = .true.
       END IF
       
       IF ( x.eq.0.0d0 ) THEN
          x = DBL_MIN
          p_ = pgamma(x, alpha, scale, lower_tail, log_p)

          IF ( (lower_tail .and. p_.gt.p*d1_p) .or. ((.not.lower_tail) .and. p_.lt.p*d1_m) ) THEN
             fn_val = 0.0d0
             RETURN
          END IF      ! else:  continue, using x = DBL_MIN instead of 0

       ELSE
          p_ = pgamma(x, alpha, scale, lower_tail, log_p)

          IF ( p_.eq.NEG_INF ) THEN
             fn_val = 0.0d0   ! PR#14710
             RETURN
          END IF

          DO i = 1, max_it_Newton
             p1 = p_ - p;
             IF ( abs(p1).lt.abs(EPS_N*p) ) EXIT

             g = dgam(x, alpha, scale, log_p)
             IF ( (log_p .and. g.eq.NEG_INF) .or. ((.not.log_p) .and. g.eq.0.0d0) ) EXIT

	    ! else :
	    ! delta x = f(x)/f'(x);
	    ! if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
	    ! ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))

             IF (log_p) THEN
                t = p1*exp(p_-g)
             ELSE
                t = p1/g     ! = "delta x"
             END IF

             IF (lower_tail) THEN
                t = x - t
             ELSE
                t = x + t
             END IF

             p_ = pgamma (t, alpha, scale, lower_tail, log_p)
             IF (abs(p_-p) .gt. abs(p1) .or. (i .gt. 1 .and. abs(p_-p) .eq. abs(p1))) EXIT   ! <- against flip-flop

             x = t
          END DO
       END IF
    
    END IF

    fn_val = x
    RETURN

END FUNCTION qgamma
