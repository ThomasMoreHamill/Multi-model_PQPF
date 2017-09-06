!-----------------------------------------------------------------------------
!
! Evaluation of the quantile function of a Gaussian distribution
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998       Ross Ihaka
!  Copyright (C) 2000--2005 The R Core Team
!  based on AS 111 (C) 1977 Royal Statistical Society
!  and   on AS 241 (C) 1988 Royal Statistical Society
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
!	double qnorm5(double p, double mu, double sigma,
!		      int lower_tail, int log_p)
!            {qnorm (..) is synonymous and preferred inside R}
!
!  DESCRIPTION
!
!	Compute the quantile function for the normal distribution.
!
!	For small to moderate probabilities, algorithm referenced
!	below is used to obtain an initial approximation which is
!	polished with a final Newton step.
!
!	For very large arguments, an algorithm of Wichura is used.
!
!  REFERENCE
!
!	Beasley, J. D. and S. G. Springer (1977).
!	Algorithm AS 111: The percentage points of the normal distribution,
!	Applied Statistics, 26, 118-121.
!
!      Wichura, M.J. (1988).
!      Algorithm AS 241: The Percentage Points of the Normal Distribution.
!      Applied Statistics, 37, 477-484.


FUNCTION qnorm (p, mu, sigma, lower_tail, log_p) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: p, mu, sigma

  LOGICAL, INTENT(IN) :: lower_tail, log_p

  DOUBLE PRECISION, PARAMETER :: POS_INF = 1.0d99
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -1.0d99
  DOUBLE PRECISION, PARAMETER :: NaN = -99.99d0

  DOUBLE PRECISION :: p_, q, r, val, expm1, fn_val

  ! #ifdef IEEE_754
  !    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
  !	return p + mu + sigma;
  ! #endif

  IF (log_p) THEN
     IF (p .gt. 0.0d0) THEN
        PRINT *, ' Error! "qnorm" is avaluated with probability outside [0,1] !'
        fn_val = NaN
        RETURN
     END IF
     IF (p .eq. 0.0d0) THEN  ! upper bound
        IF (lower_tail) THEN
           fn_val = POS_INF
        ELSE
           fn_val = NEG_INF
        END IF
        RETURN
     END IF
     IF (p .eq. NEG_INF) THEN
        IF (lower_tail) THEN
           fn_val = NEG_INF
        ELSE
           fn_val = POS_INF
        END IF
        RETURN
     END IF

  ELSE   ! !log_p
     IF (p .lt. 0.0d0 .or. p .gt. 1.0d0) THEN
        PRINT *, ' Error! "qnorm" is avaluated with probability outside [0,1] !'
        fn_val = NaN
        RETURN
     END IF
     IF (p .eq. 0.0) THEN
        IF (lower_tail) THEN
           fn_val = NEG_INF
        ELSE
           fn_val = POS_INF
        END IF
        RETURN
     END IF
     IF (p .eq. 1.0d0) THEN
        IF (lower_tail) THEN
           fn_val = POS_INF
        ELSE
           fn_val = NEG_INF
        END IF
        RETURN
     END IF
  END IF

  IF (sigma .lt. 0.0d0) THEN
     PRINT *, ' Error! "qnorm" is avaluated with sigma < 0 !'
     fn_val = NaN
     RETURN
  END IF

  IF (sigma .eq. 0.0d0) THEN
     fn_val = mu
     RETURN
  END IF

  IF (log_p) THEN            ! real lower_tail prob. p
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

  q = p_ - 0.5d0

! -- use AS 241 --- 
! double ppnd16_(double *p, long *ifault)
!      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
!
!        Produces the normal deviate Z corresponding to a given lower
!        tail area of P; Z is accurate to about 1 part in 10**16.
!
!        (original fortran code used PARAMETER(..) for the coefficients
!         and provided hash codes for checking them...)

  IF (abs(q) .le. 0.425d0) THEN  ! 0.075 <= p <= 0.925
     r = 0.180625d0 - q * q
     val = q * (((((((r * 2509.0809287301226727d0 + &
                       33430.575583588128105d0) * r + 67265.770927008700853d0) * r + &
                     45921.953931549871457d0) * r + 13731.693765509461125d0) * r + &
                   1971.5909503065514427d0) * r + 133.14166789178437745d0) * r + &
                 3.387132872796366608d0) &
            / (((((((r * 5226.495278852854561d0 + &
                     28729.085735721942674d0) * r + 39307.89580009271061d0) * r + &
                   21213.794301586595867d0) * r + 5394.1960214247511077d0) * r + &
                 687.1870074920579083d0) * r + 42.313330701600911252d0) * r + 1.0d0)

  ELSE  ! closer than 0.075 from {0,1} boundary

     ! r = min(p, 1-p) < 0.075
     IF (q .gt. 0.0d0) THEN
        IF (log_p) THEN
           IF (lower_tail) THEN
              r = -expm1(p)
           ELSE
              r = exp(p)
           END IF
        ELSE
           IF (lower_tail) THEN
              r = 0.5d0 - p + 0.5d0
           ELSE
              r = p
           END IF
        END IF
     ELSE
        r = p_   ! = R_DT_Iv(p) ^=  p
     END IF

     IF (log_p .and. ((lower_tail .and. q .le. 0.0d0) .or. (.not.lower_tail .and. q .gt. 0.0d0))) THEN
        r = sqrt(-p)
     ELSE
        r = sqrt(-log(r))
     END IF
     ! r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 )

     IF (r .le. 5.0d0) THEN  ! <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11
        r = r - 1.6d0
        val = (((((((r * 7.7454501427834140764d-4 + &
                       0.0227238449892691845833d0) * r + 0.24178072517745061177d0) * &
                     r + 1.27045825245236838258d0) * r + &
                    3.64784832476320460504d0) * r + 5.7694972214606914055d0) * &
                  r + 4.6303378461565452959d0) * r + &
                 1.42343711074968357734d0) &
                / (((((((r * &
                         1.05075007164441684324d-9 + 5.475938084995344946d-4) * &
                        r + 0.0151986665636164571966d0) * r + &
                       0.14810397642748007459d0) * r + 0.68976733498510000455d0) * &
                     r + 1.6763848301838038494d0) * r + &
                    2.05319162663775882187d0) * r + 1.0d0)

     ELSE  ! very close to  0 or 1
        r = r -5.0d0
        val = (((((((r * 2.01033439929228813265d-7 + &
                       2.71155556874348757815d-5) * r + &
                      0.0012426609473880784386d0) * r + 0.026532189526576123093d0) * &
                    r + 0.29656057182850489123d0) * r + &
                   1.7848265399172913358d0) * r + 5.4637849111641143699d0) * &
                 r + 6.6579046435011037772d0) &
                / (((((((r * &
                         2.04426310338993978564d-15 + 1.4215117583164458887d-7)* &
                        r + 1.8463183175100546818d-5) * r + &
                       7.868691311456132591d-4) * r + 0.0148753612908506148525d0) &
                     * r + 0.13692988092273580531d0) * r + &
                    0.59983220655588793769d0) * r + 1.0d0)
     END IF

     IF (q .lt. 0.0d0) val = -val    ! return (q >= 0.)? r : -r

  END IF

  fn_val = mu + sigma * val
  RETURN

END FUNCTION qnorm


