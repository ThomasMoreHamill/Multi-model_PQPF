!-----------------------------------------------------------------------------
!
! Evaluation of the PDF of a Gaussian distribution
!
! This Fortran code is translated from the C code in /src/nmath of that comes
! with version 3.1 of the statistical programming language 'R':
!
!-----------------------------------------------------------------------------
!
!  Mathlib : A C Library of Special Functions
!  Copyright (C) 1998 Ross Ihaka
!  Copyright (C) 2000-2014 The R Core Team
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
!	double dnorm4(double x, double mu, double sigma, int give_log)
!	      {dnorm (..) is synonymous and preferred inside R}
!
!  DESCRIPTION
!
!	Compute the density of the normal distribution.
!



FUNCTION dnorm (xx, mu, sigma, give_log) RESULT(fn_val)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: xx, mu, sigma

  LOGICAL, INTENT(IN) :: give_log

  DOUBLE PRECISION, PARAMETER :: DBL_MAX = 1.7976931348623157d+308
  DOUBLE PRECISION, PARAMETER :: M_LN_SQRT_2PI = 0.918938533204672741780329736406d0
  DOUBLE PRECISION, PARAMETER :: M_1_SQRT_2PI = 0.398942280401432677939946059934d0
 
  DOUBLE PRECISION, PARAMETER :: NEG_INF = -1.0d99

  DOUBLE PRECISION :: x, x1, x2, fn_val

  IF (sigma .le. 0.0d0) THEN
     PRINT *, 'Error in function dnorm! Program terminated.'
     CALL EXIT
  END IF

  x = (xx - mu) / sigma

  x = abs(x)
  IF (x .ge. 2.0d0 * sqrt(DBL_MAX)) THEN
     IF (give_log) THEN
        fn_val = NEG_INF
     ELSE
        fn_val = 0.0d0
     END IF
     RETURN
  END IF

  IF (give_log) THEN
     fn_val = -(M_LN_SQRT_2PI + 0.5d0 * x * x + log(sigma))
     RETURN
  END IF

  fn_val = M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma

  RETURN

END FUNCTION dnorm

