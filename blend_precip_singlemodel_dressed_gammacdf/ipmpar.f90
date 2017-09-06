function ipmpar ( i )

!*****************************************************************************80
!
!! IPMPAR returns integer machine constants.
!
!  Discussion:
!
!    Input arguments 1 through 3 are queries about integer arithmetic.
!    We assume integers are represented in the N-digit, base A form
!
!      sign * ( X(N-1)*A^(N-1) + ... + X(1)*A + X(0) )
!
!    where 0 <= X(0:N-1) < A.
!
!    Then:
!
!      IPMPAR(1) = A, the base of integer arithmetic;
!      IPMPAR(2) = N, the number of base A digits;
!      IPMPAR(3) = A^N - 1, the largest magnitude.
!
!    It is assumed that the single and real ( kind = 8 ) floating
!    point arithmetics have the same base, say B, and that the
!    nonzero numbers are represented in the form
!
!      sign * (B^E) * (X(1)/B + ... + X(M)/B^M)
!
!    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
!    EMIN <= E <= EMAX.
!
!    Input argument 4 is a query about the base of real arithmetic:
!
!      IPMPAR(4) = B, the base of single and real ( kind = 8 ) arithmetic.
!
!    Input arguments 5 through 7 are queries about single precision
!    floating point arithmetic:
!
!     IPMPAR(5) = M, the number of base B digits for single precision.
!     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
!     IPMPAR(7) = EMAX, the largest exponent E for single precision.
!
!    Input arguments 8 through 10 are queries about real ( kind = 8 )
!    floating point arithmetic:
!
!     IPMPAR(8) = M, the number of base B digits for real ( kind = 8 ).
!     IPMPAR(9) = EMIN, the smallest exponent E for real ( kind = 8 ).
!     IPMPAR(10) = EMAX, the largest exponent E for real ( kind = 8 ).
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable FORTRAN Subroutine Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, 1978, pages 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired constant.
!
!    Output, integer ( kind = 4 ) IPMPAR, the value of the desired constant.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imach(10)
  integer ( kind = 4 ) ipmpar
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!     data imach( 1) /   2 /
!     data imach( 2) /  31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /  16 /
!     data imach( 5) /   6 /
!     data imach( 6) / -64 /
!     data imach( 7) /  63 /
!     data imach( 8) /  14 /
!     data imach( 9) / -64 /
!     data imach(10) /  63 /
!
!     Machine constants for the AT&T 3B SERIES, AT&T
!     PC 7300, AND AT&T 6300.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the BURROUGHS 1700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   33 /
!     data imach( 3) / 8589934591 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -256 /
!     data imach( 7) /  255 /
!     data imach( 8) /   60 /
!     data imach( 9) / -256 /
!     data imach(10) /  255 /
!
!     Machine constants for the BURROUGHS 5700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /    8 /
!     data imach( 5) /   13 /
!     data imach( 6) /  -50 /
!     data imach( 7) /   76 /
!     data imach( 8) /   26 /
!     data imach( 9) /  -50 /
!     data imach(10) /   76 /
!
!     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
!
!     data imach( 1) /      2 /
!     data imach( 2) /     39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /      8 /
!     data imach( 5) /     13 /
!     data imach( 6) /    -50 /
!     data imach( 7) /     76 /
!     data imach( 8) /     26 /
!     data imach( 9) / -32754 /
!     data imach(10) /  32780 /
!
!     Machine constants for the CDC 6000/7000 SERIES
!     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS OPERATING SYSTEM).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   48 /
!     data imach( 3) / 281474976710655 /
!     data imach( 4) /    2 /
!     data imach( 5) /   48 /
!     data imach( 6) / -974 /
!     data imach( 7) / 1070 /
!     data imach( 8) /   95 /
!     data imach( 9) / -926 /
!     data imach(10) / 1070 /
!
!     Machine constants for the CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS/VE OPERATING SYSTEM).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    48 /
!     data imach( 6) / -4096 /
!     data imach( 7) /  4095 /
!     data imach( 8) /    96 /
!     data imach( 9) / -4096 /
!     data imach(10) /  4095 /
!
!     Machine constants for the CRAY 1, XMP, 2, AND 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    47 /
!     data imach( 6) / -8189 /
!     data imach( 7) /  8190 /
!     data imach( 8) /    94 /
!     data imach( 9) / -8099 /
!     data imach(10) /  8190 /
!
!     Machine constants for the data GENERAL ECLIPSE S/200.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the HARRIS 220.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   23 /
!     data imach( 3) / 8388607 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   38 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HONEYWELL 600/6000
!     AND DPS 8/70 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   63 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     3 WORD real ( kind = 8 ) OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   39 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     4 WORD real ( kind = 8 ) OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   55 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 9000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -126 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the IBM 360/370 SERIES,
!     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
!     5/7/9 AND THE SEL SYSTEMS 85/86.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the IBM PC.
!
!      data imach(1)/2/
!      data imach(2)/31/
!      data imach(3)/2147483647/
!      data imach(4)/2/
!      data imach(5)/24/
!      data imach(6)/-125/
!      data imach(7)/128/
!      data imach(8)/53/
!      data imach(9)/-1021/
!      data imach(10)/1024/
!
!     Machine constants for the MACINTOSH II - ABSOFT
!     MACFORTRAN II.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the MICROVAX - VMS FORTRAN.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the PDP-11 FORTRAN SUPPORTING
!     32-BIT integer ARITHMETIC.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the SEQUENT BALANCE 8000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the SILICON GRAPHICS IRIS-4D
!     SERIES (MIPS R3000 PROCESSOR).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
  data imach( 1) /     2 /
  data imach( 2) /    31 /
  data imach( 3) / 2147483647 /
  data imach( 4) /     2 /
  data imach( 5) /    24 /
  data imach( 6) /  -125 /
  data imach( 7) /   128 /
  data imach( 8) /    53 /
  data imach( 9) / -1021 /
  data imach(10) /  1024 /
!
!     Machine constants for the UNIVAC 1100 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   60 /
!     data imach( 9) /-1024 /
!     data imach(10) / 1023 /
!
!     Machine constants for the VAX 11/780.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
  ipmpar = imach(i)

  return
end