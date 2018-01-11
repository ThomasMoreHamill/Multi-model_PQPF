
SUBROUTINE read_ppn_analog_locns_ccpa9(cmonth, nxa, nya, &
    nsupp, xlocations, ylocations, conusmask, &
    rlons, rlats, penalty, penalty_pa, penalty_ter, penalty_facet, &
    penalty_dist)

! compile with f2py -c -m read_ppn_analog_locns_ccpa9 read_ppn_analog_locns_ccpa9.f90

CHARACTER*3, INTENT(IN) :: cmonth

INTEGER, INTENT(IN) :: nxa, nya, nsupp

INTEGER, INTENT(OUT), DIMENSION(nxa,nya,nsupp) :: xlocations, ylocations
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlons, rlats
REAL, INTENT(OUT), DIMENSION(nxa,nya,nsupp) :: penalty, penalty_pa, &
     penalty_ter, penalty_facet, penalty_dist

! f2py intent(in) cmonth, nya, nxa, nsupp
! f2py intent(out) xlocations, ylocations, conusmask
! f2py intent(out) rlons, rlats, penalty, penalty_pa, penalty_ter
! f2py intent(out) penalty_facet, penalty_dist
! f2py depend(nxa,nya,nsupp) xlocations, ylocations, penalty, 
! f2py depend(nxa,nya,nsupp) penalty_pa, penalty_ter, penalty_facet, penalty_dist
! f2py depend(nxa,nya) conusmask, rlona, rlata

CHARACTER*120 infile

infile = '/Users/thamill/refcst2/test/ppn_analog_locns_ccpa9_'//cmonth//'.dat'
PRINT *, 'reading from ', TRIM(infile)
OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='unformatted')
READ (1) nxfin, nyfin, nxain, nyain, nsuppin
print *,'nxain, nyain = ',nxain, nyain
print *,'nxa, nya     = ',nxa, nya
print *,'nsuppin, nsupp = ',nsuppin, nsupp
READ (1) xlocations
READ (1) ylocations
READ (1) conusmask
READ (1) rlons
READ (1) rlats
READ (1) penalty
READ (1) penalty_pa
READ (1) penalty_ter
READ (1) penalty_facet
READ (1) penalty_dist
CLOSE (1)

RETURN
END SUBROUTINE read_ppn_analog_locns_ccpa9
