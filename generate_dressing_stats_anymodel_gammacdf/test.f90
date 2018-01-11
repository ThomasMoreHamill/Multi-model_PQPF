program test
character*5 cthresh

cthresh = '1.0'
READ (TRIM(cthresh),*) rthresh
PRINT cthresh,' ', rthresh

cthresh = '10.0'
READ (TRIM(cthresh),*) rthresh
PRINT cthresh,' ', rthresh

cthresh = '100.0'
READ (TRIM(cthresh),*) rthresh
PRINT cthresh,' ', rthresh

stop
end