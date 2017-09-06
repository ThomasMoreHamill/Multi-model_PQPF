subroutine check(status)
use netcdf
integer, intent ( in) :: status
if (status /= nf90_noerr) then
  print *,'netcdf error ', (nf90_strerror(status))
  stop
end if
end subroutine check 
