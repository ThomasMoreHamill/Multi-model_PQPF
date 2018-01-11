program test
real a(3)
data a /1., 2., 3./
integer, dimension(1) :: i
i = minloc(a)
print *,i 
stop
end program test