program main
! work-in-progress
! input rho
! input times
! input tend
! input initial step size
! determine if step size is automatic
! determine reactor type
use iso_fortran_env
use neudens

implicit none
character(1)         :: rtype, Y_N     ! Reactor type and answer variable
character(100)       :: filename       ! name of input file
real(real64)         :: stepsize       ! time step size (h)
integer              :: length         ! length of input file


write (*,10) "Please input reactor type ('f' for fast, 't' for thermal)"
read *, rtype

write (*, 10) "PLease input filename:"
read *, filename

write (*, 10) "Please enter the number of time steps used:"
read *, length

write (*,10) "Please input step size, input 0 for auto stepsize:"
read *, stepsize


! Format:
10 format (/, A)

call neuden(rtype, stepsize, filename, length)


end program main
