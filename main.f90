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
real(real64)         :: pt             ! Rho value
real(real64)         :: stepsize       ! time step size (h)
real(real64)         :: tstart, tend   ! start and end time

write (*,10) "Please input reactor type ('f' for fast, 't' for thermal)"
read *, rtype

write (*,10) "Would you like to load input for rho? (Y/N)"
read *, Y_N
! if yes, load file
if (Y_N == 'N') then
  write (*,10) "Would you like to input constant rho? (Y/N)"
  read *, Y_N
  if (Y_N == 'Y') then
    write(*,10) "Please input constant rho:"
    read *, pt
  else if (Y_N == 'N') then
    write (*,10) "Would you like to use default rho? (Y/N)"
    read *, Y_N
    if (Y_N == 'Y') then
      pt = 0.0022
    else
      stop "No rho value determined"
    end if
  end if
end if    

write (*,10) "Please input step size:"
read *, stepsize

write (*,10) "Please input start and end time"
read *, tstart, tend 


! Format:
10 format (/, A)

call neuden(rtype, stepsize, tstart, tend, pt)


end program main
