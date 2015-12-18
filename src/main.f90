program main

!---------- program main --------------
! Not sure what to put here:
! under construction
!                  :
!                \ _ /
!              -= (_) =-
!        ,________|_]
!        | // // // |
!        |//_//_//_/|
!        ||""""""""||
!        ||________||
!        |.========.|
!        ||        ||     
!--------------------------------------

use iso_fortran_env
use inputinterp
use neudens
implicit none

character(100)     :: filename     ! name of input file
integer            :: nCLP         ! number of command line parameter

! Get the number of command line arguments
nCLP = command_argument_count()
if (nCLP.NE.1) call code_usage()   ! Code needs exactly one argument    

! Get file name
call get_command_argument(1, filename) 

! Initialize input
call init_input_data(TRIM(filename))

! Run the calculation
call neuden()

contains
subroutine code_usage() 
! Prints how to run the code
print *, " This program does  !!TODO!! "
print *, " It expects a single argument: the input file name"
print *, " The format of the input file name is: !! TODO!!!"
end subroutine code_usage

end program main
