!------------------- program main --------------------------------------
! ======================================================================
! * Point kinetics equation solver for education and research purposes * 
! ======================================================================
! This code implements PKE solver expanding on 
!   Yang Xue, Jevremović Tatjana: "Revisiting the Rosenbrock numerical 
!   solutions of the reactor point kinetics equation with numerous examples",
!   Nuclear Technology and Radiation Protection 24, p. 3-12, (2009)
!   doi:10.2298/NTRP0901003Y
! by adding neutron source and thermal-reactivity feedback. 
!
! Authors: 
!   Dallas Moser <dmoser4@vols.utk.edu> 
!   Ondrej Chvala <ochvala@utk.edu>
! 
! License: GNU/GPL
!-----------------------------------------------------------------------
program main
use iso_fortran_env
use inputinterp
use neudens
implicit none
!
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
