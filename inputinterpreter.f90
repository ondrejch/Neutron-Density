module inputinterp
use iso_fortran_env
implicit none
real(real64), common :: inputdata(:,:)
inputdata = init_input_data(filename, n)

contains
subroutine init_input_data(filename, n)
character, intent(in)      :: filename  ! name of input file
real(real64), allocatable  :: init_input_data ! input data from file
integer                    :: n         ! data length/ number of rows
integer                    :: i
allocate(inputdata(n,2))
open(unit=1, file=filename)
do i = 1,n
  read(1, *) init_input_data(i,1), init_input_data(i,2)
end do         
 

end module inputinterp

