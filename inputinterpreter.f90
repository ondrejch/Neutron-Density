module inputinterp
use iso_fortran_env
implicit none

real(real64), protected :: inputdata(10,2) ! input array
integer, protected      :: n              ! length of input array


contains
  subroutine init_input_data(filename, n)
    character(100), intent(in)      :: filename        ! name of input file
    real(real64), allocatable       :: inputdata(:,:)  ! input data from file
    integer                         :: n               ! data length/ number of rows
    integer                         :: i
    allocate(inputdata(n,2))
    open(unit=1, file=filename)
    do i = 1,n
      read(1, *) inputdata(i,1), inputdata(i,2)
    end do         
  end subroutine init_input_data

  function get_reactivity(t)
    real(real64), intent(in) :: t               ! desired time
    real(real64)             :: get_reactivity  ! reactivity at desired time
    integer                  :: i               ! counting variable

    do i = 1,n-1
      if (t >= inputdata(i,1) .and. t < inputdata(i+1,1)) then
        get_reactivity = inputdata(i,2)
        exit
      else if (t > inputdata(n,1)) then
        get_reactivity = inputdata(n,1)
      else
        cycle
      end if
    end do 
  end function get_reactivity
end module inputinterp

