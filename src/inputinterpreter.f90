module inputinterp
use iso_fortran_env
implicit none

real(real64), protected, allocatable :: inputdata(:,:) ! input array
integer, protected                   :: n              ! length of input array


contains
  subroutine init_input_data(filename, length)
    character(100), intent(in)      :: filename        ! name of input file
    integer                         :: length               ! data length/ number of rows
    integer                         :: i
    n = length
    allocate(inputdata(n,2))
    open(unit=10, file=filename)    
    do i = 1,n
      read(10, *) inputdata(i,1), inputdata(i,2)
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

  function get_start_time()
    real(real64) ::  get_start_time ! starting time 
    get_start_time = inputdata(1,1)
  end function get_start_time

  function get_end_time()
    real(real64) :: get_end_time  ! ending time
    get_end_time = inputdata(n,1)
  end function get_end_time

end module inputinterp

