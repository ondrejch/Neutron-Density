module inputinterp
! Initializes ... 
! TODO comments
!
use iso_fortran_env
implicit none

integer, parameter                   :: fDebug = 3     ! debugging level
real(real64), protected, allocatable :: inputdata(:,:) ! input array
integer, protected                   :: nRecords = -1  ! length of input array
logical, protected                   :: isThermal      ! reactor type

contains
  subroutine init_input_data(filename)
  ! Reads in input data file
    character(100), intent(in)      :: filename        ! name of input file
    integer                         :: i, ioerr
    real(real64)                    :: tmp1, tmp2      ! temporary vars for input validation
    if (fDebug>1) print *, "[DEBUG] Reading input from file: ", filename
    open(unit=10, file=filename, status="old", action="read", iostat=ioerr)
    if (ioerr/=0) stop "Input data file does not exist, bailing out!"
    read(10, *,iostat=ioerr) isThermal                 ! first record: is the core thermal?
    if (ioerr/=0) stop "Input data file reading error, bailing out!"
    if (fDebug>1) print *, "[DEBUG] Thermal spectrum problem? ", isThermal
    ! Figure out the number of records
    nRecords = 0
    do
     read(10, *,iostat=ioerr) tmp1, tmp2
     if (ioerr>0) stop "Input data file reading error, bailing out!"
     if (ioerr<0) exit                                 ! reached EOF
     nRecords = nRecords + 1
    end do
    if (fDebug>1) print *, "[DEBUG] There appears to be ",nRecords," state points in the input file"
    if (nRecords<2) stop "Not enough input data, bailing out!"
    ! Now  read n records
    rewind(10)
    read(10,*) ! skip the first line
    allocate(inputdata(nRecords,2))
    do i = 1, nRecords
      read(10, *) inputdata(i,1), inputdata(i,2)
    end do
  end subroutine init_input_data

  function get_reactivity(t)
  ! Returns reactivity at time t
    real(real64), intent(in) :: t               ! desired time
    real(real64)             :: get_reactivity  ! reactivity at desired time
    integer                  :: i               ! counting variable
    do i = 1, nRecords-1
      if (t >= inputdata(i,1) .and. t < inputdata(i+1,1)) then
        get_reactivity = inputdata(i,2)
        exit
      else if (t > inputdata(nRecords,1)) then
        get_reactivity = inputdata(nRecords,1)
      else
        cycle
      end if
    end do 
  end function get_reactivity

  function nearest_time_step(t)
  ! Returns distance to the next time step specified in the input file
    real(real64), intent(in) :: t                 ! current time
    real(real64)             :: nearest_time_step ! distance to the next time step
    integer                  :: i
    do i = 1, nRecords-1
       if (t >= inputdata(i,1) .and. t < inputdata(i+1,1)) then
          nearest_time_step = inputdata(i+1,1) -t
          exit
       endif
    end do
  end function nearest_time_step
            
  function get_start_time()
    real(real64) ::  get_start_time ! starting time 
    get_start_time = inputdata(1,1)
  end function get_start_time

  function get_end_time()
    real(real64) :: get_end_time  ! ending time
    get_end_time = inputdata(nRecords,1)
  end function get_end_time

end module inputinterp

