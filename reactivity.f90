module reactivity
use iso_fortran_env
implicit none

contains 
subroutine react(in, tin, pout)
real(real64), intent(in)   :: in(:,:)   ! input matrix with columns of time and reactivity
real(real64), intent(in)   :: tin       ! input desired time
real(real64), intent(out)  :: pout      ! reactivity at input time
real(real64), dimension(:) :: t
real(real64), dimension(:) :: p
integer                    :: i


do i = 1, size(in,1)
  t(i) = in(i,1)
  p(i) = in(i,2)
end do



