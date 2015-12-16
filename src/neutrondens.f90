module neudens
use iso_fortran_env
use inputinterp
implicit none
real(real64) :: beta(7)      ! beta values for each decay group
real(real64) :: lambda(6)    ! half life constants
real(real64) :: pt           ! reactivity value
real(real64) :: ngen         ! neutron generation time

contains

subroutine neuden()
integer                 :: i                    ! counting variables
integer                 :: counter              ! iteration counter
integer                 :: info                 ! llapack error variable
real(real64)            :: h                    ! step size
real(real64)            :: y(7)                 ! matrix containing n(t) and c(t)
real(real64)            :: y0(7)                ! matrix containing initial values for n(t) and c(t)
real(real64)            :: yscale(7)            ! truncation error scaling value
real(real64)            :: g1(7)                ! variable of first equation
real(real64)            :: g2(7)                ! variable of second equation
real(real64)            :: g3(7)                ! variable of fourth equation
real(real64)            :: g4(7)                ! variable of fifth equation
real(real64)            :: dfdt(7)         
real(real64)            :: fyt(7)
real(real64)            :: RHS1(7)              ! right-hand-side of equation 1
real(real64)            :: RHS2(7)              ! right-hand-side of equation 2
real(real64)            :: RHS3(7)              ! right-hand-side of equation 3
real(real64)            :: RHS4(7)              ! right-hand-side of equation 4
real(real64)            :: nt                   ! neutron density
real(real64)            :: ptprev               ! reactivity from the previous iteration
real(real64)            :: t                    ! time
real(real64)            :: Ct(6)                ! delayed neutron precursers
real(real64)            :: identity(7,7)        ! identity matrix 
real(real64)            :: LHS(7,7)             ! left-hand-side of all linear equations
real(real64)            :: dfdy(7,7)
real(real64)            :: ipiv(7)              ! pivot vector used in llapack subroutines
real(real64), parameter :: gma = 0.5_real64        
real(real64), parameter :: a21 = 2.0_real64
real(real64), parameter :: a31 = 1.92_real64
real(real64), parameter :: a32 = 0.24_real64
real(real64), parameter :: c21 = -8.0_real64
real(real64), parameter :: c31 = 14.88_real64
real(real64), parameter :: c32 = 2.4_real64
real(real64), parameter :: c41 = -0.869_real64
real(real64), parameter :: c42 = -0.432_real64
real(real64), parameter :: c43 = -0.4_real64
real(real64), parameter :: b1 = 19.0_real64/9.0_real64
real(real64), parameter :: b2 = 0.5_real64
real(real64), parameter :: b3 = 25.0_real64/108.0_real64
real(real64), parameter :: b4 = 125.0_real64/108.0_real64
real(real64), parameter :: e1 = 17.0_real64/54.0_real64
real(real64), parameter :: e2 = 7.0_real64/36.0_real64
real(real64), parameter :: e3 = 0.0_real64
real(real64), parameter :: e4 = 125.0_real64/108.0_real64
real(real64), parameter :: c1 = 0.5_real64
real(real64), parameter :: c2 = -1.5_real64
real(real64), parameter :: c3 = 2.42_real64
real(real64), parameter :: c4 = 0.116_real64
real(real64), parameter :: a2 = 1.0_real64
real(real64), parameter :: a3 = 0.6_real64
real(real64)            :: err(7)
real(real64), parameter :: eps = 1E-2_real64       ! accepted error value
real(real64)            :: hretry
real(real64)            :: hnext
real(real64)            :: havg
real(real64)            :: errmax

external dgetrf, dgetrs
 

! for thermal neutrons
if (isThermal) then
  beta(1) = 0.000285_real64   ! beta of group 1
  beta(2) = 0.0015975_real64  ! beta of group 2
  beta(3) = 0.00141_real64    ! beta of group 3
  beta(4) = 0.0030525_real64  ! beta of group 4
  beta(5) = 0.00096_real64    ! beta of group 5
  beta(6) = 0.000195_real64   ! beta of group 6
  beta(7) = 0.0075_real64     ! Total Beta

  lambda(1) = 0.0127_real64   ! decay constant of group 1
  lambda(2) = 0.0317_real64   ! decay constant of group 2
  lambda(3) = 0.115_real64    ! decay constant of group 3
  lambda(4) = 0.311_real64    ! decay constant of group 4
  lambda(5) = 1.4_real64      ! decay constant of group 5
  lambda(6) = 3.87_real64     ! decay constant of group 6
  ngen      = 0.0005_real64   ! average neutron generation time
! for fast neutrons
else
  beta(1) = 0.0001672_real64   ! beta of group 1     
  beta(2) = 0.001232_real64    ! beta of group 2
  beta(3) = 0.0009504_real64   ! beta of group 3
  beta(4) = 0.001443_real64    ! beta of group 4
  beta(5) = 0.0004534_real64   ! beta of group 5
  beta(6) = 0.000154_real64    ! beta of group 6
  beta(7) = 0.0044_real64      ! Total Beta

  lambda(1) = 0.0129_real64    ! decay constant of group 1
  lambda(2) = 0.0311_real64    ! decay constant of gourp 2
  lambda(3) = 0.134_real64     ! decay constant of group 3
  lambda(4) = 0.331_real64     ! decay constant of group 4
  lambda(5) = 1.26_real64      ! decay constant of group 5
  lambda(6) = 3.21_real64      ! decay constant of group 6
  ngen      = 1E-7_real64      ! average neutron generation time  
end if

! initial values for nt and ct
nt = 1.0_real64
Ct(1) = (beta(1)/(ngen*lambda(1)))*nt 
Ct(2) = (beta(2)/(ngen*lambda(2)))*nt 
Ct(3) = (beta(3)/(ngen*lambda(3)))*nt 
Ct(4) = (beta(4)/(ngen*lambda(4)))*nt 
Ct(5) = (beta(5)/(ngen*lambda(5)))*nt 
Ct(6) = (beta(6)/(ngen*lambda(6)))*nt 

! initial y values
y0(1) = nt
y0(2) = Ct(1)
y0(3) = Ct(2)
y0(4) = Ct(3)
y0(5) = Ct(4)
y0(6) = Ct(5)
y0(7) = Ct(6)

counter = 1
t = get_start_time()
h = inputdata(2,1) - inputdata(1,1)
open(unit=50, file="nt.out")
open(unit=60, file="ct.out")

y = y0
! write initial values to file
write(50,'(ES15.3, ES15.6)') t, y(1)
  
write(60, '(ES10.3)', advance='no') t
do i = 1,6
  write(60, '(ES15.6)', advance='no') y(i+1)
end do
write(60,*)
! calculate y

do  ! Main loop
 ! assign values to matrix dfdy
  pt = get_reactivity(t) 
  dfdy(1,1) = (pt - beta(7))/ngen
  do i = 2,7
    dfdy(i,1) = beta(i-1)/ngen
    dfdy(1,i) = lambda(i-1)
  end do

  dfdy(2:7,2:7)  = 0.0_real64
  do i = 2,7
    dfdy(i,i) = -lambda(i-1)
  end do
  
  ! Build vector fyt
!!!  fyt = matmul(dfdy, y) !!!
  fyt = get_fyt(y,t)
  ! Calculate yscale
  yscale = abs(y) + abs(h * fyt) + 1E-30_real64
 
  ! create matrix dfdt
  if (counter == 1) then
    dfdt(1) = (y(1)/ngen)*(pt)
  else
    dfdt(1) = (y(1)/ngen)*(pt - ptprev)
  end if  
  dfdt(2:7) = 0.0_real64
      
  ! start building left hand side matrix
  
  ! Build identity matrix
  identity = 0.0_real64
  do i = 1,7
    identity(i,i) = 1.0_real64
  end do
  
  LHS = ((1.0_real64/(gma*h))*identity - dfdy)
  ! build first right hand side matrix
  RHS1 = fyt + h*c1*dfdt
  ! LU decomposition using LAPACK
  ! reference material for dgetrf can be found at: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html 
  call dgetrf(7, 7, LHS, 7, ipiv, info)
  if (info > 0) stop "Matrix is singular"
  if (info < 0) stop "Illegal value"
  
  ! Solve for g1 using LAPACK
  ! reference material for dgetrs can be found at: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS1, 7, info)
  if (info /= 0) stop "Solution of linear system failed"
  g1 = RHS1 
  
  ! Solve for g2
  fyt = get_fyt((y + a21*g1), (t + a2*h))
  RHS2 = fyt + h*c2*dfdt + (c21*g1)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS2, 7, info)
  g2 = RHS2

  ! Solve for g3
  fyt = get_fyt((y + a31*g1 + a32*g2), (t + a3*h))
  RHS3 = fyt + h*c3*dfdt + (c31*g1 + c32*g2)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS3, 7, info)
  g3 = RHS3

  ! Solve for g4
  fyt = get_fyt((y + a31*g1 + a32*g2), (t + a3*h))
  RHS4 = fyt + h*c4*dfdt + (c41*g1 + c42* g2 + c43*g3)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS4, 7, info)
  g4 = RHS4

!--------------------------------------------------------------------------
! Automatic step size calculation
  err = e1*g1+e2*g2+e3*g3+e4*g4
  errmax = maxval(abs(err/yscale))
   if (errmax > eps) then
     hretry = max(0.9_real64*h/(errmax**(1.0/3.0)),0.5_real64*h) ! note that x^-a == 1/x^a
     h = hretry
     cycle
   else
     if (errmax>0.1296) then
       hnext = 0.9_real64 / errmax**0.25
     else 
       hnext = 1.5_real64*h
     endif
     h = hnext
   end if
!----------------------------------------------------------------------    

! Check if the input file specifies a shorter time step
  if (nearest_time_step(t) < h) h = nearest_time_step(t)

  write(50,'(ES15.3, ES15.6)') t, y(1)
  
  write(60, '(ES10.3)', advance='no') t
  do i = 1,6
    write(60, '(ES15.6)', advance='no') y(i+1)
  end do
  write(60,*)


  ! Calculate next y
  y = y + (b1*g1 + b2*g2 + b3*g3 + b4*g4)


  havg = (havg + h)/counter
  counter = counter + 1

  if (t >= get_end_time()) then
    write(6, *)"havg = ", havg 
    exit ! exit main do loop
  end if
  
  ! save previous reactivity
  ptprev = get_reactivity(t)

  ! find next time value  
  t = t + h
  
end do

close(50)
close(60)

end subroutine neuden

function get_fyt(y_in,t_in)
! function to calculate the value of fyt
  real(real64), intent(in) :: y_in(7)
  real(real64), intent(in) :: t_in 
  real(real64)             :: df(7,7)
  real(real64)             :: get_fyt(7)
  integer                  :: i  
  pt = get_reactivity(t_in) 
  df(1,1) = (pt - beta(7))/ngen
  do i = 2,7
    df(i,1) = beta(i-1)/ngen
    df(1,i) = lambda(i-1)
  end do

  df(2:7,2:7)  = 0.0_real64
  do i = 2,7
    df(i,i) = -lambda(i-1)
  end do

  get_fyt = matmul(df,y_in)
end function get_fyt



end module neudens
