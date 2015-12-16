module neudens
use iso_fortran_env
use inputinterp
implicit none

contains

subroutine neuden()
integer                          :: i, j            ! counting variables
integer                          :: counter         ! iteration counter
integer                          :: info            ! llapack error variable
real(real64)                     :: h               ! step size
real(real64), dimension(7)       :: y               ! matrix containing n(t) and c(t)
real(real64), dimension(7)       :: y0              ! matrix containing initial values for n(t) and c(t)
real(real64), dimension(7)       :: yscale          ! truncation error scaling value
real(real64), dimension(7)       :: g1              ! variable of first equation
real(real64), dimension(7)       :: g2              ! variable of second equation
real(real64), dimension(7)       :: g3              ! variable of fourth equation
real(real64), dimension(7)       :: g4              ! variable of fifth equation
real(real64), dimension(7)       :: dfdt         
real(real64), dimension(7)       :: fyt
real(real64), dimension(7)       :: RHS1            ! right-hand-side of equation 1
real(real64), dimension(7)       :: RHS2            ! right-hand-side of equation 2
real(real64), dimension(7)       :: RHS3            ! right-hand-side of equation 3
real(real64), dimension(7)       :: RHS4            ! right-hand-side of equation 4
real(real64)                     :: nt              ! neutron density
real(real64)                     :: pt              ! reactivity value
real(real64)                     :: t               ! time
real(real64), dimension(6)       :: Ct              ! delayed neutron precursers
real(real64), dimension(7,7)     :: identity        ! identity matrix 
real(real64), dimension(7,7)     :: LHS             ! left-hand-side of all linear equations
real(real64), dimension(7,7)     :: dfdy
real(real64), dimension(7)       :: beta            ! beta values for each decay group
real(real64), dimension(7)       :: ipiv            ! pivot vector used in llapack subroutines
real(real64), dimension(6)       :: lambda          ! half life constants
real(real64)                     :: ngen            ! neutron generation time
real(real64), parameter          :: gma = 0.5_real64        
real(real64), parameter          :: a21 = 2.0_real64
real(real64), parameter          :: a31 = 1.92_real64
real(real64), parameter          :: a32 = 0.24_real64
real(real64), parameter          :: c21 = -8.0_real64
real(real64), parameter          :: c31 = 14.88_real64
real(real64), parameter          :: c32 = 2.4_real64
real(real64), parameter          :: c41 = -0.869_real64
real(real64), parameter          :: c42 = -0.432_real64
real(real64), parameter          :: c43 = -0.4_real64
real(real64), parameter          :: b1 = 19.0_real64/9.0_real64
real(real64), parameter          :: b2 = 0.5_real64
real(real64), parameter          :: b3 = 25.0_real64/108.0_real64
real(real64), parameter          :: b4 = 125.0_real64/108.0_real64
real(real64), parameter          :: e1 = 17.0_real64/54.0_real64
real(real64), parameter          :: e2 = 7.0_real64/36.0_real64
real(real64), parameter          :: e3 = 0.0_real64
real(real64), parameter          :: e4 = 125.0_real64/108.0_real64
real(real64), parameter          :: c1 = 0.5_real64
real(real64), parameter          :: c2 = -1.5_real64
real(real64), parameter          :: c3 = 2.42_real64
real(real64), parameter          :: c4 = 0.116_real64
!real(real64), parameter          :: a2 = 1.0_real64
!real(real64), parameter          :: a3 = 0.6_real64
real(real64), dimension(7)       :: err
real(real64), parameter          :: eps = 10.0_real64**(-6)       ! accepted error value
real(real64)                     :: hretry, hnext, havg, errmax
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
  ngen      = 10.0E-7_real64   ! average neutron generation time  
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

!havg = h
counter = 1
t = get_start_time()
h = inputdata(3,1) - inputdata(2,1)
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
do  

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
  fyt = matmul(dfdy, y)

  ! Calculate yscale
  yscale = abs(y) + abs(h * fyt) + 10.0_real64**(-30)
 
  ! create matrix dfdt
  dfdt(1) = (y(1)/ngen)*(0.0_real64) ! for non-constant rho dfdt = (nt/ngen)(dpt/dt)
  dfdt(2:7) = 0.0_real64
      
  ! start building left hand side matrix
  
  ! Build identity matrix
  do i = 1,7
    do j = 1,7
      identity(i,j) = 0.0_real64
    end do
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
  fyt = matmul(dfdy,(y + a21*g1))
  RHS2 = fyt + h*c2*dfdt + (c21*g1)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS2, 7, info)
  g2 = RHS2

  ! Solve for g3
  fyt = matmul(dfdy,(y + a31*g1 + a32*g2))
  RHS3 = fyt + h*c3*dfdt + (c31*g1 + c32*g2)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS3, 7, info)
  g3 = RHS3

  ! Solve for g4
  fyt = matmul(dfdy,(y + a31*g1 + a32*g2))
  RHS4 = fyt + h*c4*dfdt + (c41*g1 + c42* g2 + c43*g3)/h
  call dgetrs('N', 7, 1, LHS, 7, ipiv, RHS4, 7, info)
  g4 = RHS4

!--------------------------------------------------------------------------
! automatic step size stuff
  err = e1*g1+e2*g2+e3*g3+e4*g4
  errmax = maxval(abs(err/yscale))
!  print *, "errmax =", err
!  print *, "yscale =", yscale
   if (errmax > eps) then
!     the paper uses this retry but it returns an error:
!     Floating-point exception - erroneous arithmetic operation

!     hretry = max(0.9_real64*h*(errmax)**(-1.0_real64/3.0_real64),5.0_real64*h)
     hretry = h*.1
     h = hretry
!     print *, h
     cycle
   else
     hnext = 1.5*h
     h = hnext
   end if
!----------------------------------------------------------------------    
  write(50,'(ES15.3, ES15.6)') t, y(1)
  
  write(60, '(ES10.3)', advance='no') t
  do i = 1,6
    write(60, '(ES15.6)', advance='no') y(i+1)
  end do
  write(60,*)


  ! Calculate next y
  y = y + (b1*g1 + b2*g2 + b3*g3 + b4*g4)


  havg = (havg + h)/counter
!  print *," counter =", counter, "reactivity =", get_reactivity(t), "step = ", h
  counter = counter + 1

  if (t >= get_end_time()) then
    write(6, *)"havg = ", havg 
    exit ! exit the do loop
  end if
  ! find next time value  
  t = t + h
!  h = inputdata(counter,1) - inputdata(counter-1,1) 
  
  
end do

close(50)
close(60)

end subroutine neuden

end module neudens
