module matrices
use iso_fortran_env
implicit none

contains

subroutine mat(rtype, h, tstart, tend, pt)
character(1), intent(in)         :: rtype           ! reactor type
real(real64), intent(in)         :: h               ! step size
real(real64), intent(in)         :: tstart, tend    ! start and end time
real(real64), intent(in)         :: pt              ! rho
integer                          :: i, j, counter, info   ! counting variables
real(real64), dimension(7)       :: y, y0, yscale, g1, g2, g3, g4
real(real64), dimension(7)       :: dfdt, fyt, RHS1, RHS2, RHS3, RHS4
real(real64)                     :: nt, t           ! neutron density
real(real64), dimension(6)       :: Ct              ! delayed neutron precursers
real(real64), dimension(7,7)     :: identity, dfdy, LHS
real(real64), dimension(7)       :: beta, ipiv      ! Beta values
real(real64), dimension(6)       :: lambda          ! half life constants
real(real64)                     :: ngen            ! neutron generation time
real(real64)                     :: gamma           ! constant
real(real64)                     :: a21, a31, a32
real(real64)                     :: c21, c31, c32, c41, c42, c43
real(real64)                     :: b1, b2, b3, b4
real(real64)                     :: e1, e2, e3, e4
real(real64)                     :: c1, c2, c3, c4
real(real64)                     :: a2, a3

external dgetrf, dgetrs
gamma = 0.5

a21 = 2.0
a31 = 1.92
a32 = 0.24

c21 = -8.0
c31 = 14.88
c32 = 2.4
c41 = -0.896
c42 = -0.432
c43 = -0.4

b1 = 19.0/9.0
b2 = 0.5
b3 = 25.0/108.0
b4 = 125.0/108.0

e1 = 17.0/54.0
e2 = 7.0/36.0
e3 = 0.0
e4 = 125.0/108.0

c1 = 0.5
c2 = -1.5
c3 = 2.42
c4 = 0.116

a2 = 1.0
a3 = 0.6 

! for thermal neutrons
if (rtype == 't') then
  beta(1) = 0.000285   ! beta of group 1
  beta(2) = 0.0015975  ! beta of group 2
  beta(3) = 0.00141    ! beta of group 3
  beta(4) = 0.0030525  ! beta of group 4
  beta(5) = 0.00096    ! beta of group 5
  beta(6) = 0.000195   ! beta of group 6
  beta(7) = 0.0075     ! Total Beta

  lambda(1) = 0.0127   ! decay constant of group 1
  lambda(2) = 0.0317   ! decay constant of group 2
  lambda(3) = 0.115    ! decay constant of group 3
  lambda(4) = 0.311    ! decay constant of group 4
  lambda(5) = 1.4      ! decay constant of group 5
  lambda(6) = 3.87     ! decay constant of group 6
  ngen      = 0.0005   ! average neutron generation time
! for fast neutrons
else if (rtype == 'f') then
  beta(1) = 0.0001672   ! beta of group 1     
  beta(2) = 0.001232    ! beta of group 2
  beta(3) = 0.0009504   ! beta of group 3
  beta(4) = 0.001443    ! beta of group 4
  beta(5) = 0.0004534   ! beta of group 5
  beta(6) = 0.000154    ! beta of group 6
  beta(7) = 0.0044      ! Total Beta

  lambda(1) = 0.0129    ! decay constant of group 1
  lambda(2) = 0.0311    ! decay constant of gourp 2
  lambda(3) = 0.134     ! decay constant of group 3
  lambda(4) = 0.331     ! decay constant of group 4
  lambda(5) = 1.26      ! decay constant of group 5
  lambda(6) = 3.21      ! decay constant of group 6
  ngen      = 10.0E-7   ! average neutron generation time  
end if

! initial values for nt and ct
nt = 1.0
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

! Assume initial g values to be 0
counter = 0
t = tstart
open(unit=5, file="nt.out")
open(unit=6, file="ct.out")

do  

  ! Calculate y   
  if (counter == 0) then
    y = y0
  end if

  write(5,'(ES15.3, ES15.6)') t, y(1)
  
  write(6, '(ES10.3)', advance='no') t
  do i = 1,6
    write(6, '(ES15.6)', advance='no') y(i+1)
  end do
  write(6,*)

  ! assign values to matrix dfdy
  dfdy(1,1) = (pt - beta(7))/ngen
  do i = 2,7
    dfdy(i,1) = beta(i-1)/ngen
    dfdy(1,i) = lambda(i-1)
  end do
  do i = 2,7
    do j = 2,7
      dfdy(i,j) = 0.0
    end do
    dfdy(i,i) = -lambda(i-1)
  end do
  
  ! Build vector fyt
  fyt = matmul(dfdy, y)

  ! Calculate yscale
  yscale = abs(y) + abs(h * fyt) + 10**(-30)
 
  ! create matrix dfdt
  dfdt(1) = (nt/ngen)*(0.0) ! for non-constant rho dfdt = (nt/ngen)(dpt/dt)
  do i = 2,7
    dfdt(i) = 0.0
  end do
  
  ! start building left hand side matrix
  
  ! Build identity matrix
  do i = 1,7
    do j = 1,7
      identity(i,j) = 0.0
    end do
    identity(i,i) = 1.0
  end do
  
  LHS = ((1.0/(gamma*h))*identity - dfdy)

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

  ! Calculate next y
  y = y + (b1*g1 + b2*g2 + b3*g3 + b4*g4)


  
  
  counter = counter + 1
  t = t + h
  
  if (t > tend) then
    exit
  else
    continue
  end if

end do

close(5)
close(6)

end subroutine mat

end module matrices
