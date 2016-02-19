!----------------------- module neudens -----------------------
! Module responsible for calculating the neutron density and  | 
! delayed neutron precursor concentrations given reactivity   | 
! and time.                                                   |
!                                                             |
! Authors:                                                    |
!   Dallas Moser <dmoser4@vols.utk.edu>                       |
!   Ondrej Chvala <ochvala@utk.edu>                           |
!                                                             |
! License: GNU/GPL                                            |
!--------------------------------------------------------------
module neudens
use iso_fortran_env
use inputinterp
use feedback
implicit none

real(real64) :: beta(7)      ! beta values for each decay group
real(real64) :: lambda(6)    ! half life constants
real(real64) :: pt           ! reactivity value
real(real64) :: ngen         ! neutron generation time

contains

!-------------------- neuden --------------------
!  Subroutine that calculates neutron density   |
!------------------------------------------------
subroutine neuden()
integer                 :: i, j                           ! Counting variables
integer                 :: counter                        ! Iteration counter
integer                 :: info                           ! llapack error variable
real(real64), parameter :: gma = 0.5_real64               ! Constants for differentiation                         
real(real64), parameter :: a21 = 2.0_real64               !             |
real(real64), parameter :: a31 = 1.92_real64              !             |
real(real64), parameter :: a32 = 0.24_real64              !             |
real(real64), parameter :: c21 = -8.0_real64              !             | 
real(real64), parameter :: c31 = 14.88_real64             !             |
real(real64), parameter :: c32 = 2.4_real64               !             |
real(real64), parameter :: c41 = -0.869_real64            !             |
real(real64), parameter :: c42 = -0.432_real64            !             |
real(real64), parameter :: c43 = -0.4_real64              !             | 
real(real64), parameter :: b1 = 19.0_real64/9.0_real64    !             |
real(real64), parameter :: b2 = 0.5_real64                !             |
real(real64), parameter :: b3 = 25.0_real64/108.0_real64  !             |
real(real64), parameter :: b4 = 125.0_real64/108.0_real64 !             |
real(real64), parameter :: c1 = 0.5_real64                !             |
real(real64), parameter :: c2 = -1.5_real64               !             |
real(real64), parameter :: c3 = 2.42_real64               !             |
real(real64), parameter :: c4 = 0.116_real64              !             |
real(real64), parameter :: a2 = 1.0_real64                !             |
real(real64), parameter :: a3 = 0.6_real64                !    ________\_/_________
real(real64), parameter :: e1 = 17.0_real64/54.0_real64   ! Constants used in error calculations
real(real64), parameter :: e2 = 7.0_real64/36.0_real64    !             |
real(real64), parameter :: e3 = 0.0_real64                !             |
real(real64), parameter :: e4 = 125.0_real64/108.0_real64 !    ________\_/_________
real(real64), parameter :: identity(8,8) = RESHAPE([(1.0_real64,(0.0_real64,i=1,8),j=1,8),1.0_real64],[8,8]) ! identity matrix 
real(real64), parameter :: eps = 1E-5_real64    ! accepted error value
real(real64)            :: h                    ! Time step size
real(real64)            :: nearest_h            ! Time step size to next input data
real(real64)            :: y(8)                 ! Matrix containing n(t) and c(t)
real(real64)            :: y0(8)                ! Matrix containing initial values for n(t) and c(t)
real(real64)            :: yscale(8)            ! Truncation error scaling value
real(real64)            :: g1(8)                ! Variable of first equation
real(real64)            :: g2(8)                ! Variable of second equation
real(real64)            :: g3(8)                ! Variable of fourth equation
real(real64)            :: g4(8)                ! Variable of fifth equation
real(real64)            :: dfdt(8)              ! Derivative of f(y,t) with respect to t
real(real64)            :: fyt(8)               ! Function f(y,t) used to represent the PRK equations in matrix form
real(real64)            :: RHS1(8)              ! Right-hand-side of equation 1
real(real64)            :: RHS2(8)              ! Right-hand-side of equation 2
real(real64)            :: RHS3(8)              ! Right-hand-side of equation 3
real(real64)            :: RHS4(8)              ! Right-hand-side of equation 4
real(real64)            :: nt                   ! Neutron density as a function of time
real(real64)            :: t                    ! Time
real(real64)            :: Ct(6)                ! Delayed neutron precursor concentrations
real(real64)            :: LHS(8,8)             ! Left-hand-side of all linear equations
real(real64)            :: dfdy(8,8)            ! Derivative of f(y,t) with respect to y
real(real64)            :: ipiv(8)              ! Pivot vector used in llapack subroutines
real(real64)            :: err(8)               ! Array of error values
real(real64)            :: hnext                ! Next time step size if small error
real(real64)            :: havg                 ! Average time step size
real(real64)            :: errmax               ! Max error in y
external dgetrf, dgetrs

! initialize the delayed neutron constants
call init_delayed_consts()

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
y0(8) = 1.0_real64
counter = 0
t = get_start_time()                ! store the starting time
h = inputdata(2,1) - inputdata(1,1) ! find the first time step size

! open files for writing
open(unit=50, file="nt.out") 
open(unit=60, file="ct.out")
open(unit=70, file="T.out")
y = y0
do  ! Main loop
  pt = get_reactivity(t)                               ! reactivity at current time
  pt = pt + get_feedback(y(1)+get_source(t)*h, t, h)   ! get temperature feedback - off by default

  ! assign values to matrix dfdy
  dfdy = 0.0_real64
  dfdy(1,1) = (pt - beta(7))/ngen
  dfdy(1,8) = get_source(t)
  do i = 2,7
    dfdy(i,1) = beta(i-1)/ngen
    dfdy(1,i) = lambda(i-1)
    dfdy(i,i) = -lambda(i-1)
  end do
  
  ! build vector fyt
  fyt = get_fyt(y,t)

  ! calculate yscale
  yscale = abs(y) + abs(h * fyt) + 1E-30_real64
 
  ! build matrix dfdt
  dfdt(1) = (y(1)/ngen)*get_reactivity_slope(t) + get_source_slope(t)
  dfdt(2:8) = 0.0_real64
  
  ! build left hand side of the equation.. this only has to be calculated once
  LHS = ((1.0_real64/(gma*h))*identity - dfdy)

  ! build first right hand side matrix
  RHS1 = fyt + h*c1*dfdt

  ! LU decomposition using LAPACK
  ! reference material for dgetrf can be found at: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html 
  call dgetrf(8, 8, LHS, 8, ipiv, info)
  if (info > 0) stop "Matrix is singular"
  if (info < 0) stop "Illegal value"
  
  ! solve for g1 using LAPACK
  ! reference material for dgetrs can be found at: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
  call dgetrs('N', 8, 1, LHS, 8, ipiv, RHS1, 8, info)
  if (info /= 0) stop "Solution of linear system failed"
  g1 = RHS1 
  
  ! solve for g2
  fyt = get_fyt((y + a21*g1), (t + a2*h))
  RHS2 = fyt + h*c2*dfdt + (c21*g1)/h
  call dgetrs('N', 8, 1, LHS, 8, ipiv, RHS2, 8, info)
  g2 = RHS2

  ! solve for g3
  fyt = get_fyt((y + a31*g1 + a32*g2), (t + a3*h))
  RHS3 = fyt + h*c3*dfdt + (c31*g1 + c32*g2)/h
  call dgetrs('N', 8, 1, LHS, 8, ipiv, RHS3, 8, info)
  g3 = RHS3

  ! solve for g4
  fyt = get_fyt((y + a31*g1 + a32*g2), (t + a3*h))
  RHS4 = fyt + h*c4*dfdt + (c41*g1 + c42* g2 + c43*g3)/h
  call dgetrs('N', 8, 1, LHS, 8, ipiv, RHS4, 8, info)
  g4 = RHS4

  ! shortest time step we should take
  nearest_h = nearest_time_step(t)
  
  ! Automatic step size calculation
  err = e1*g1+e2*g2+e3*g3+e4*g4                 ! calculate error values
  errmax = 0.0
  do i = 1, 8
    errmax = max(errmax, abs(err(i)/yscale(i))) ! calculate max error including truncation
  end do
  errmax = errmax/eps
  if (errmax < 1.0) then
    if (errmax > 0.1296_real64) then
      hnext = 0.9_real64*h/(errmax**0.25)
      h = hnext
    else
      hnext = 1.5_real64*h
      h = hnext
    end if
  else
    hnext = 0.9_real64*h/(errmax**(1.0/3.0))
    h = max(abs(hnext), 0.5_real64*abs(h))
    cycle
  end if

  ! Check if the input file specifies a shorter time step
  if (nearest_h < h) h = nearest_h

  ! write current values to file
  write(50,51) t, y(1)
  write(60, 61, advance='no') t, (y(i), i=2,7)
  write(60,*)

  ! calculate next y
  y = y + (b1*g1 + b2*g2 + b3*g3 + b4*g4)

  ! calculate average time step size  
  counter = counter + 1

  if (t >= get_end_time()) then
    exit ! exit main do loop
  end if
  
  ! find next time value  
  t = t + h
end do
havg = (get_end_time()-get_start_time())/counter
if(fDebug>0) print *, " havg = ", havg 

! close files
close(50)
close(60)
close(70)

! formats
51 FORMAT (ES13.6, ES25.16)
61 FORMAT (ES13.6,6ES25.16)

end subroutine neuden


!---------- get_fyt(y_in, t_in) ----------
! Function to calculate the value of fyt |
!-----------------------------------------
function get_fyt(y_in,t_in)
  real(real64), intent(in) :: y_in(8)
  real(real64), intent(in) :: t_in 
  real(real64)             :: df(8,8)
  real(real64)             :: get_fyt(8)
  integer                  :: i  
  pt = get_reactivity(t_in) 
  df = 0.0_real64
  df(1,1) = (pt - beta(7))/ngen
  df(1,8) = get_source(t_in)
  do i = 2,7
    df(i,1) = beta(i-1)/ngen
    df(1,i) = lambda(i-1)
    df(i,i) = -lambda(i-1)
  end do
  get_fyt = matmul(df,y_in)
end function get_fyt

!----------- init_delayed_consts() ----------
! Initializes constants of delayed neutrons |
!--------------------------------------------
subroutine init_delayed_consts()
! TODO: generalize to support other inputs
if (isThermal) then ! for thermal neutrons
  beta(1) = 0.000285_real64    ! beta of group 1
  beta(2) = 0.0015975_real64   ! beta of group 2
  beta(3) = 0.00141_real64     ! beta of group 3
  beta(4) = 0.0030525_real64   ! beta of group 4
  beta(5) = 0.00096_real64     ! beta of group 5
  beta(6) = 0.000195_real64    ! beta of group 6
  beta(7) = 0.0075_real64      ! Total Beta

  lambda(1) = 0.0127_real64    ! decay constant of group 1
  lambda(2) = 0.0317_real64    ! decay constant of group 2
  lambda(3) = 0.115_real64     ! decay constant of group 3
  lambda(4) = 0.311_real64     ! decay constant of group 4
  lambda(5) = 1.4_real64       ! decay constant of group 5
  lambda(6) = 3.87_real64      ! decay constant of group 6
  ngen      = 0.0005_real64    ! average neutron generation time
else ! for fast neutrons
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
end subroutine init_delayed_consts

end module neudens
