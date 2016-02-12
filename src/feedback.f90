!------------------------ module feedback ------------------------
! Module containing functions responsible for simulating thermal 
! reactivity feedback                                            
!
! Authors: 
!   Dallas Moser <dmoser4@vols.utk.edu> 
!   Ondrej Chvala <ochvala@utk.edu>
! 
! License: GNU/GPL
!-----------------------------------------------------------------
module feedback
use iso_fortran_env
use inputinterp
implicit none
!
real(real64), parameter :: temp_equil    = 300.0     ! 300K of equilibrium temperature
real(real64)            :: reactor_temp = temp_equil ! Reactor temperature, starts at equilibrium
!
contains

!-------------- get_feedback(flux,dt) ---------------
! Returns reactivity feedback based on reactor flux |
! and time step size                                |
!----------------------------------------------------
function get_feedback(flux, t, dt)
! This really should use its own time step...
!
! Constants are made up ... fixme units
real(real64)            :: get_feedback                ! reactivity feedback
real(real64)            :: flux                        ! neutron flux [nt][cm-1 s-1]
real(real64)            :: t                           ! time [s]
real(real64)            :: dt                          ! time step [s]
real(real64), parameter :: macro_cross   = 585E-24     ! Macroscopic cross-section [cm-1]
! real(real64), parameter :: energy_abs = 200.0          ! Energy absorbed in the fuel [ MeV ]
real(real64), parameter :: heat_per_fission = .35E-10  ! Heat per fission [W s]
real(real64), parameter :: volume = 100.0**3           ! Volume of one cubic meter core [cm3]
real(real64), parameter :: fuel_specific_heat = 120    ! Specific heat of U-235 [W s kg-1 K-1]
real(real64), parameter :: fuel_density = 0.01905      ! Density of U-235 [kg cm-3]
real(real64), parameter :: alpha_temp    = -3.2969E-4  ! Temperature coefficient of reactivity [K-1], note pcm = 1E-5
real(real64), parameter :: heat_capacity = 4.18E3      ! Heat capacity for water @ 300K [J kg-1 K-1]
real(real64), parameter :: ht_coeff      = .62         ! Thermal conductivity of water @ 300K [W m-1 K-1]
real(real64), parameter :: sb            = 5.670367E-8 ! Stefan–Boltzmann constant [W m−2 K−4]
real(real64), parameter :: emisivity     = 0.94        ! Rough Concrete 
real(real64)            :: fuel_heat_capacity          ! Heat capacity of the fuel [W s K-1 Volume-1] 
real(real64)            :: power                       ! Power generated by fission [W]
real(real64)            :: energy_temp_conversion      ! Conversion factor for energy to temperature [K W-1 s-1]
real(real64)            :: fission_heating             ! Heating by fission
real(real64)            :: ht_conduction               ! Conduction heat loss
real(real64)            :: ht_radiation                ! Radiation heat loss
real(real64)            :: delta_reactor_temp          ! Reactor temperature change

! Heating up by neutrons 
power = macro_cross*flux*volume*heat_per_fission

fuel_heat_capacity = fuel_specific_heat*fuel_density

energy_temp_conversion = (1.0_real64/fuel_heat_capacity)*(1.0_real64/volume)

fission_heating = power*energy_temp_conversion


! Heat loss by conduction (natural convection)
ht_conduction = ht_coeff*(reactor_temp - temp_equil) 

! Heat loss by radiation 
ht_radiation = emisivity*sb*(reactor_temp**4 - temp_equil**4)

! Calculate new reactor temperature
delta_reactor_temp = dt*(fission_heating - ht_conduction - ht_radiation)/heat_capacity
reactor_temp = reactor_temp + delta_reactor_temp
write(70,*) t, reactor_temp

! Feedback
get_feedback = delta_reactor_temp * alpha_temp * get_reactivity(t)

end function get_feedback

end module feedback
