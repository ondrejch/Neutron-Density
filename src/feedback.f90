!------------------------- module feedback ---------------------------
! Module containing functions responsible for simulating thermal     |
! reactivity feedback. The normalization is done by a 1m x 1m block  |
! of U-235 submerged in an infinite pool of water. The block has an  |
! array of 50 x 50 holes with a diameter of 1cm and with 2cm spacing.|                                             
!                                                                    |
! Cross section data can be found here: http://atom.kaeri.re.kr/     |
!                                                                    |
! Authors:                                                           |
!   Dallas Moser <dmoser4@vols.utk.edu>                              |
!   Ondrej Chvala <ochvala@utk.edu>                                  |
!                                                                    |
! License: GNU/GPL                                                   |
!---------------------------------------------------------------------
module feedback
use iso_fortran_env
use inputinterp
implicit none

real(real64), parameter :: temp_equil   = 300.0      ! Temperature of water coolant / equilibrium (300K)
real(real64)            :: reactor_temp = temp_equil ! Reactor temperature, starts at equilibrium

contains

!------------------ get_feedback(flux,dt) -------------------
! Returns thermal reactivity feedback based on reactor flux |
! and time step size. The equations for core heat           |
! generation, energy to temperature conversion and related  |
! properties or parameters used were found here:            |               
!                                                           |
!   Karl O. Ott and Robert J. Neuhold, "Indrocutory Nuclear |
!   Reactor Dynamics", p. 228-234 (1985)                    |
!   ISBN: 0-89448-029-4                                     |
!------------------------------------------------------------
function get_feedback(n_density, t, dt)
real(real64), intent(in) :: n_density                                   ! Neutron density n(t) [neutron cm-3]
real(real64), intent(in) :: t                                           ! Time [s]
real(real64), intent(in) :: dt                                          ! Time step [s]
real(real64), parameter  :: fast_v             = 2.0E9                  ! Velocity of fast neutrons [cm s-1]
real(real64), parameter  :: thermal_v          = 2.2E5                  ! Velocity of thermal neutrons [cm s-1]            
real(real64), parameter  :: micro_cross_f      = 1.28672E-24            ! Fast microscopic cross-section @ 2 MeV [cm2]
real(real64), parameter  :: micro_cross_t      = 585E-24                ! Thermal microscopic cross-section @ 0.0253 eV [cm2]
real(real64), parameter  :: a_density          = (6.022E23/235.0)*19.05 ! Atomic density [atoms cm-3]
real(real64), parameter  :: heat_per_fission   = .35E-10                ! Heat per fission [J per fission]
real(real64), parameter  :: volume             = 100.0**3               ! Volume of one cubic meter core [cm3]
real(real64), parameter  :: fuel_specific_heat = 120                    ! Specific heat of U-235 [W s kg-1 K-1]
real(real64), parameter  :: fuel_density       = 0.01905                ! Density of U-235 [kg cm-3]
!real(real64), parameter  :: alpha_temp         = -3.2969E-4             ! Temperature coefficient of reactivity [K-1], note pcm = 1E-5
real(real64), parameter  :: alpha_temp         = -3.0            ! Temperature coefficient of reactivity [K-1], note pcm = 1E-5
real(real64), parameter  :: sb                 = 5.670367E-8            ! Stefan–Boltzmann constant [W m−2 K−4]
real(real64), parameter  :: emisivity          = 0.94                   ! Rough Concrete 

real(real64)             :: get_feedback                  ! reactivity feedback
real(real64)             :: micro_cross                   ! The determined microscopic cross section [cm]
real(real64)             :: n_velocity                    ! The determined neutron velocity [cm s-1]
real(real64)             :: fuel_heat_capacity            ! Heat capacity of the fuel [J K-1 cm-3] 
real(real64)             :: power                         ! Power generated by fission [W]
real(real64)             :: energy_temp_conversion        ! Conversion factor for energy to temperature [K W-1 s-1]
real(real64)             :: ht_conduction                 ! Conduction heat loss
real(real64)             :: ht_radiation                  ! Radiation heat loss
real(real64)             :: delta_reactor_temp            ! Reactor temperature change

! Determine neutron velocity through reactor type
if (isThermal) then
  n_velocity = thermal_v
  micro_cross = micro_cross_t
else
  n_velocity = fast_v
  micro_cross = micro_cross_f
end if

! Heat generation due to fissions 
power = micro_cross*a_density*n_density*n_velocity*volume*heat_per_fission

! Heat loss by conduction (forced convection)
ht_conduction = get_forced_convection()

! Heat loss by radiation 
if (reactor_temp >= 6000.0) stop "Reactor is hotter than Sun and already vaporized, stopping program!"
ht_radiation = emisivity*sb*(reactor_temp**4 - temp_equil**4)

! Calculate the energy to temperature conversion factor
fuel_heat_capacity = fuel_specific_heat*fuel_density

energy_temp_conversion = (1.0_real64/fuel_heat_capacity)*(1.0_real64/volume)

! Calculate the change in nreactor temperature
delta_reactor_temp = (power - ht_conduction - ht_radiation)*energy_temp_conversion*dt

reactor_temp = reactor_temp + delta_reactor_temp

! Output change in reactivity due to thermal feedback  
!get_feedback = delta_reactor_temp * alpha_temp
get_feedback = (reactor_temp - temp_equil) * alpha_temp

! Write results to file named T.out
write(70,*) t, n_density, reactor_temp, power, ht_conduction, ht_radiation, get_feedback

end function get_feedback



!------------------ get_forced_convection ------------------
! Function to calculate the heat loss due to forced        |
! convection. All properties and constants used were       |
! found here:                                              |
!                                                          |
!   Frank W Schmidt, Robert E. Henderson, and              |
!   Carl H. Wolgenmuth, "Introduction to Thermal Sciences" | 
!   2nd edition, p. 216-222 & 425 (1984)                   |
!   ISBN: 0-471-54939-8                                    |
!-----------------------------------------------------------
function get_forced_convection()
real(real64) :: get_forced_convection
real(real64), parameter :: pi = 4.0_real64*atan(1.0) ! Pi
real(real64), parameter :: mu = 0.008536_real64      ! Dynamic viscosity of water [g cm-1 s-1]
real(real64), parameter :: k = 0.0061056_real64      ! Thermal conductivity of water [W cm-1 K-1] 
real(real64), parameter :: w_density = .9966_real64  ! Density of water [g cm-3]
real(real64), parameter :: pr = 5.846                ! Pradntl number 
real(real64), parameter :: velocity = 100.0_real64   ! Velocity of water through core [cm s-1]
real(real64), parameter :: ch_length = pi/2.0        ! Characteristc length for one pin [cm]
real(real64), parameter :: pins = 2500               ! Number of pins in 50 x 50 array
real(real64), parameter :: area_ht = pi*(1.0)/100.0  ! Area of heat transfer of one pin [cm2]
real(real64)            :: re_lc                     ! Reynolds number at the characteristic length
real(real64)            :: nu_lam                    ! Nusselt number for laminar flow
real(real64)            :: nu_turb                   ! Nusselt number for turbulent flow
real(real64)            :: nu_avg                    ! Average Nusselt number
real(real64)            :: ht_coef                   ! Average heat transfer coefficient

! Calculate the Reynolds number at the characteristic length
re_lc = (w_density*velocity*ch_length)/mu

! Calculate the laminar flow Nusselt number
nu_lam = 0.664_real64*sqrt(re_lc)*pr**(1.0/3.0)

! Caluclate the turbulent flow Nusselt number
nu_turb = (0.037*re_lc**(0.8)*pr)/(1.0+2.443*re_lc**(-0.1)*(pr**(2.0/3.0)-1.0))

! Calculate the average Nusselt number ( 0.3 the nu_0 for pipes/tubes/cylinders )
nu_avg = 0.3 + sqrt(nu_lam**2 + nu_turb**2)

! Calculate the heat transfer coefficient 
ht_coef = nu_avg*k/ch_length

! Calculate the heat loss due to forced convection in all 2500 pins
get_forced_convection = pins*(ht_coef*area_ht*(reactor_temp - temp_equil))

end function get_forced_convection

end module feedback
