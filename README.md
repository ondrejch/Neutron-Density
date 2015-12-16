# Neutron-Density
A program that takes in values for reactivty and outputs the neutrion density at given time intervals.

The file main.f90 acts as the main program file that will handle the input and pass it to the subroutine neuden() 
neutrondens.f90. 

Example inputs:

input1 : fast/ constant reactivity of 0.0022/ t = 0 to 10

input2 : fast/ ramp reactivity of 0.0044t/ t = 0 to 1 

input3 : thermal/ constant reactivity of -0.00375/ t = 0 to 10

input_zigzag : thermal/ zigzag reactivity described on page 9 of article/ t = 0 to 10

Link to the article by Xue Yang and Tatjana Jevremovic:
http://www.doiserbia.nb.rs/img/doi/1451-3994/2009/1451-39940901003Y.pdf
