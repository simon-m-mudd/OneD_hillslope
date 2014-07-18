OneD_hillslope
==============

A numerical model of a one dimensional hillslope written in c++

This repository includes a c++ object, OneDHillslope, and various driver files for running hillslope simulations.

The object solves sediment transport folliwing the Roering et al. (1999) sediment flux law:
q(s) = -DS/( 1 - (S/S_c)^2 )
This flux is inserted into the continuity equation and then nondimensionalised using the scheme of 
Roering et al. (2008). 
The resulting governing equation is then solved implicitly. 
