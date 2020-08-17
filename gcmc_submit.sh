#!/bin/bash

# use current working directory for input and output
# default is to use the users home directory
# print date and time
date
~/julia-1.4.2/bin/julia -p 4 cof_isotherm_sim.jl $xtal
