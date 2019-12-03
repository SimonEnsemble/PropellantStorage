#!/bin/bash

# use current working directory for input and output
# default is to use the users home directory
#$ -cwd

# name this job
#$ -N cof_103

#$ -pe thread 4 # use 4 threads/cores

# send stdout and stderror to this file
#$ -o cof_103.o
#$ -e cof_103.e

# select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.
##$ -q mime5

# print date and time
date
~/julia-1.1.0/bin/julia cof_103_isotherm_sim.jl 
