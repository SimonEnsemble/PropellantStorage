using PorousMaterials

# read in crystal structure name from command line arguments
if length(ARGS) != 1
    error("pass the crystal structure name as a command line argument, such as:
    julia cof_isotherm_sim.jl COF-102.cif")
end
crystal = ARGS[1]
println("running mol sim in ", crystal)

frame = Framework(crystal)
strip_numbers_from_atom_labels!(frame)

temp =  298.0 # K
forcefield = LJForceField("Dreiding.csv", cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")
# forcefield = LJForceField("UFF.csv", cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")

xe = Molecule("Xe")

pmin  = -2  # log10 of minimum pressure in bar
pmax  = 300 # actual value of maximum pressure in bar
nstep = 15  # number of steps from pmin to pmax  
pressures = 10 .^ range(pmin, stop=log10(pmax), length=nstep) # bar

n_sample_cycles = 25000
n_burn_cycles = 25000

adsorption_data = adsorption_isotherm(frame, xe, temp, pressures, forcefield,
		                              n_burn_cycles=n_burn_cycles, 
                                      n_sample_cycles=n_sample_cycles,
                                      eos=:PengRobinson)
