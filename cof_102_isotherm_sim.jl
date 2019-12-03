using PorousMaterials

crystal = "102.cif"
frame = Framework(crystal)
strip_numbers_from_atom_labels!(frame)

temp =  298.0
ljff = LJForceField("UFF.csv")

xe = Molecule("Xe")

pressures = 3 * 10 .^ range(-2, stop=2, length= 100) # bar

n_sample_cycles = 10000
n_burn_cycles = 5000

adsorption_data = adsorption_isotherm(frame, xe, temp, pressures,
                    ljff, n_burn_cycles = n_burn_cycles, n_sample_cycles = n_sample_cycles,
                    eos = :PengRobinson)

