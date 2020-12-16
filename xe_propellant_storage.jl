### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 6091fca2-0b62-11eb-0351-93554cbc8d52
using PorousMaterials, PyPlot, DataFrames, CSV, Interpolations, Optim, Printf, PyCall, Roots

# ╔═╡ 54e4a2c4-0b62-11eb-3e1d-017a70ae43d7
md"
# Xe propellant storage
"

# ╔═╡ 1cffb8d4-18a3-11eb-087f-83393730c57b
adjustText = pyimport("adjustText")

# ╔═╡ 68ff0da6-0b66-11eb-1911-bf3faa3ceb7a
if ! isdir("figz")
	mkdir("figz")
end

# ╔═╡ 422c02b2-1802-11eb-31a7-c3671345f577
begin
	PyPlot.matplotlib.style.use("seaborn-notebook")

	rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcParams["font.size"] = 14

	const figsize = (5, 4)
end

# ╔═╡ 0c7cac6e-1802-11eb-122a-0180ae055dbc
function despine()
	ax = gca()
	ax.spines["top"].set_visible(false)
	ax.spines["right"].set_visible(false)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
end

# ╔═╡ 6f90c580-0b62-11eb-092a-3f99accec091
md"### define constants
1. Tam,  W.;  Jackson,  A.;  Nishida,  E.;  Ka-sai, Y.; Tsujihata, A.; Kajiwara, K. Designand  manufacture  of  the  ETS  VIII  xenontank.36thAIAA/ASME/SAE/ASEEJoint  Propulsion  Conference  and  Exhibit.2000; p 3677.(https://doi.org/10.2514/6.2000-3677)

2. Welsch,  G.;  Boyer,  R.;  Collings,  E.Materials properties handbook: titanium alloys;ASM international, 1993. 

3. Niinomi,   M.   Mechanical   properties   of biomedical  titanium  alloys. Materials Science and Engineering: A1998,243, 231–236.

4. propellant storage considerations [link](https://erps.spacegrant.org/uploads/images/images/iepc_articledownload_1988-2007/1991index/IEPC1991-107.pdf)

"

# ╔═╡ 6b759e8a-0b62-11eb-059d-8325bf28a129
begin
	# density of walls of pressure vessel
	const ρ_t = 4428.785      # kg/m³ (converted from 0.16 lb/in³ listed in [2])
	
	# yield strength of walls of pressure vessel
	const σ_y = 8250.0        # bar   (lower-limit of reported values [3])
	
	# temperature
	const temperature = 298.0 # K
	
	# safety factor
	const β = 0.5             # unitless
	
	# Universal Gas Constant:
	const R = 8.3144598e-5; # m³-bar/(K-mol)
	
	const xe_molar_mass = 131.293 / 1000    # kg Xe/ mol Xe 
	# common amount of xenon gas required to bring into space [1].
	const mass_desired_xe_propellant = 100.0 # kg Xe
	const mol_desired_xe_propellant = mass_desired_xe_propellant / xe_molar_mass # mol Xe
end

# ╔═╡ e4f7db60-0b65-11eb-1dda-37f708ccc230
md"## bulk Xe data"

# ╔═╡ ed862216-0b65-11eb-188a-ab1107ab8e24
# critical pressure
const Pc = 58.420 # bar, from NIST

# ╔═╡ fdd5dc94-0b65-11eb-1907-47bbff769958
begin
	df_xe = vcat(CSV.read(joinpath("data", "NIST_data", "low_pressure_xenon_NIST_data.txt"), DataFrame)[2:end, 2:3],
	             CSV.read(joinpath("data", "NIST_data", "xenon_NIST_data.txt"), DataFrame)[:, 2:3])
	df_xe[:, Symbol("Density (mol/m³)")] = df_xe[:, Symbol("Density (mol/l)")] * 1000 # convert L to m³
	sort!(df_xe, Symbol("Pressure (bar)"))
	
	first(df_xe, 3)
end

# ╔═╡ a4a6f3be-0f65-11eb-0c80-53046792dd1f
md"
numerical implementation of $\rho_{Xe}(P)$
"

# ╔═╡ 2c25ab68-0b66-11eb-0f41-a91ffdf1e371
ρ_xe = LinearInterpolation(df_xe[:, Symbol("Pressure (bar)")], 
	                       df_xe[:, Symbol("Density (mol/m³)")]
	                       )
# example use: 
#   ρ_xe(59.0)

# ╔═╡ 5bbe5418-0b66-11eb-1d57-798a7b5b78de
begin
	function viz_bulk_ρ()
		P_range_bulk = range(0.0, 350.0, length=350) # bar
		# ideal gas density
		ρ_ideal_gas = P_range_bulk / (R * temperature) # mol / m³
		
		figure(figsize=figsize)
		plot(P_range_bulk, ρ_xe.(P_range_bulk), color="C0", 
			label="xenon (NIST data)", lw=3, zorder=4)
		vlines(Pc, 0.0, ρ_xe(Pc), linestyle="-.", 
			color="gray", lw=2, label="critical pressure")
		plot(P_range_bulk, ρ_ideal_gas, color="C2", label="ideal gas", 
			lw=3, linestyle="--")
		xlabel(L"pressure, $P$ [bar]")
		ylabel(L"density of bulk xenon, $\rho_{Xe}$ [mol/m$^3$]")
		ylim([0.0, 17500*1.025])
		xlim([0.0, 351])
		legend()
		tight_layout()
		savefig("figz/xenon_gas_density.pdf", bbox_inches="tight")
		gcf()
	end
	
	viz_bulk_ρ()
end

# ╔═╡ 29239f5e-0b63-11eb-1b90-67892ca5c215
md"## materials data"

# ╔═╡ 7cc89760-0b64-11eb-2360-c97a464d908f
const xtal_names = ["SBMOF-1", "CC3", "Ni-MOF-74", "HKUST-1", "SBMOF-2", 
	                "Co-formate", "MOF-505", "Activated-Carbon", "NiPyC2", 
	                "COF-103 (simulated)"]

# ╔═╡ b5cef332-0f59-11eb-22c4-cf85220a7af5
md"##### get xtal densities"

# ╔═╡ be69b9dc-0f59-11eb-2361-69b0f7e1ccf3
begin
	xtal_to_ρ = Dict{String, Float64}()
	
	# find xtal density using PorousMaterials.jl
	for xtal_name in xtal_names
	    if xtal_name == "Activated-Carbon"
			# use density from pamphlet
			xtal_to_ρ[xtal_name] = 500.0 # kg/m³
		else 
			# use crystal density
			xtal_filename = xtal_name
	        if xtal_name in ["NiPyC2", "COF-103 (simulated)"]
	            xtal_filename *= ".cif"
	        else    
	            xtal_filename *= ".cssr"
	        end
			xtal = Crystal(xtal_filename)
			strip_numbers_from_atom_labels!(xtal)
	        xtal_to_ρ[xtal_name] = crystal_density(xtal) # kg/m³
	    end
	end
	
	xtal_to_ρ # kg/m³	
end

# ╔═╡ b853e33e-0f35-11eb-2027-9d77b137546e
md"
##### Fit Langmuir adsorption isotherm models to adsorption data in the materials

well, first, read them in!
"

# ╔═╡ eb9159f2-0f35-11eb-156e-d3762896ec0d
begin
	const PATH_TO_EXP_ISOTHERMS = joinpath("data", "exptl_xe_isotherms")
	const PATH_TO_SIM_ISOTHERMS = joinpath("data", "sim_xe_isotherms")
end

# ╔═╡ 0be36a3a-0f36-11eb-0d80-d11b529da9ad
begin
	xe_isotherms = Dict{String, DataFrame}()
	
	for xtal_name in xtal_names
	    isotherm_filename = ""
		if xtal_name == "COF-103 (simulated)"
			isotherm_filename = joinpath(PATH_TO_SIM_ISOTHERMS, xtal_name * ".csv")
		else
			isotherm_filename = joinpath(PATH_TO_EXP_ISOTHERMS, xtal_name, "Xe.csv")
		end
		
		xe_isotherms[xtal_name] = CSV.read(isotherm_filename, DataFrame; copycols=true)
	    # remove hysteresis branch manually
	    if xtal_name == "FMOF-Cu" 
	        deleterows!(xe_isotherms[xtal_name], 15:24)
	    elseif xtal_name == "SBMOF-2"
	        deleterows!(xe_isotherms[xtal_name], 30:39)
	    elseif xtal_name == "COF-103 (simulated)"
	        xe_isotherms[xtal_name] = xe_isotherms[xtal_name][:, [Symbol("pressure (bar)"), Symbol("⟨N⟩ (mmol/g)")]]
	    end
	end
	
	xe_isotherms["SBMOF-1"]
end

# ╔═╡ eed22218-0f5b-11eb-15a1-a75120aaf145
md"convert units."

# ╔═╡ 40efd702-0f5b-11eb-2970-c1ee6ac4eedf
begin
	# Define what the desired common units are for the data
	common_pressure_units = Symbol("Pressure (bar)")
	common_loading_units = Symbol("Loading (mol/kg)")
	
	# Define a dictionary with conversion factors.
	pressure_conversion = Dict{Symbol, Float64}()
	loading_conversion = Dict{Symbol, Float64}()
	
	# pressure conversions to bar
	pressure_conversion[Symbol("pressure (bar)")] = 1.0 # (1 bar) / (1 bar)
	pressure_conversion[Symbol("P(bar)")]  = 1.0 # (1 bar) / (1 bar)
	# pressure_conversion[Symbol("fugacity (bar)")] = 1.0 # (1 bar) / (1 bar)
	pressure_conversion[Symbol("P(mbar)")] = 1 / 1000 # (1 bar) / (1000 mbar)
	pressure_conversion[Symbol("P(kPa)")]  = 1 / 100 # (1 bar) / (100 kPa)
	pressure_conversion[Symbol("P(torr)")] = 1 / 750.062 # (1 bar) / (750.062 torr)
	pressure_conversion[Symbol("P(atm)")]  = 1 / 0.986923 # (1 bar) / (0.986923 atm)
	pressure_conversion[Symbol("P(mmHg)")] = 1 / 750.062 # (1 bar) / (750.026 mmHg)
	
	# loading conversions to mol/kg
	# these conversion factors will put the quantity into mmol/g
	loading_conversion[Symbol("L(mmol/g)")]    = 1.0 # (1 mol/kg) / (1 mmol/g)
	loading_conversion[Symbol("⟨N⟩ (mmol/g)")] = 1.0 # (1 mol/kg) / (1 mmol/g)
	loading_conversion[Symbol("L(ccSTP/g)")]   = 1 / 22.4 # (cc STP /g) (1000 g /kg) (1 mol/ 22.4 L STP) (1 L / 1000 cc)
	loading_conversion[Symbol("L(cm3STP/g)")]  = 1 / 22.4
	# [(% mass) / 100 g Xe / g MOF](1 mol / MW_Xe g) (1000 g /1 kg)
	loading_conversion[:PercentMass] = 1000.0 / 131.1 / 100.0
	# this one is an exception where xtal density not needed
	# loading_conversion[Symbol("L(mol/L)")] = 1000.0 # (mol / L)(1000 L / m3)
end

# ╔═╡ 0fdab85c-0f5c-11eb-2cec-69e1a95c6e45
for xtal_name in xtal_names
    # loop over columns in the adsorption isotherm DataFrame
    for col_name in names(xe_isotherms[xtal_name])
        # convert pressure units to bar
        if col_name in keys(pressure_conversion)
            xe_isotherms[xtal_name][!, common_pressure_units] = xe_isotherms[xtal_name][:, col_name] * pressure_conversion[col_name]
        # convert loading units to mol/m3
        elseif col_name in keys(loading_conversion)
            # if col_name == Symbol("L(mol/L)")
            #     xe_isotherms[xtal_name][!, common_loading_units] = xe_isotherms[xtal_name][!, col_name] * loading_conversion[col_name]
            # else 
			xe_isotherms[xtal_name][!, common_loading_units] = xe_isotherms[xtal_name][:, col_name] * loading_conversion[col_name]
            # end 
		end
    end
	
	@assert common_loading_units in names(xe_isotherms[xtal_name])
	@assert common_pressure_units in names(xe_isotherms[xtal_name])
end

# ╔═╡ fb3bb86c-0f5c-11eb-0438-57bba5356391
md"finally, fit Langmuir models."

# ╔═╡ 06ad0a18-0f5d-11eb-32bb-874fc19bb525
begin
	xtal_to_K = Dict{String, Float64}()
	xtal_to_M = Dict{String, Float64}()
	for xtal_name in xtal_names
	    lang_params = fit_adsorption_isotherm(xe_isotherms[xtal_name],
			                                  common_pressure_units,
	                                          common_loading_units, 
			                                  :langmuir)
		xtal_to_K[xtal_name] = lang_params["K"]
		xtal_to_M[xtal_name] = lang_params["M"]
	end
end

# ╔═╡ d336fcca-0f5f-11eb-3cb2-7f5c3f040f7b
md"for plotting..."

# ╔═╡ 331c7dda-0dcc-11eb-1658-d7c4d23b788e
begin
	xtal_to_label = Dict(xtal_name => xtal_name for xtal_name in xtal_names)
	xtal_to_label["NiPyC2"] = L"Ni(PyC)$_2$"
	xtal_to_label["Activated-Carbon"] = "activated carbon"
	xtal_to_label["Co-formate"] = L"Co$_3$(HCOO)$_6$"
	xtal_to_label
end

# ╔═╡ feaf132e-0b63-11eb-28ff-7f327d815001
begin
	# for plotting
	_markers = ["o", ">", "<", "*", "H", "^", "v",  "d", "s", "8"]
	
	const xtal_to_marker = Dict(zip(xtal_names, _markers[1:length(xtal_names)]))
		
	sborn = pyimport("seaborn")
	colorz = sborn.color_palette("husl", length(xtal_names))
	const xtal_to_color =  Dict(zip(xtal_names, colorz))
end

# ╔═╡ fa9fe136-0b64-11eb-3718-3527c3599f56
md"
our $\rho_{Xe}^{ads}(P)$ [mol/m³]
"

# ╔═╡ ef807360-0b64-11eb-0f1b-cdd324b8b1bb
# P [=] bar
function ρ_xe_ads(P::Float64, xtal_name::String)
    M = xtal_to_M[xtal_name]
    K = xtal_to_K[xtal_name]
	ρ = xtal_to_ρ[xtal_name]
    return M * ρ * K * P / (1 + K * P) # mol/m³
end

# ╔═╡ e7a49568-1893-11eb-321c-69b330620a80
function plot_langmuir_fits(;gravimetric::Bool=false)    
	figure(figsize=(8, 5))
	xlabel(L"pressure, $P$ [bar]")
	if gravimetric
    	ylabel(L"adsorbed Xe density, $\rho_{Xe}$ [mol/kg]")
	else
		ylabel(L"adsorbed Xe density, $\rho_{Xe}$ [mol/m$^3$]")
	end

	# Langmuir fits
    Ps = range(0.0, stop=1.2, length=500)
	for xtal_name in xtal_names	
		if xtal_name == "COF-103 (simulated)"
			continue
		end
		
		if gravimetric		
			plot(Ps, ρ_xe_ads.(Ps, xtal_name) / xtal_to_ρ[xtal_name], color=xtal_to_color[xtal_name])

			scatter(xe_isotherms[xtal_name][:, common_pressure_units], 
					xe_isotherms[xtal_name][:, common_loading_units],
					color=xtal_to_color[xtal_name],
					edgecolor="k", zorder=400, label=xtal_to_label[xtal_name])
		else	
			plot(Ps, ρ_xe_ads.(Ps, xtal_name), color=xtal_to_color[xtal_name])

			scatter(xe_isotherms[xtal_name][:, common_pressure_units], 
					xe_isotherms[xtal_name][:, common_loading_units] * xtal_to_ρ[xtal_name],
					color=xtal_to_color[xtal_name],
					edgecolor="k", zorder=400, label=xtal_to_label[xtal_name])
		end
	end
	if ! gravimetric
	    # bulk gas density
   		plot(Ps, ρ_xe.(Ps), label="bulk gas", linestyle="--", color="gray", lw=1)
	end
	
    legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
	ylim(ymin=0.0)
	xlim([0, 1.05])
    tight_layout()
    if gravimetric
    	savefig("figz/langmuir_fits_gravimetric.pdf", format="pdf", bbox_inches="tight")
	else
		savefig("figz/langmuir_fits_volumetric.pdf", format="pdf", bbox_inches="tight")
	end
	gcf()
end

# ╔═╡ f5711acc-1893-11eb-1937-337d2b173321
plot_langmuir_fits()

# ╔═╡ 2d85d0dc-0b66-11eb-1a99-7fee2143ef57
md"our $\rho_{Xe}(P)$ [mol/m³]"

# ╔═╡ 83c86fb4-3a7c-11eb-281d-0b0fac556a7c
plot_langmuir_fits(gravimetric=true)

# ╔═╡ aa5c1ad4-0b67-11eb-37b5-9746aa76c7f8
md"
## bulk Xe storage
"

# ╔═╡ 218243b0-0dc8-11eb-01ed-258156f8aa9d
Ps = range(20.0, 125.0, length=300)

# ╔═╡ 9c3c5810-0b67-11eb-0141-8bfa134546a6
begin
	# radius of vessel
	# input:
	#       P [=] bar
	# output:
	#       r [=] meters
	function r(P::Float64) # bulk Xe storage
	    return (3.0 * mol_desired_xe_propellant / (4.0 * π * ρ_xe(P))) ^ (1/3)
	end
	
	# thickness of vessel walls
	# input:
	#       P [=] bar
	# output:
	#       t [=] meters
	function t(P::Float64) # bulk Xe storage
	    return P * r(P) / (2.0 * β * σ_y)
	end
	
	# mass of pressure vessel (tank)
	# input:
	#       P [=] bar
	# output:
	#       mₜ [=] kg
	function mₜ(P::Float64)  # bulk Xe storage
	    return 4.0 * π * r(P) ^ 2.0 * t(P) * ρ_t
	end
	
	function tankage_fract(P::Float64)
		return mₜ(P) / mass_desired_xe_propellant
	end
	
	# optimum pressure and mass of vessel
	# input:
	#       none
	# output:
	#       .minimizer [=] bar
	#       .minimum [=] kg  
	function optimum_storage() # bulk Xe storage
	    return optimize(P -> mₜ(P[1]), 0.00001, 100.0)
	end
end

# ╔═╡ d90fc4f0-0b66-11eb-33da-372fc8695a3c
md"
#### adsorbed Xe storage

if the function takes in a crystal name, then it corresponds to adsorbed xenon storage; otherwise, it corresponds to bulk storage.

###### adsorbed storage
"

# ╔═╡ e1dbf234-0b66-11eb-244a-33506ae7db4a
begin
	# radius of vessel
	# input:
	#       P [=] bar
	# output:
	#       r [=] meters
	function r(P::Float64, xtal::String) # ads Xe storage
	    return (3.0 * mol_desired_xe_propellant / (4.0 * π * ρ_xe_ads(P, xtal))) ^ (1/3)
	end
	
	# thickness of vessel walls
	# input:
	#       P [=] bar
	# output:
	#       t [=] meters
	function t(P::Float64, xtal::String) # ads Xe storage
	    return P * r(P, xtal) / (2.0 * β * σ_y)
	end
	
	# mass of vessel material
	# input:
	#       P [=] bar
	# output:
	#       mₜ [=] kg
	function mₜ(P::Float64, xtal::String) # ads Xe storage
	    return 4.0 * π * r(P, xtal) ^ 2.0 * t(P, xtal) * ρ_t
	end
	
	# mass of adsorbent
	# input:
	#       P [=] bar
	# output:
	#       mₐ [=] kg
	function mₐ(P::Float64, xtal::String) # ads Xe storage (of course)
	    return xtal_to_ρ[xtal] * mol_desired_xe_propellant / ρ_xe_ads(P, xtal)
	end
	
	function tankage_fract(P::Float64, xtal::String)
		return (mₜ(P, xtal) + mₐ(P, xtal))/ mass_desired_xe_propellant
	end
	
	# optimum pressure and mass of vessel
	# input:
	#       none
	# output:
	#       .minimizer [=] opt. pressure, bar
	#       .minimum   [=] opt. mass of storage media, kg  
	function optimum_storage(xtal::String) # ads Xe storage
	    # use Optim.jl to find the minimum
	    res = optimize(P -> (mₜ(P, xtal) + mₐ(P, xtal)), 0.00001, 100.0)
	    @assert res.converged "Optimization not successful."
	    return res
	end
end

# ╔═╡ 2f55d9e6-0b68-11eb-271a-b502d0bc5a33
begin
	bulk_opt = Dict()
	bulk_opt["P (bar)"] = optimum_storage().minimizer
	bulk_opt["tf"] = optimum_storage().minimum / mass_desired_xe_propellant
	bulk_opt["t [m]"] = t(bulk_opt["P (bar)"])
	bulk_opt["r [m]"] = r(bulk_opt["P (bar)"])
	bulk_opt
end

# ╔═╡ 33ed9f0e-0dc8-11eb-1b06-89178c499937
begin
	function viz_bulk_storage()
	    figure(figsize=(5, 10))
		ax1 = plt.subplot2grid((6,1), (0,0))
		ax2 = plt.subplot2grid((6,1), (1,0))
		ax3 = plt.subplot2grid((6,1), (2,0), rowspan=2)
		
		ax1.set_title("bulk Xe storage")

		# tank radius
		ax1.plot(Ps, r.(Ps), color="C0", lw=3)
		ax1.set_ylabel(L"radius, $r$ [m]")
		ax1.set_ylim(ymin=0.0)
		ax1.set_xticks([])

		# tank thickness
		ax2.plot(Ps, t.(Ps) * 1000, color="C1", lw=3)
		ax2.set_ylabel(L"thickness, $t$ [mm]")
		ax2.set_ylim(ymin=0.0)
		ax2.set_xticks([])
		# # density of adsorbed xe
		# ax3.plot(Ps, ρ_xe.(Ps), color="C2", clip_on=false)
		# ax3.set_ylabel("ρₓₑ [mol/m\$^3\$]")
		
		# tankage fraction
		ax3.plot(Ps, mₜ.(Ps) / mass_desired_xe_propellant, color="C3", 
			linestyle="-", label="tank walls", lw=3)
		ax3.set_ylabel(L"tankage fraction, $m_v/(n_{Xe}w_{Xe})$")
		ax3.set_xlabel(L"pressure, $P$ [bar]")
		ax3.set_ylim(ymin=0.0)
		ax3.annotate(#L"$(P_{opt}=$" * 
	         @sprintf("(%.1f bar, %.2f)", bulk_opt["P (bar)"], bulk_opt["tf"]),
	         xy=(bulk_opt["P (bar)"], bulk_opt["tf"]),
		      xytext=(100, 0.15), ha="center",
			 fontsize=12,
	         arrowprops=Dict(:arrowstyle=>"-|>", 
			                 :connectionstyle=>"arc3,rad=-0.2", 
			                 :fc=>"k", 
			                 :ec=>"k", 
			                 :linewidth=>2), 
	         bbox=Dict("boxstyle"=>"round", "fc"=>"C6", "alpha"=>0.3),
		)
		ax3.scatter([bulk_opt["P (bar)"]], [bulk_opt["tf"]], 
			marker="x", color="black", label="optimum", zorder=500, s=50)

		for ax in [ax1, ax2, ax3]
			ax.set_xlim([20, 125])
	        # ax.axhline(y=0, color="0.6", zorder=1, clip_on=false)
		    ax.axvline(x=bulk_opt["P (bar)"], 
				       linestyle="--", color="0.3", zorder=1, clip_on=false)
		end
		
		

	    savefig("figz/bulk_storage_opt.pdf", bbox_inches="tight")
		gcf()
	end
	
	viz_bulk_storage()
end

# ╔═╡ 11185356-0b67-11eb-1bb3-e15f94234be3
md"
test analytical solutions from paper
"

# ╔═╡ 1721c25a-0b67-11eb-1b47-9bad136bd63a
begin
	function P_opt_analytical(xtal::String) 
	    return sqrt(2 * σ_y * xtal_to_ρ[xtal] * β / (3 * xtal_to_K[xtal] * ρ_t))
	end
	
	function tankage_fraction_opt_analytical(xtal::String)
		return 1 / (xtal_to_M[xtal] * xtal_to_ρ[xtal] * xe_molar_mass) * (sqrt(xtal_to_ρ[xtal]) + sqrt(3 * ρ_t / (2 * σ_y * β * xtal_to_K[xtal]))) ^ 2
	end
	
	function test_analytical_eqns()
		for xtal in xtal_names
			# numerical
			res = optimum_storage(xtal)
			# analytical
			P_opt = P_opt_analytical(xtal)
			tf_opt = tankage_fraction_opt_analytical(xtal)

			@assert isapprox(res.minimizer, P_opt, atol=1e-6)
			@assert res.minimum / mass_desired_xe_propellant ≈ tf_opt
			@assert tankage_fract(P_opt, xtal) ≈ tf_opt
		end
		return "tests pass"
	end
	
	test_analytical_eqns()
end

# ╔═╡ f4cb6ea8-0da6-11eb-120e-8d67d8ca58b7
md"### adsorbed Xe storage
"

# ╔═╡ 7fe2a354-0f33-11eb-2f83-e3b66bbdb333
function scatter_kwargs(xtal::String)
	kws = Dict(:label     => xtal_to_label[xtal],
			   :marker    => xtal_to_marker[xtal],
			   :s         => 75,
			   :clip_on   => false,
			   :color     => xtal_to_color[xtal],
			   :edgecolor => xtal_to_color[xtal],
	           :lw        => 3)
	# make COF-103 hollow to emphasize it's simulated.
	if xtal == "COF-103 (simulated)"
		kws[:color] = "None"
	end
	return kws
end

# ╔═╡ 018c2b64-0da7-11eb-100a-4d1ee4e512f3
begin
	ads_opt = Dict()
	
	for xtal in xtal_names
	    opt = optimum_storage(xtal)
		
	    ads_opt[xtal] = Dict()
	    ads_opt[xtal]["P [bar]"] = opt.minimizer
	    ads_opt[xtal]["tf"] = opt.minimum / mass_desired_xe_propellant
	    ads_opt[xtal]["t [m]"] = t(opt.minimizer, xtal)
	    ads_opt[xtal]["r [m]"] = r(opt.minimizer, xtal)
	    ads_opt[xtal]["m_t [kg]"] = mₜ(opt.minimizer, xtal)
	    ads_opt[xtal]["m_a [kg]"] = mₐ(opt.minimizer, xtal)
	    ads_opt[xtal]["m [kg]"] = mₐ(opt.minimizer, xtal) + mₜ(opt.minimizer, xtal)
		
		@assert ads_opt[xtal]["tf"] ≈ tankage_fraction_opt_analytical(xtal)
		@assert isapprox(ads_opt[xtal]["P [bar]"], P_opt_analytical(xtal), atol=0.001)
	end
end

# ╔═╡ d182f874-18bb-11eb-28b8-ad20c019015e
md"
summarize
"

# ╔═╡ c5a0c9e2-18a1-11eb-17a4-8572d491905e
md"should this be volume at same pressure?
or volume at optimal pressure?
both seem a valuable comparison.
"

# ╔═╡ 1015e5fa-0dcd-11eb-3352-09bc6056a02f
begin
	function performance_plot()
		figure(figsize=figsize)
		xlabel(L"optimal storage pressure, $P_{opt}$ [bar]")
		ylabel("optimal tankage fraction")
		texts = []
		for xtal_name in xtal_names
			p_opt = ads_opt[xtal_name]["P [bar]"]
			tf_opt = ads_opt[xtal_name]["tf"]
			scatter(p_opt, tf_opt ; scatter_kwargs(xtal_name)...)
			
    		push!(texts,
        		  annotate(xtal_to_label[xtal_name], (p_opt, tf_opt), fontsize=10)
        	)
		end
		
		p_opt_b = bulk_opt["P (bar)"]
		tf_opt_b = bulk_opt["tf"]
		scatter([p_opt_b], [tf_opt_b], 
				marker="x", color="black", label="bulk storage")
		push!(texts,
        	  annotate("bulk storage", (p_opt_b, tf_opt_b), fontsize=10)
       	)
		adjustText.adjust_text(texts, autoalign="xy")#, force_points=(1.0, 2.0))
		
		ylim(ymin=0.0)
		xlim(xmin=0.0)
		tight_layout()
		# legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
		savefig("figz/tf_vs_P_opt.pdf")
		gcf()
	end
	
	performance_plot()
end

# ╔═╡ 95054c18-2b93-11eb-0d17-3f8166d00547
ads_opt["SBMOF-1"]["P [bar]"]

# ╔═╡ af8c02fc-2b93-11eb-3e3e-8f4695ff5bbf
ads_opt["SBMOF-1"]["tf"]

# ╔═╡ ce3da226-18a3-11eb-25df-572b24f5d3f0
close("all")

# ╔═╡ a8b0c9c4-0dcd-11eb-2924-91b2875da4ea
begin
	ids = 1:length(xtal_names)
	ids_sorted = sortperm([ads_opt[xtal]["tf"] for xtal in xtal_names])                      
	figure()
	xlabel("material")
	ylabel("mass [kg]")
	
	m_xtal = [ads_opt[xtal]["m_a [kg]"] for xtal in xtal_names]
	m_tank_walls = [ads_opt[xtal]["m_t [kg]"] for xtal in xtal_names]
	
	bar(ids, m_xtal[ids_sorted], label=L"adsorbent, $m_{ads}$")
	bar(ids, m_tank_walls[ids_sorted], bottom=m_xtal[ids_sorted],
		label=L"vessel walls, $m_v$", color="C1")
	
	legend()
	xticks(ids, [xtal_to_label[xtal] for xtal in xtal_names[ids_sorted]], rotation="vertical")
	tight_layout()
	savefig("figz/m_t_m_a.pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ e4ed9b6c-18bb-11eb-2796-8fe1689d5e35
begin
	function viz_mt()
		ids = 1:length(xtal_names)
		colors = [xtal_to_color[xtal] for xtal in xtal_names][ids_sorted]
		
		fig, axs = plt.subplots(4, 1, figsize=(5, 7), sharex=true)
		xlabel("material")
		xticks(ids, 
			   [xtal_to_label[xtal_name] for xtal_name in xtal_names][ids_sorted], 
			   rotation="vertical")
		
		# m_v
		axs[1].set_ylabel(L"$m_v(P_{opt})$ [kg]")
		m_ts = [ads_opt[xtal]["m_t [kg]"] for xtal in xtal_names][ids_sorted]
		axs[1].bar(ids, m_ts, color=colors)
		m_t_b = mₜ(bulk_opt["P (bar)"])
		axs[1].axhline(y=m_t_b, linestyle="--", color="gray")
		
		# m_a
		axs[2].set_ylabel(L"$m_{ads}(P_{opt})$ [kg]")
		m_as = [ads_opt[xtal]["m_a [kg]"] for xtal in xtal_names][ids_sorted]
		axs[2].bar(ids, m_as, color=colors)
		
		# r
		axs[3].set_ylabel(L"$r(P_{opt})$ [m]")
		rs = [ads_opt[xtal]["r [m]"] for xtal in xtal_names][ids_sorted]
		axs[3].bar(ids, rs, color=colors)
		r_b = r(bulk_opt["P (bar)"])
		axs[3].axhline(y=r_b, linestyle="--", color="gray")
		
		# t
		axs[4].set_ylabel(L"$t(P_{opt})$ [mm]")
		ts = [ads_opt[xtal]["t [m]"] for xtal in xtal_names][ids_sorted] * 1000
		axs[4].bar(ids, ts, color=colors)
		t_b = t(bulk_opt["P (bar)"]) * 1000
		axs[4].axhline(y=t_b, linestyle="--", color="gray")
		
		tight_layout()
		savefig("figz/ads_summary.pdf", bbox_inches="tight")
		gcf()
	end
	
	viz_mt()
end

# ╔═╡ fb16b62c-0f32-11eb-05fd-e72d057be910
md"
how does tankage fraction relate to materials properties?
* crystal density, $\rho_m$
* Langmuir $K$
* Langmuir $M$
"

# ╔═╡ 37cf428e-0dce-11eb-154b-25742d645575
begin
	figure(figsize=figsize)
	xlabel(L"density of material, $\rho_{ads}$ [kg/m$^3$]")
	ylabel("optimal tankage fraction")
	texts = []
	for xtal_name in xtal_names
		ρ = xtal_to_ρ[xtal_name]
		tf_opt = ads_opt[xtal_name]["tf"]
	    scatter(ρ, tf_opt; scatter_kwargs(xtal_name)...)
		push!(texts,
        	  annotate(xtal_to_label[xtal_name], (ρ, tf_opt), fontsize=10)
        )
	end
	
	axhline(y=[bulk_opt["tf"]], linestyle="--", color="gray", label="bulk storage")
	adjustText.adjust_text(texts, autoalign="xy")#, force_points=(1.0, 2.0))

	ylim(ymin=0.0)
	xlim(xmin=0.0)
	tight_layout()
	# legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
	savefig("figz/tf_vs_crystal_density.pdf", format="pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ 105f4f12-0f33-11eb-161c-a1d26651a48c
begin
	figure(figsize=figsize)
	xlabel(L"Langmuir $K$ [1/bar]")
	ylabel("optimal tankage fraction")
	for xtal in xtal_names
	    scatter(xtal_to_K[xtal], ads_opt[xtal]["tf"],
	   		    ; scatter_kwargs(xtal)...)
	end
	
	axhline(y=[bulk_opt["tf"]], linestyle="--", color="black", label="bulk storage")
	
	ylim(ymin=0.0)
	xlim(xmin=0.0)
	legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
	savefig("figz/tf_vs_K.pdf", format="pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ 11b18cca-0f33-11eb-2556-f946a9ac7d0b
begin
	figure(figsize=figsize)
	xlabel(L"Langmuir $M$ [mol/kg]")
	ylabel("optimal tankage fraction")
	
	for xtal in xtal_names
	    scatter(xtal_to_M[xtal], ads_opt[xtal]["tf"],
	   		    ; scatter_kwargs(xtal)...)
	end
	
	axhline(y=[bulk_opt["tf"]], linestyle="--", color="black", label="bulk storage")
	
	ylim(ymin=0.0)
	xlim(xmin=0.0)
	legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
	savefig("figures/tf_vs_M.pdf",                    
		    bbox_inches="tight")
	gcf()
end

# ╔═╡ e9c0dacc-189e-11eb-0a0b-49ea73af6ef7
function materials_space_viz()
# 	fig = plt.figure()
# 	ax = fig.add_subplot(111, projection="3d")
# 	title("material space")

# 	for xtal_name in xtal_names
# 		x = [xtal_to_K[xtal_name], xtal_to_M[xtal_name], xtal_to_ρ[xtal_name]]
		
# 		scatter3D(x...,
# 			label=xtal_to_label[xtal_name], 
# 			marker=xtal_to_marker[xtal_name], 
# 			s=75, color=xtal_to_color[xtal_name])
		
# 		for k = 1:3
# 			x_projected = copy(x)
# 			if k == 2
# 				x_projected[k] = 40
# 			else
# 				x_projected[k] = 0
# 			end
# 			plot(
# 				[x_projected[1], x[1]], 
# 				[x_projected[2], x[2]], 
# 				[x_projected[3], x[3]], 
# 				color=xtal_to_color[xtal_name], alpha=0.7
# 				)
# 		end
# 	end

# 	xlabel(L"Langmuir $K$ [bar$^{-1}$]", labelpad=10, fontsize=12)
# 	ylabel(L"Langmuir $M$ [mol/kg]", labelpad=10, fontsize=12)
# 	zlabel(L"density, ρ$_{ads}$ [kg/m³]", labelpad=10, fontsize=12)

# 	tick_params(labelsize=10)
# 	ylim(ymin=0.0, ymax=40)
# 	xlim(left=0.0)
# 	zlim(bottom=0.0)
# 	legend(bbox_to_anchor=(1.05, 1), borderaxespad=5)
# 	# ax.view_init(elev=20.0, azim=45)
# 	tight_layout()
# 	savefig("figz/material_space.pdf", bbox_inches="tight")
# 	gcf()
	
	x = zeros(3, length(xtal_names))
	props = [L"$M$ [mol/kg]", L"$K$ [bar$^{-1}$]", L"$\rho_{ads}$ [kg/m³]"]
	for (c, xtal_name) in enumerate(xtal_names)
		x[1, c] = xtal_to_M[xtal_name]
		x[2, c] = xtal_to_K[xtal_name]
		x[3, c] = xtal_to_ρ[xtal_name]
	end
	
	fig, axs = subplots(3, 3, figsize=(8, 8))
	for i = 1:3
		for j = 1:3
			if i < j
				axs[i, j].axis("off")
				continue
			end
			if i == j
				axs[i, j].hist(x[i, :], color="k")
				axs[i, j].set_ylabel("# adsorbents")
				axs[i, j].set_xlabel(props[i])
				axs[i, j].set_xlim(xmin=0)
			else
				for c = 1:length(xtal_names)
					axs[i, j].scatter(x[i, c], x[j, c];
						scatter_kwargs(xtal_names[c])...)
				end
				axs[i, j].set_xlabel(props[i])
				axs[i, j].set_ylabel(props[j])
				axs[i, j].set_xlim(xmin=0)
				axs[i, j].set_ylim(ymin=0)
				if i == 3
					handles, labels = axs[i, j].get_legend_handles_labels()
					fig.legend(handles, labels)
				end
			end
		end
	end
	tight_layout()
	savefig("figz/material_space.pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ 53ea8fba-189f-11eb-381f-293a22ca726a
materials_space_viz()

# ╔═╡ 47eef0be-1973-11eb-399a-a5657c838de7
function silver_lining()
	
	Ps = range(0.0, 50.0, length=100)
	fig, axs = subplots(2, 1, sharex=true)
	axs[2].set_xlabel(L"storage pressure, $P$ [bar]")
	
	# tf 
	axs[1].set_ylabel("tankage fraction")
	for xtal_name in xtal_names
		tf_ads = tankage_fract.(Ps, xtal_name)
		axs[1].plot(Ps, tf_ads, color=xtal_to_color[xtal_name])
	end
	axs[1].plot(Ps, tankage_fract.(Ps), linestyle="--", color="gray")
	axs[1].set_yscale("log")
	
	# r
	axs[2].set_ylabel(L"radius, $r$, [m]")
	for xtal_name in xtal_names
		r_ads = r.(Ps, xtal_name)
		axs[2].plot(Ps, r_ads, color=xtal_to_color[xtal_name])
	end
	axs[2].plot(Ps, r.(Ps), linestyle="--", color="gray")
	
	axs[1].set_xlim(minimum(Ps), maximum(Ps))
	gcf()
end

# ╔═╡ aa7d343e-1973-11eb-094e-735edc33dd81
silver_lining()

# ╔═╡ 992aaf76-2f77-11eb-0d5f-3f19e5205ab6
md"### Compare crystal and bulk densities

compare xtal density with reported tap/bulk density for materials for which we have it.
source: https://pubs.acs.org/doi/full/10.1021/acsami.0c11200 pg. S-13 of SI has table.
"

# ╔═╡ c8a8abb4-35ac-11eb-3552-2fa75cb611f9
xtal_to_bulk_ρ = Dict("CC3"       => 0.139 * 1000,
					  "HKUST-1"   => 0.35 * 1000,
					  "Ni-MOF-74" => 0.539 * 1000)

# ╔═╡ fd5c8c0e-35ac-11eb-040c-1301adb58152
begin
	function compare_bulk_xtal_ρ()
		w = 0.35 # width
		bulk_ρ_xtal_names = keys(xtal_to_bulk_ρ)
		ids = 1:length(bulk_ρ_xtal_names)
		
		ρ_xtal = [xtal_to_ρ[xtal_name] for xtal_name in bulk_ρ_xtal_names]
		ρ_bulk = [xtal_to_bulk_ρ[xtal_name] for xtal_name in bulk_ρ_xtal_names]
		colors = [xtal_to_color[xtal_name] for xtal_name in bulk_ρ_xtal_names]
		
		figure(figsize=(4.4, 4.8))
		xlabel("adsorbents")
		ylabel(L"adsorbent density, $\rho_{ads}$ [kg/m$^3$]")
		
		bar(ids, ρ_xtal, 
			w, color=colors, edgecolor=colors, hatch="")
		bar(ids .+ w, ρ_bulk, 
			w, fill=false, edgecolor=colors, hatch="////", linewidth=1.0)
		xticks(ids .+ w/2, 
			[xtal_to_label[xtal_name] for xtal_name in bulk_ρ_xtal_names],  
			rotation="vertical", fontsize=10)
		
		
		# dummy plot for legend
		xlims = gca().get_xlim() # store x-lims to put back later
		bar([-2], [0], w, color="k", edgecolor="k", hatch="", label="crystal")
		bar([-2], [0], w, fill=false, edgecolor="k", hatch="////", linewidth=1.0, label="bulk")
		legend()
		xlim(xlims)

		tight_layout()
		savefig("figz/xtal_vs_bulk_density.pdf", format="pdf")
		gcf()
	end
	
	compare_bulk_xtal_ρ()
end

# ╔═╡ 32913bc4-35d3-11eb-0923-950b040de4f0
md"### Redo previous analysis using new densities
"

# ╔═╡ ed2754f6-3f46-11eb-0039-2b0419e7f318
begin
	for xtal_name in keys(xtal_to_bulk_ρ)
		xtal_name_bulk = xtal_name * "_bulk"
		
		xtal_to_M[xtal_name_bulk] = xtal_to_M[xtal_name]
		xtal_to_K[xtal_name_bulk] = xtal_to_K[xtal_name]
		xtal_to_ρ[xtal_name_bulk] = xtal_to_bulk_ρ[xtal_name]
		
		opt = optimum_storage(xtal_name_bulk)
			
		ads_opt[xtal_name_bulk] = Dict()
		ads_opt[xtal_name_bulk]["P [bar]"] = opt.minimizer
		ads_opt[xtal_name_bulk]["tf"] = opt.minimum / mass_desired_xe_propellant
		
		@assert ads_opt[xtal_name_bulk]["tf"] ≈ tankage_fraction_opt_analytical(xtal_name_bulk)
		@assert isapprox(ads_opt[xtal_name_bulk]["P [bar]"], P_opt_analytical(xtal_name_bulk), atol=0.001)
	end
end

# ╔═╡ 1fdb33fe-3f56-11eb-1864-c142c3315bfd
function bulk_vs_xtal_ρ_on_tf()
	figure(figsize=figsize)
	xlabel(L"density of material, $\rho_{ads}$ [kg/m$^3$]")
	ylabel("optimal tankage fraction")
	texts = []
	for xtal_name in keys(xtal_to_bulk_ρ)
		xtal_name_bulk = xtal_name * "_bulk"
		
		ρ_xtal = xtal_to_ρ[xtal_name]
		ρ_bulk = xtal_to_ρ[xtal_name_bulk]
		
		tf_opt_xtal = ads_opt[xtal_name]["tf"]
		tf_opt_bulk = ads_opt[xtal_name_bulk]["tf"]
		scatter(ρ_xtal, tf_opt_xtal; scatter_kwargs(xtal_name)..., zorder=10)
		scatter(ρ_bulk, tf_opt_bulk; scatter_kwargs(xtal_name)..., color="k", edgecolor="k")
		plot([ρ_xtal, ρ_bulk], [tf_opt_xtal, tf_opt_bulk], linestyle="--", color="0.1")
		push!(texts,
			  annotate(xtal_to_label[xtal_name], (ρ_xtal, tf_opt_xtal), fontsize=10)
		)
	end
	
	axhline(y=[bulk_opt["tf"]], linestyle="--", color="gray", label="bulk storage")
	adjustText.adjust_text(texts, force_points=(4.0, 2.0))
	
	ylim(ymin=0.0)
	xlim(xmin=0.0)
	tight_layout()
	# legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
	savefig("figz/tf_vs_crystal_and_bulk_density.pdf", format="pdf", bbox_inches="tight")
	gcf()
end

# ╔═╡ 5a9aff8e-3f57-11eb-127a-e934dbed5245
bulk_vs_xtal_ρ_on_tf()

# ╔═╡ Cell order:
# ╟─54e4a2c4-0b62-11eb-3e1d-017a70ae43d7
# ╠═6091fca2-0b62-11eb-0351-93554cbc8d52
# ╠═1cffb8d4-18a3-11eb-087f-83393730c57b
# ╠═68ff0da6-0b66-11eb-1911-bf3faa3ceb7a
# ╠═422c02b2-1802-11eb-31a7-c3671345f577
# ╠═0c7cac6e-1802-11eb-122a-0180ae055dbc
# ╟─6f90c580-0b62-11eb-092a-3f99accec091
# ╠═6b759e8a-0b62-11eb-059d-8325bf28a129
# ╟─e4f7db60-0b65-11eb-1dda-37f708ccc230
# ╠═ed862216-0b65-11eb-188a-ab1107ab8e24
# ╠═fdd5dc94-0b65-11eb-1907-47bbff769958
# ╟─a4a6f3be-0f65-11eb-0c80-53046792dd1f
# ╠═2c25ab68-0b66-11eb-0f41-a91ffdf1e371
# ╠═5bbe5418-0b66-11eb-1d57-798a7b5b78de
# ╟─29239f5e-0b63-11eb-1b90-67892ca5c215
# ╠═7cc89760-0b64-11eb-2360-c97a464d908f
# ╟─b5cef332-0f59-11eb-22c4-cf85220a7af5
# ╠═be69b9dc-0f59-11eb-2361-69b0f7e1ccf3
# ╟─b853e33e-0f35-11eb-2027-9d77b137546e
# ╠═eb9159f2-0f35-11eb-156e-d3762896ec0d
# ╠═0be36a3a-0f36-11eb-0d80-d11b529da9ad
# ╟─eed22218-0f5b-11eb-15a1-a75120aaf145
# ╠═40efd702-0f5b-11eb-2970-c1ee6ac4eedf
# ╠═0fdab85c-0f5c-11eb-2cec-69e1a95c6e45
# ╟─fb3bb86c-0f5c-11eb-0438-57bba5356391
# ╠═06ad0a18-0f5d-11eb-32bb-874fc19bb525
# ╟─d336fcca-0f5f-11eb-3cb2-7f5c3f040f7b
# ╠═331c7dda-0dcc-11eb-1658-d7c4d23b788e
# ╠═feaf132e-0b63-11eb-28ff-7f327d815001
# ╟─fa9fe136-0b64-11eb-3718-3527c3599f56
# ╠═ef807360-0b64-11eb-0f1b-cdd324b8b1bb
# ╠═e7a49568-1893-11eb-321c-69b330620a80
# ╠═f5711acc-1893-11eb-1937-337d2b173321
# ╟─2d85d0dc-0b66-11eb-1a99-7fee2143ef57
# ╠═83c86fb4-3a7c-11eb-281d-0b0fac556a7c
# ╟─aa5c1ad4-0b67-11eb-37b5-9746aa76c7f8
# ╠═218243b0-0dc8-11eb-01ed-258156f8aa9d
# ╠═9c3c5810-0b67-11eb-0141-8bfa134546a6
# ╠═2f55d9e6-0b68-11eb-271a-b502d0bc5a33
# ╠═33ed9f0e-0dc8-11eb-1b06-89178c499937
# ╟─d90fc4f0-0b66-11eb-33da-372fc8695a3c
# ╠═e1dbf234-0b66-11eb-244a-33506ae7db4a
# ╟─11185356-0b67-11eb-1bb3-e15f94234be3
# ╠═1721c25a-0b67-11eb-1b47-9bad136bd63a
# ╟─f4cb6ea8-0da6-11eb-120e-8d67d8ca58b7
# ╠═7fe2a354-0f33-11eb-2f83-e3b66bbdb333
# ╠═018c2b64-0da7-11eb-100a-4d1ee4e512f3
# ╟─d182f874-18bb-11eb-28b8-ad20c019015e
# ╠═e4ed9b6c-18bb-11eb-2796-8fe1689d5e35
# ╠═c5a0c9e2-18a1-11eb-17a4-8572d491905e
# ╠═1015e5fa-0dcd-11eb-3352-09bc6056a02f
# ╠═95054c18-2b93-11eb-0d17-3f8166d00547
# ╠═af8c02fc-2b93-11eb-3e3e-8f4695ff5bbf
# ╠═ce3da226-18a3-11eb-25df-572b24f5d3f0
# ╠═a8b0c9c4-0dcd-11eb-2924-91b2875da4ea
# ╟─fb16b62c-0f32-11eb-05fd-e72d057be910
# ╠═37cf428e-0dce-11eb-154b-25742d645575
# ╠═105f4f12-0f33-11eb-161c-a1d26651a48c
# ╠═11b18cca-0f33-11eb-2556-f946a9ac7d0b
# ╠═e9c0dacc-189e-11eb-0a0b-49ea73af6ef7
# ╠═53ea8fba-189f-11eb-381f-293a22ca726a
# ╠═47eef0be-1973-11eb-399a-a5657c838de7
# ╠═aa7d343e-1973-11eb-094e-735edc33dd81
# ╟─992aaf76-2f77-11eb-0d5f-3f19e5205ab6
# ╠═c8a8abb4-35ac-11eb-3552-2fa75cb611f9
# ╠═fd5c8c0e-35ac-11eb-040c-1301adb58152
# ╟─32913bc4-35d3-11eb-0923-950b040de4f0
# ╠═ed2754f6-3f46-11eb-0039-2b0419e7f318
# ╠═1fdb33fe-3f56-11eb-1864-c142c3315bfd
# ╠═5a9aff8e-3f57-11eb-127a-e934dbed5245
