global PATH_TO_DATA = joinpath(pwd(), "data/NIST_data/experimental_isotherms/renamed_isotherm_files")

function isotherm_data(isotherm_filename::AbstractString)
    file = open(joinpath(PATH_TO_DATA, isotherm_filename))
    lines = readlines(file)
    close(file)
    data = JSON.parse(join(lines))
    
    pressure_pts = Array{Float64,1}()
    adsorption_pts = Array{Float64,1}()

    for datapoints_dict in data["isotherm_data"]
        append!(pressure_pts, datapoints_dict["pressure"])
        append!(adsorption_pts, datapoints_dict["total_adsorption"])
    end

    df_isotherm = Dict{Symbol, Array{Float64}}()
    df_isotherm[Symbol("P(bar)")] = pressure_pts
    df_isotherm[Symbol("Adsorption(mmol/g)")] = adsorption_pts
    DataFrame(df_isotherm)
end

dict_xe_isotherms[:COFORMATE] = isotherm_data("COFORMATE.json")

# df_xe_isotherms[Symbol("MOF-505")] = isotherm_data("MOF-505.json")

# I can't figure out why this isn't working

for mof in readdir(PATH_TO_DATA)
    destination = joinpath(PATH_TO_DATA,mof)
    name = Base.basename(destination, ".json")
    file = string(mof)
    df_xe_isotherms= isotherm_data(String(file))
end
