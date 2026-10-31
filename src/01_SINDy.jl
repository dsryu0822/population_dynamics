include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2,6]])
include(joinpath(@__DIR__, "06_analyze.jl"))
@load joinpath(@__DIR__, "../data.jld2") ISO3 iso3

vrbl = data_config.vars
diso3 = population_data(iso3; data_config...)
data1 = deepcopy(diso3) # retain the exact observations used by f1
data = vcat([add_diff(diso3[k][:, vrbl]) for k in sort(collect(keys(diso3)))]...)
cnfg = cook(half(names(data)), poly = 0:5)
f1 = SINDy(data, half(names(data)), cnfg; λ = 1e-5)
print(f1)

fixedpoints1 = find_fixedpoints(f1, data; vars = vrbl)
display(fixedpoints1)
