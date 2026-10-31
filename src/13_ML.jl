include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join = true)[[1,2,6]])
include(joinpath(@__DIR__, "06_analyze.jl")); include(joinpath(@__DIR__, "12_aNGRC.jl"))
using Flux, NearestNeighbors
@load joinpath(@__DIR__, "../data.jld2") iso3
vrbl = data_config.vars; cnfg2 = bake(vrbl; poly = 0:3, delay = 0:1)
EPOCH, nleaves, width_scale, λ2, λ4 = 5000, 100, 2.8, 1.95e-3, 1e-5
countries = population_data(iso3; data_config...)
states = vcat([insertcols(copy(countries[k]), :ISO3 => k) for k in sort(collect(keys(countries)))]...)
X, Y = map_data(Dict(k => v[maximum(metadata(cnfg2, "delay"))+1:end, :] for (k,v) in countries); vars = vrbl)
f2 = aNGRC(states, vrbl, cnfg2; tag = :ISO3, λ = λ2)
f3 = mlp_map(fit_mlp(X, Y; epochs = EPOCH, seed = 0))
balls = treeregion.(ball_leaves(BallTree(X; leafsize = 1, reorder = false), nleaves))
cnfg4 = cook_rbf(("next_" .* string.(vrbl), string.(vrbl)), stack(b.center for b in balls), width_scale .* max.([b.r for b in balls], 1e-6))
f4 = SINDy(hcat(DataFrame(X', string.(vrbl)), DataFrame(Y', "next_" .* string.(vrbl))), ("next_" .* string.(vrbl), string.(vrbl)), cnfg4; λ = λ4)
models = (; f2, f3, f4)
comparison = vcat([insertcols(filter(:stability => ==("stable"), find_fixedpoints(f, states; vars = vrbl, discrete = true, lags, ngrid = name == "RBF" ? 13 : 7)), 1, :model => name) for (name,f,lags) in (("aNGRC",f2,maximum(f2.delay)+1),("MLP",f3,1),("RBF",f4,1))]...)
display(comparison)

scatter(comparison.x, comparison.y, comparison.z, group = comparison.model)
