include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2,5]])
include(joinpath(@__DIR__, "06_analyze.jl"))
@load joinpath(@__DIR__, "../data.jld2") ISO3 iso3

vrbl = data_config.vars
diso3 = population_data(iso3; data_config...)
data2 = deepcopy(diso3) # retain the exact observations used by f2
R, S = delayed(diso3; vars = vrbl)

@time f2_ = tikhonov_λ(R, S)
λ_ = exp10.(-10:0.5:3)
L_ = [mean(abs2, R * f2_(λ) - S) for λ in λ_]
λ = λ_[argmin(L_)]
f2 = NGRC(f2_(λ), λ)

pfmc = DataFrame(variable = vrbl,
    training_rmse = vec(sqrt.(mean(abs2, R * f2.Wout - S, dims = 1))))
display(pfmc) # fit on all selected observations
plt_loss = plot(λ_, L_, xscale = :log10, yscale = :log10,
    xlabel = "λ", ylabel = "training MSE", label = :none)

pred = Dict(k => forecast(f2, Matrix(v[1:2, vrbl]), nrow(v)-2; discrete = true)
    for (k, v) in diso3)
pfmc_rollout = DataFrame(ISO3 = sort(collect(keys(diso3))))
pfmc_rollout.rmse = [sqrt(mean(abs2, pred[k] - Matrix(diso3[k][3:end, vrbl])))
    for k in pfmc_rollout.ISO3]
display(pfmc_rollout)
k = haskey(diso3, "KOR") ? "KOR" : first(pfmc_rollout.ISO3)
plt_prediction = prediction_plot(diso3[k][3:end, :], ["NGRC" => pred[k]]; vars = vrbl, title = k)

data = vcat([diso3[k] for k in sort(collect(keys(diso3)))]...)
fixedpoints2 = find_fixedpoints(f2, data; vars = vrbl, discrete = true, residual_tol = 1e-10)
display(fixedpoints2)
