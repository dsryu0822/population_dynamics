include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2]])
if !isdefined(@__MODULE__, :STLSQresult)
    include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_SINDy.jl")
end
include(joinpath(@__DIR__, "06_analyze.jl"))
using Flux
@load joinpath(@__DIR__, "../data.jld2") ISO3 iso3

KFOLD, EPOCH = 10, 5000
fold_seed, model_seed = 0, 0 # hold initialization fixed to isolate fold sensitivity
λ_sindy, poly_degree = 1e-5, 5
match_tol = 0.15 # Euclidean distance in the original log10(x,y,z) coordinates
vars = data_config.vars
countries = population_data(iso3; data_config...)
data = transition_data(countries; vars)
data.fold = zeros(Int, nrow(data))
add_fold!(data; k = KFOLD, seed = fold_seed)
states = vcat([countries[k] for k in sort(collect(keys(countries)))]...)
vrbl = ("d" .* string.(vars), string.(vars))
cnfg = cook(vrbl; poly = 0:poly_degree)

MLP_, SINDy_, mlp_runs, roots_ = [], [], [], []
metrics = DataFrame(fold = Int[], ntrain = Int[], nvalid = Int[],
    mlp_train_next_rmse = Float64[], mlp_valid_next_rmse = Float64[],
    sindy_train_rate_rmse = Float64[], sindy_valid_rate_rmse = Float64[])
output_dir = joinpath(@__DIR__, "../archive/SINDy_MLP")
mkpath(output_dir)
CSV.write(joinpath(output_dir, "fold_assignments.csv"), data)
blas_threads = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    for fold in 1:KFOLD
        @info "Starting paired fold" fold KFOLD
        train_indices, valid_indices = findall(data.fold .!= fold), findall(data.fold .== fold)
        trng, valid = data[train_indices, :], data[valid_indices, :]
        X, Y = Matrix(trng[:,vars])', Matrix(trng[:,"next_" .* string.(vars)])'
        Xv, Yv = Matrix(valid[:,vars])', Matrix(valid[:,"next_" .* string.(vars)])'
        mlp = fit_mlp(X, Y; epochs = EPOCH, seed = model_seed)
        MLP = mlp_map(mlp)
        sindy = SINDy(trng, vrbl, cnfg; λ = λ_sindy)
        push!(MLP_, MLP); push!(SINDy_, sindy); push!(mlp_runs, mlp)
        push!(metrics, (fold, nrow(trng), nrow(valid),
            sqrt(mean(abs2, MLP(X)-Y)), sqrt(mean(abs2, MLP(Xv)-Yv)),
            sqrt(mean(abs2, residual(sindy,trng))), sqrt(mean(abs2, residual(sindy,valid)))))
        fp_mlp = find_fixedpoints(MLP, states; vars, discrete = true, lags = 1)
        fp_sindy = find_fixedpoints(sindy, states; vars)
        insertcols!(fp_mlp, 1, :fold => fold, :model => "MLP")
        insertcols!(fp_sindy, 1, :fold => fold, :model => "SINDy")
        push!(roots_, fp_mlp, fp_sindy)
        checkpoint = joinpath(output_dir, "fold_$(lpad(fold,2,'0')).jld2")
        JLD2.@save checkpoint mlp sindy train_indices valid_indices fp_mlp fp_sindy
        @info "Completed paired fold" fold mlp_roots = nrow(fp_mlp) sindy_roots = nrow(fp_sindy)
        CSV.write(joinpath(output_dir,"metrics.csv"), metrics)
    end
finally
    BLAS.set_num_threads(blas_threads)
end

fixedpoints_folds = vcat(roots_...)
insertcols!(fixedpoints_folds, 1, :id => 1:nrow(fixedpoints_folds))
stable_folds = fixedpoints_folds[(fixedpoints_folds.stability .== "stable") .& fixedpoints_folds.in_data_bounds, :]
matches = match_fixedpoints(stable_folds; tolerance = match_tol)
counts = DataFrame(fold = 1:KFOLD,
    mlp_stable = [count((fixedpoints_folds.fold .== k) .& (fixedpoints_folds.model .== "MLP") .& (fixedpoints_folds.stability .== "stable")) for k in 1:KFOLD],
    sindy_stable = [count((fixedpoints_folds.fold .== k) .& (fixedpoints_folds.model .== "SINDy") .& (fixedpoints_folds.stability .== "stable")) for k in 1:KFOLD],
    mlp_in_bounds = [count((stable_folds.fold .== k) .& (stable_folds.model .== "MLP")) for k in 1:KFOLD],
    sindy_in_bounds = [count((stable_folds.fold .== k) .& (stable_folds.model .== "SINDy")) for k in 1:KFOLD],
    matched = [count(matches.fold .== k) for k in 1:KFOLD])
plt_folds, plt_projections = fold_fixedpoint_plots(states, fixedpoints_folds, matches; vars, kfold = KFOLD)
display(counts)
display(plt_folds)
display(plt_projections)
savefig(plt_folds, joinpath(output_dir,"stable_by_fold.png"))
savefig(plt_projections, joinpath(output_dir,"stable_projections.png"))
CSV.write(joinpath(output_dir,"fixedpoints.csv"), select(fixedpoints_folds, Not(:eigenvalues)))
CSV.write(joinpath(output_dir,"stable_counts.csv"), counts)
CSV.write(joinpath(output_dir,"matches.csv"), matches)
JLD2.@save joinpath(output_dir,"results.jld2") data_config fold_seed model_seed KFOLD EPOCH λ_sindy poly_degree match_tol data metrics fixedpoints_folds matches counts
