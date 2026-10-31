include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2]])
if !isdefined(@__MODULE__, :STLSQresult)
    include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_SINDy.jl")
end
include(joinpath(@__DIR__, "06_analyze.jl"))
using NearestNeighbors
@load joinpath(@__DIR__, "../data.jld2") ISO3 iso3

vrbl = data_config.vars
diso3 = population_data(iso3; data_config...)
data4 = deepcopy(diso3)
X, Y = map_data(diso3; vars = vrbl)

nleaves = 128 # exact number of terminal regions / RBF centers
width_scale = 1.0 # Gaussian sigma = leaf radius in the original x, y, z coordinates
linear = true # include constant and linear terms alongside the Gaussian candidates
λ4 = 1e-5 # STLSQ sparsity threshold, as in 01_SINDy.jl
# Build a fine tree, then cut it to the requested number of terminal regions.
tree4 = BallTree(X; leafsize = 1, reorder = false)
leaves4 = ball_leaves(tree4, nleaves)
counts4 = length.(ball_members.(leaves4))
balls4 = treeregion.(leaves4)
centers4 = stack(ball.center for ball in balls4)
radii4 = [ball.r for ball in balls4]
widths4 = width_scale .* max.(radii4, 1e-6) # zero-radius duplicate-point leaves
vrbl4 = ("next_" .* string.(vrbl), string.(vrbl))
data = hcat(DataFrame(X', last(vrbl4)), DataFrame(Y', first(vrbl4)))
cnfg4 = cook_rbf(vrbl4, centers4, widths4; poly = linear ? (0:1) : (0:0))
f4 = SINDy(data, vrbl4, cnfg4; λ = λ4)
print(f4)

pfmc4 = DataFrame(variable = vrbl,
    training_rmse = vec(sqrt.(mean(abs2, residual(f4, data), dims = 1))))
@show length(leaves4) extrema(counts4) extrema(radii4)
display(pfmc4)
pred4 = Dict(k => forecast(f4, Matrix(v[1:2, vrbl]), nrow(v)-2; discrete = true, lags = 1)
    for (k, v) in diso3)
k = haskey(diso3, "KOR") ? "KOR" : first(sort(collect(keys(diso3))))
plt_prediction4 = prediction_plot(diso3[k][3:end, :], ["RBF" => pred4[k]]; vars = vrbl, title = k)
display(plt_prediction4)
states4 = vcat([diso3[k] for k in sort(collect(keys(diso3)))]...)
fixedpoints4 = find_fixedpoints(f4, states4; vars = vrbl, discrete = true, lags = 1)
display(fixedpoints4)

output_dir = joinpath(@__DIR__, "../archive")
mkpath(output_dir)
rbf_discrete = true
JLD2.@save joinpath(output_dir, "RBF.jld2") f4 tree4 data4 fixedpoints4 pfmc4 nleaves counts4 width_scale λ4 rbf_discrete
savefig(plt_prediction4, joinpath(output_dir, "RBF_prediction.png"))
CSV.write(joinpath(output_dir, "RBF_fixedpoints.csv"), select(fixedpoints4, Not(:eigenvalues)))
