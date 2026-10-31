include("01_SINDy.jl")
include("02_NGRC.jl")
include("03_MLP.jl")
include("04_RBF.jl")

@assert all(isdefined(@__MODULE__, name) for name in (:f1, :f2, :f3, :f4)) "Run the model blocks in 01–04 first."
k = "KOR"
comparison_data = population_data(iso3; data_config...)
snapshot_names = (:data1, :data2, :data3, :data4)
for name in snapshot_names
    if isdefined(@__MODULE__, name)
        @assert isequal(getfield(@__MODULE__, name), comparison_data) "$name differs from data_config. Check that model's training data."
    end
end
missing_snapshots = [name for name in snapshot_names if !isdefined(@__MODULE__, name)]
if !isempty(missing_snapshots)
    @info "Using the same evaluation data; past training data cannot be verified for these models." missing_snapshots
end
v = comparison_data[k]
vrbl = setdiff(names(v), ["t"])
U = Matrix(v[1:2, vrbl])
actual = Matrix(v[3:end, vrbl])
pred1 = forecast(f1, U, size(actual, 1))
pred2 = forecast(f2, U, size(actual, 1); discrete = true)
pred3_k = forecast(f3, U, size(actual, 1); discrete = true, lags = 1)

# Same evaluation window for all models; saved training data checked when available.
pred4_k = forecast(f4, U, size(actual, 1); discrete = true, lags = 1)
predictions = ["SINDy" => pred1, "NGRC" => pred2, "MLP" => pred3_k, "RBF" => pred4_k]
comparison = DataFrame(model = first.(predictions),
    rmse = [sqrt(mean(abs2, p - actual)) for (_, p) in predictions])
display(comparison) # log10-population units
plt_comparison = prediction_plot(v[3:end, :],
    predictions; vars = vrbl, title = k)
display(plt_comparison)

if all(isdefined(@__MODULE__, name) for name in (:fixedpoints1, :fixedpoints2, :fixedpoints3, :fixedpoints4))
    fixedpoints_comparison = vcat(
        insertcols(copy(fixedpoints1), 1, :model => "SINDy"),
        insertcols(copy(fixedpoints2), 1, :model => "NGRC"),
        insertcols(copy(fixedpoints3), 1, :model => "MLP"),
        insertcols(copy(fixedpoints4), 1, :model => "RBF"))
    display(select(fixedpoints_comparison,
        :model, :x, :y, :z, :residual, :in_data_bounds, :growth_rate, :stability))
end
