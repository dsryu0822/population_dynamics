include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join = true)[[1,2]])
isdefined(@__MODULE__, :STLSQresult) || include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_SINDy.jl")
include(joinpath(@__DIR__, "06_analyze.jl"))
isdefined(@__MODULE__, :Reservoir) || include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_RC.jl")
include(joinpath(@__DIR__, "model_tools.jl"))
@load joinpath(@__DIR__, "../data.jld2") iso3

vrbl = [:x, :y, :z]
train_years, test_years = 1950:2019, 2020:2024
country_ids = readlines(joinpath(@__DIR__, "14_country_ids.txt"))
train14, test14 = model_periods(iso3, country_ids; train_years, test_years, vars = vrbl)
states14 = vcat([train14[k] for k in sort(country_ids)]...)
output14 = joinpath(@__DIR__, "../archive/models_2000_2019")
epochs14, seed14, rate14 = 5000, 0, 1e-3

rc_specs = [(ID = "RC_n$(n)_a$(α)", n_nodes = n, alpha = α)
    for n in (100, 150, 200, 250) for α in (0.5, 0.7, 1.0)]
mlp_specs = [(ID = "M_sig_h$(h)_w$(w)", hidden_layers = h, width = w, activation = "sigmoid",
    epochs = epochs14, seed = seed14, rate = rate14) for h in 1:3 for w in (32, 64, 128, 256)]
rbf_specs = [(ID = "R_g$(n)_s$(Float64(s))", n_basis = n, width_scale = Float64(s),
    constant = false, linear = false) for n in 25:25:100 for s in 1:5]

# Train: 19 one-step transitions. Test: five recursive steps, with no test inputs.
# RC resets per country and uses all training observations to initialize its memory.
result_rc = run_model_grid14(:RC, rc_specs, train14, test14, states14, output14;
    filename = "result_rc.csv", vars = vrbl)
result_mlp = run_model_grid14(:MLP, mlp_specs, train14, test14, states14, output14;
    filename = "result_mlp.csv", vars = vrbl)
result_rbf = run_model_grid14(:RBF, rbf_specs, train14, test14, states14, output14;
    filename = "result_rbf.csv", vars = vrbl)

JLD2.jldsave(joinpath(output14, "configuration.jld2"); train14, test14, country_ids, train_years, test_years,
    rc_specs, mlp_specs, rbf_specs, fp_rule = "nearest [5,5,5] or [7,7,7]; ties to fp5", metric = "train one-step; test recursive log10 RMSE")

f2, f3, f4 = best_model14.((:RC, :MLP, :RBF), (result_rc, result_mlp, result_rbf), Ref(output14))
CSV.write("result_rc.csv", result_rc); sort(result_rc, :test_rmse)
CSV.write("result_mlp.csv", result_mlp); sort(result_mlp, :test_rmse)
CSV.write("result_rbf.csv", result_rbf); sort(result_rbf, :test_rmse)

fixedpoints = vcat(map(("RC", "MLP", "RBF"), (result_rc, result_mlp, result_rbf)) do model, result
    valid = filter(:status => ==("ok"), result)
    id = valid.ID[argmin(valid.test_rmse)]
    points = select(JLD2.load(joinpath(output14, "weights", id * ".jld2"), "points"),
        :fp_group, :x, :y, :z, :stability, :spectral_radius, :growth_rate, :residual, :in_data_bounds)
    insertcols(points, 1, :ID => id, :model => model)
end...)
fps = filter(:stability => ==("stable"), fixedpoints)

CSV.write("fps.csv", fps)

tempcolor = ([:white, :black])[1 .+ (fixedpoints.stability .== "stable")]
scatter(fixedpoints.x, fixedpoints.y, fixedpoints.z, group = fixedpoints.model, msc = tempcolor)
scatter(fps.x, fps.y, fps.z, group = fps.model, xlabel = "x", lims = [5, Inf])

f2([4.0, 4.0, 4.0])
f3([4.0, 4.0, 4.0])
f4([4.0, 4.0, 4.0])

f5 = ensemble(f2, f3, f4) # Equal-weight mean of the three annual maps.
fps5 = findfp(f5, bounds = ((5, 9), (5, 9), (5, 9)), ngrid = 5) # Stable observed states; RC stability includes its reservoir.

efps = DataFrame(stack(fps5, dims = 1), [:x, :y, :z])
scatter!(efps.x, efps.y, efps.z, label = "ensemble")

plt_basic = plot(legend = :none, lims = [2, Inf], size = [400, 400])
for (k, v) in iso3
    plot!(plt_basic, v.x, v.y, v.z, xlabel = "x", ylabel = "y", zlabel = "z", )
    scatter!(plt_basic, v.x[[end]], v.y[[end]], v.z[[end]], label = "", msc = :white)
end
plt_basic

basins2024 = country_basins(f5, iso3, country_ids, fps5; year = 2024, warmup_start = first(train_years))
basin_counts2024 = combine(groupby(basins2024, [:point, :point_name, :status]), nrow => :countries)
CSV.write("country_basins_2024.csv", basins2024, bom = true)
sort(basin_counts2024, :point)

sort(basins2024, :point)


extrema(Matrix(basins2024[:, 3:5]))

scatter(basins2024.x2024, basins2024.y2024, basins2024.z2024, group = basins2024.point, ms = 3, msw = 0)
scatter!(efps.x, efps.y, efps.z, color = 1:3, label = :none, shape = :square)

x_ = [basins2024.x2024[basins2024.point .== p] for p in 1:3]
y_ = [basins2024.y2024[basins2024.point .== p] for p in 1:3]
z_ = [basins2024.z2024[basins2024.point .== p] for p in 1:3]
i_ = [basins2024.ISO3[basins2024.point .== p] for p in 1:3]

using MATLAB
mat"""
figure()

scatter3($(x_[1]), $(y_[1]), $(z_[1]), 20, 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter3($(x_[2]), $(y_[2]), $(z_[2]), 20, 'filled', 'MarkerEdgeColor', 'k')
scatter3($(x_[3]), $(y_[3]), $(z_[3]), 20, 'filled', 'MarkerEdgeColor', 'k')

scatter3($(efps.x[1]), $(efps.y[1]), $(efps.z[1]), 100, 'k', 'square', 'filled')
scatter3($(efps.x[2]), $(efps.y[2]), $(efps.z[2]), 100, 'k', 'square', 'filled')
scatter3($(efps.x[3]), $(efps.y[3]), $(efps.z[3]), 100, 'k', 'square', 'filled')

xlabel('x'); ylabel('y'); zlabel('z'); xlim([4.5 9.5]); ylim([4.5 9.5]); zlim([4.5 9.5])
view(3); grid on; box on
"""



x_ = [basins2024.x2024[basins2024.point .== p] for p in 1:3]
y_ = [basins2024.y2024[basins2024.point .== p] for p in 1:3]
z_ = [basins2024.z2024[basins2024.point .== p] for p in 1:3]
i_ = [String.(basins2024.ISO3[basins2024.point .== p]) for p in 1:3]
color_ = [0.12 0.56 1.0; # dodgerblue
          1.0 0.39 0.28; # tomato
          0.2  0.8  0.2] # limegreen
begin
mat"""
figure(); hold on
"""
for i in 1:3
    c = color_[i:i, :]   # 1×3 행렬로 유지
    mat"""
    s = scatter3($(x_[i]), $(y_[i]), $(z_[i]), 20, $c, 'filled', 'MarkerEdgeColor', 'none');
    s.DataTipTemplate.DataTipRows = [
        dataTipTextRow('국가', $(i_[i]))
        dataTipTextRow('x', 'XData')
        dataTipTextRow('y', 'YData')
        dataTipTextRow('z', 'ZData')
    ];
    """
end
for i in 1:3
    mat"""
    scatter3($(efps.x[i]), $(efps.y[i]), $(efps.z[i]), 100, 'k', 'square', 'filled');
    """
end
mat"""
xlabel('x'); ylabel('y'); zlabel('z');
xlim([4.5 9.5]); ylim([4.5 9.5]); zlim([4.5 9.5])
view(3); grid on; box on
"""
end



"""''''''''''''''''''''''''''''''''''''''''''''''''''

                          worldmap

''''''''''''''''''''''''''''''''''''''''''''''''''"""

iso = basins2024.ISO3
pt  = basins2024.point

locs = "[" * join(("\"$c\"" for c in iso), ",") * "]"
zs   = "[" * join(pt, ",") * "]"

html = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
</head>
<body>
  <div id="map" style="width:100%;height:95vh"></div>
  <script>
    Plotly.newPlot("map", [{
      type: "choropleth",
      locationmode: "ISO-3",
      locations: $locs,
      z: $zs,
      zmin: 0.5,
      zmax: 3.5,
      colorscale: [
        [0, "dodgerblue"], [1/3, "dodgerblue"],
        [1/3, "tomato"],   [2/3, "tomato"],
        [2/3, "limegreen"], [1, "limegreen"]
      ],
      showscale: false
    }], {geo: {projection: {type: "equal earth"}}});
  </script>
</body>
</html>
"""
write("worldmap.html", html)

