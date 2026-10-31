core = normpath(joinpath(@__DIR__, "../../DataDrivenModel/core"))
include(joinpath(core, "0_header.jl"))
include(joinpath(core, "1_utils.jl"))
isdefined(@__MODULE__, :STLSQresult) || include(joinpath(core, "3_SINDy.jl"))
isdefined(@__MODULE__, :Reservoir) || include(joinpath(core, "3_RC.jl"))
include(joinpath(@__DIR__, "06_analyze.jl"))
@load joinpath(@__DIR__, "../data.jld2") iso3 ISO3 gdppc landarea

# plt_basic = plot(legend = :none, lims = [2, Inf], size = [400, 400])
# for (k, v) in iso3
#     plot!(plt_basic, v.x, v.y, v.z, xlabel = "x", ylabel = "y", zlabel = "z", )
#     scatter!(plt_basic, v.x[[end]], v.y[[end]], v.z[[end]], label = "", msc = :white)
# end
# plt_basic

vrbl = [:x, :y, :z]
train_years, test_years = 1950:2019, 2020:2024
train, test = model_periods(iso3, ISO3; train_years, test_years, vars = vrbl)
output = joinpath(@__DIR__, "../archive/models_$(first(train_years))_$(last(train_years))")

rbf_specs = [(ID = "R_g$(n)_s$(Float64(s))", n_basis = n, width_scale = Float64(s))
    for n in 25:25:100 for s in 1:5]
rc_specs = [(ID = "RC_n$(n)_a$(α)", n_nodes = n, alpha = α)
    for n in (100, 150, 200, 250) for α in (0.5, 0.7, 1.0)]
mlp_specs = [(ID = "M_sig_h$(h)_w$(w)", hidden_layers = h, width = w)
    for h in 1:3 for w in (32, 64, 128, 256, 512, 1024)]

result_rbf = run_model_grid(:RBF, rbf_specs, train, test, output; vars = vrbl)
result_rc = run_model_grid(:RC, rc_specs, train, test, output; vars = vrbl)
result_mlp = run_model_grid(:MLP, mlp_specs, train, test, output; vars = vrbl)

f1, f2, f3 = best_model.((:RBF, :RC, :MLP), (result_rbf, result_rc, result_mlp), Ref(output))
CSV.write("result_rc.csv", result_rc); sort(result_rc, :test_rmse)
CSV.write("result_mlp.csv", result_mlp); sort(result_mlp, :test_rmse)
CSV.write("result_rbf.csv", result_rbf); sort(result_rbf, :test_rmse)
f4 = ensemble(f1, f2, f3) # Equal-weight mean of the three annual maps.
f_ = Dict("RBF" => f1, "RC" => f2, "MLP" => f3, "ENS" => f4)

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                fixed point analysis

''''''''''''''''''''''''''''''''''''''''''''''''''"""

result_fps = DataFrame(model = [], x = [], y = [], z = [])
for (m, f) in f_
    fps = findfp(f)
    @info "Found $(length(fps)) fixed points for model $m"
    for fp in fps
        push!(result_fps, (m, fp[1], fp[2], fp[3]))
    end
end
m_ = ["RBF", "RC", "MLP", "ENS"]
plt_fp = plot(lims = [5, 9], camera = [15, 30])
for k in 1:4
    rf = result_fps[result_fps.model .== m_[k], :]
    scatter!(plt_fp, rf.x, rf.y, rf.z, color = [:white, :black, :black, :black][k], label = m_[k], shape = [:star6, :+, :x, :rect][k])
end
plt_fp


# Basin numbers follow increasing z, matching country_basins and the world map.
efps = sort(Float64.(result_fps[result_fps.model .== "ENS", vrbl]), :z)
basins2024 = country_basins(f_["ENS"], iso3, ISO3, eachrow(efps);
    year = 2024, warmup_start = first(train_years))
basin_counts2024 = sort(combine(groupby(basins2024, [:point, :point_name, :status]), nrow => :countries), :point)
CSV.write("country_basins_2024.csv", basins2024[:, [1,3,4,5,6]], bom = true)
basins2024 = CSV.read("country_basins_2024.csv", DataFrame, stringtype = String)

x_ = [basins2024.x2024[basins2024.point .== p] for p in 1:nrow(efps)]
y_ = [basins2024.y2024[basins2024.point .== p] for p in 1:nrow(efps)]
z_ = [basins2024.z2024[basins2024.point .== p] for p in 1:nrow(efps)]
i_ = [String.(basins2024.ISO3[basins2024.point .== p]) for p in 1:nrow(efps)]
color_ = [0.12 0.56 1.0; # dodgerblue
          1.0 0.39 0.28; # tomato
          0.2 0.8 0.2]  # limegreen

using MATLAB
begin
mat"""
figure(); ax_basin_tips = axes(); hold(ax_basin_tips, 'on'); drawnow
"""
for i in eachindex(x_)
    c = color_[i:i, :]
    mat"""
    s = scatter3(ax_basin_tips, $(x_[i]), $(y_[i]), $(z_[i]), 20, $c, 'filled', 'MarkerEdgeColor', 'none');
    s.DataTipTemplate.DataTipRows = [
        dataTipTextRow('국가', $(i_[i]))
        dataTipTextRow('x', 'XData')
        dataTipTextRow('y', 'YData')
        dataTipTextRow('z', 'ZData')
    ];
    drawnow
    """
end
mat"""
scatter3(ax_basin_tips, $(efps.x), $(efps.y), $(efps.z), 100, 'k', 'square', 'filled');
xlabel(ax_basin_tips, 'x'); ylabel(ax_basin_tips, 'y'); zlabel(ax_basin_tips, 'z');
xlim(ax_basin_tips, [4.5 9.5]); ylim(ax_basin_tips, [4.5 9.5]); zlim(ax_basin_tips, [4.5 9.5]);
view(ax_basin_tips, 3); grid(ax_basin_tips, 'on'); box(ax_basin_tips, 'on'); drawnow
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


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                auxiliary data

''''''''''''''''''''''''''''''''''''''''''''''''''"""

_basins2024 = leftjoin(basins2024, gdppc, on = :ISO3)
_basins2024 = leftjoin(_basins2024, landarea, on = :ISO3)
dropmissing!(_basins2024)

ecnm = combine(groupby(_basins2024, :point), :g2024 => mean => :gmean)
geom = combine(groupby(_basins2024, :point), :a2024 => mean => :amean)

bar(ecnm.gmean, color = [:dodgerblue, :tomato, :limegreen], xlabel = "type", ylabel = "average GDP per capita (current US\$)", legend = :none, lw = 0, xticks = 1:3)
bar(geom.amean, color = [:dodgerblue, :tomato, :limegreen], xlabel = "type", ylabel = "average land area (sq. km)", legend = :none, lw = 0, xticks = 1:3)


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    plane analysis

''''''''''''''''''''''''''''''''''''''''''''''''''"""

