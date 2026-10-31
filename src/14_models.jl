core = normpath(joinpath(@__DIR__, "../../DataDrivenModel/core"))
include(joinpath(core, "0_header.jl"))
include(joinpath(core, "1_utils.jl"))
isdefined(@__MODULE__, :STLSQresult) || include(joinpath(core, "3_SINDy.jl"))
isdefined(@__MODULE__, :Reservoir) || include(joinpath(core, "3_RC.jl"))
include(joinpath(@__DIR__, "06_analyze.jl"))
@load joinpath(@__DIR__, "../data.jld2") iso3 ISO3 gdppc landarea
using MATLAB, Colors, LaTeXStrings
include(joinpath(@__DIR__, "01_plotlyutils.jl"))

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

# result_rbf = run_model_grid(:RBF, rbf_specs, train, test, output; vars = vrbl)
# result_rc = run_model_grid(:RC, rc_specs, train, test, output; vars = vrbl)
# result_mlp = run_model_grid(:MLP, mlp_specs, train, test, output; vars = vrbl)
# CSV.write("result_rc.csv", result_rc); sort(result_rc, :test_rmse)
# CSV.write("result_mlp.csv", result_mlp); sort(result_mlp, :test_rmse)
# CSV.write("result_rbf.csv", result_rbf); sort(result_rbf, :test_rmse)

result_rbf = CSV.read("result_rbf.csv", DataFrame)
result_rc = CSV.read("result_rc.csv", DataFrame)
result_mlp = CSV.read("result_mlp.csv", DataFrame)

f1, f2, f3 = best_model.((:RBF, :RC, :MLP), (result_rbf, result_rc, result_mlp), Ref(output))
f4 = ensemble(f1, f2, f3) # Equal-weight mean of the three annual maps.
f_ = Dict("RBF" => f1, "RC" => f2, "MLP" => f3, "ENS" => f4)
validate(f1); validate(f2); validate(f3); validate(f4);
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


efps = sort(Float64.(result_fps[result_fps.model .== "ENS", vrbl]), :z)
# basins2024 = country_basins(f_["ENS"], iso3, ISO3, eachrow(efps);
#     year = 2024, warmup_start = first(train_years))
# basin_counts2024 = sort(combine(groupby(basins2024, [:point, :point_name, :status]), nrow => :countries), :point)
# CSV.write("country_basins_2024.csv", basins2024[:, [1,3,4,5,6]], bom = true)
basins2024 = CSV.read("country_basins_2024.csv", DataFrame, stringtype = String)

x_ = [basins2024.x2024[basins2024.point .== p] for p in 1:nrow(efps)]
y_ = [basins2024.y2024[basins2024.point .== p] for p in 1:nrow(efps)]
z_ = [basins2024.z2024[basins2024.point .== p] for p in 1:nrow(efps)]
i_ = [String.(basins2024.ISO3[basins2024.point .== p]) for p in 1:nrow(efps)]
color_ = [0.12 0.56 1.0; # dodgerblue
          1.0 0.39 0.28; # tomato
          0.2 0.8 0.2]  # limegreen

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

ecnm = combine(groupby(_basins2024, :point), :g2024 => median => :gm)
geom = combine(groupby(_basins2024, :point), :a2024 => median => :am)

bar(ecnm.gm, color = [:dodgerblue, :tomato, :limegreen], xlabel = "type", ylabel = "median GDP per capita (current US\$)", legend = :none, lw = 0, xticks = 1:3)
bar(geom.am, color = [:dodgerblue, :tomato, :limegreen], xlabel = "type", ylabel = "median land area (sq. km)", legend = :none, lw = 0, xticks = 1:3)


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    plane analysis

''''''''''''''''''''''''''''''''''''''''''''''''''"""

plane = getplane(eachrow(efps)...)

# Overlay on the current MATLAB axes without changing the country markers or view.
mat"""
ax_plane = gca;
plane_limits = [xlim(ax_plane); ylim(ax_plane); zlim(ax_plane)];
plane_washeld = ishold(ax_plane);
plane_coeff = $plane;
[~, plane_axis] = max(abs(plane_coeff(1:3)));
plane_free = setdiff(1:3, plane_axis);
[plane_u, plane_v] = meshgrid(plane_limits(plane_free(1), :), plane_limits(plane_free(2), :));
plane_xyz = cell(1, 3);
plane_xyz{plane_free(1)} = plane_u;
plane_xyz{plane_free(2)} = plane_v;
plane_xyz{plane_axis} = -(plane_coeff(plane_free(1))*plane_u + ...
    plane_coeff(plane_free(2))*plane_v + plane_coeff(4)) / plane_coeff(plane_axis);
hold(ax_plane, 'on');
delete(findall(ax_plane, 'Type', 'surface', 'Tag', 'fixedpoint_plane'));
surf(ax_plane, plane_xyz{1}, plane_xyz{2}, plane_xyz{3}, ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off', 'HitTest', 'off', 'PickableParts', 'none', 'Tag', 'fixedpoint_plane');
xlim(ax_plane, plane_limits(1, :)); ylim(ax_plane, plane_limits(2, :)); zlim(ax_plane, plane_limits(3, :));
if ~plane_washeld, hold(ax_plane, 'off'); end
drawnow
"""

plane_data = deepcopy(_basins2024)
projections = [proj(row, plane) for row in eachrow(plane_data[:, [:x2024, :y2024, :z2024]])]
plane_data.l = - (Matrix(plane_data[:, [:x2024, :y2024, :z2024]]) * plane[1:3] .+ plane[4])
CSV.write("l2plane.csv", sort(plane_data, :l; rev=true)[:, [:ISO3, :l]]; bom = true)
plane_data.ud = (plane_data.l .> 0) .+ 1

scatter(plane_data.l, plane_data.g2024, yscale = :log10, size = [400, 400], legend = :none, color = ["#443FAD", "#EB3819"][plane_data.ud], msc = :white, xlabel = "distance to plane", ylabel = "GDP per capita (current US\$)")
bar(combine(groupby(plane_data, :ud), :g2024 => median => :gm).gm, color = ["#443FAD", "#EB3819"], ylabel = "median GDP per capita (current US\$)", legend = :none, lw = 0, xticks = (1:2, ["under", "over"]), size = [400, 400])


n      = length(plane_data.l)
idx    = [1:10; n-9:n]
ko = Dict("QAT" => "카타르", "ARE" => "아랍에미리트", "KWT" => "쿠웨이트",
          "BHR" => "바레인", "SGP" => "싱가포르", "OMN" => "오만",
          "MUS" => "모리셔스", "KOR" => "대한민국", "HKG" => "홍콩",
          "JAM" => "자메이카", "TCD" => "차드", "NGA" => "나이지리아",
          "TZA" => "탄자니아", "MLI" => "말리", "MOZ" => "모잠비크",
          "CAF" => "중앙아프리카공화국", "AGO" => "앙골라", "SOM" => "소말리아",
          "NER" => "니제르", "COD" => "콩고민주공화국")
labels = [ko[c] for c in plane_data.ISO3[idx]]
vals   = Float64.(plane_data.l[idx])
colors = [fill("#443FAD", 10); fill("#EB3819", 10)]
m      = maximum(abs, vals)
plotly_html(
    """[{
      type: "bar", orientation: "h",
      x: $(js(vals)),
      y: $(js(labels)),
      marker: {color: $(js(colors)), line: {width: 0}}
    }]""",
    """{
      width: 600, height: 600,
      showlegend: false,
      xaxis: {title: {text: "<i>l</i>"}, range: [$(-1.1m), $(1.1m)]},
      yaxis: {title: {text: "country"}, autorange: "reversed", type: "category"}
    }""",
)


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    military power

''''''''''''''''''''''''''''''''''''''''''''''''''"""