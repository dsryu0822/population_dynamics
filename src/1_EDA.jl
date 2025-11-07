include("0_datacall.jl")
using Plots
histogram(log10.(groupby(POP, :t)[end].x))
histogram(log10.(POP.x))
histogram(POP.x)

POP[POP.y .< 1, :]

# POP[POP.ISO3 .== "KOR", :]
plot(
    plot(data[data.ISO3 .== "KOR", :t], data[data.ISO3 .== "KOR", :x], ylabel = L"x", color = :black),
    plot(data[data.ISO3 .== "KOR", :t], data[data.ISO3 .== "KOR", :y], ylabel = L"y", color = :black),
    xlabel = L"t", legend = :none, layout = (2, 1), plot_title = "KOR: $(iso3["KOR"]), t ∈ [$(data[data.ISO3 .== "KOR", :t][1]), $(data[data.ISO3 .== "KOR", :t][end])] ",
    lw = 2, size = [600, 600]
)
png("ts")

DAT = data[data.ISO3 .== "KOR", :]
plot(plot(DAT.t, DAT.x, xlabel = "t", ylabel = "x", color = :black), 
     plot(DAT.t, DAT.y, xlabel = "t", ylabel = "y", color = :black), 
     plot(DAT.t, DAT.z, xlabel = "t", ylabel = "z", color = :black), 
     layout = (3, 1), size = (500, 500), legend = :none,
     plot_title = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]")
png("temp.png")
extrema(data.x)
scatter(combine(groupby(data, :ISO3), :z => last).z_last, assessment, lims = [2, 5.5], size = [400, 400], label = :none, color = :black)
sum(abs2.(combine(groupby(data, :ISO3), :z => last).z_last - assessment))
data[data.y .< 0, :]

count(isnan.(assessment))


function projection(x0, y0, a, b)
    x1 = (a/(1 + a^2))*(x0 + y0 - b)
    y1 = (1/(1 + a^2))*((a^2)*(x0 + y0) + b)
    return [x1, y1]
end

p2022 = stack([[df.x[end], df.y[end], df.z[end]] for df in groupby(data, :ISO3)])
prj2022 = stack(projection.(p2022[1, :], p2022[2, :], 1, 1/2))
s01 = plot(legend = :none, lims = [0, 6], xlabel = "x", ylabel = "y", size = [600, 600])
for k in axes(prj2022, 2)
    plot!(s01, [p2022[1, k], prj2022[1, k]], [p2022[2, k], prj2022[2, k]], color = :black, label = :none)
end
scatter!(p2022[1, :], p2022[2, :], ms = p2022[3, :], color = :white, size = [600, 600], legend = :none)
scatter!(s01, prj2022[1, :], prj2022[2, :], color = :black, ms = 1, lims = [0, 6])
s01

linedist = sqrt.(vec(sum(abs2, prj2022, dims = 1)))
GDPpc = p2022[3, :]
s02 = scatter(linedist, GDPpc, lims = [0, 7], size = [400, 400], smooth = true, legend = :none, color = :black)

out = lm(@formula(GDPpc ~ linedist), DataFrame(; linedist, GDPpc))
r2(out)
Rsq(fitted(out), GDPpc)

using SimpleWeightedGraphs
G = erdos_renyi(10, 0.1, is_directed=true)
WG = SimpleWeightedGraph(G)
WG.weights = WG.weights .* randn(10, 10)
adj = adjacency_matrix(WG)

eigenvector_centrality(WG)
