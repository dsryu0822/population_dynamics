include("../../DataDrivenModel/core/header.jl")
include("1_datacall.jl")

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                    This is 3D analysis

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""

dropmissing!(data)
add_fold!(data, seed = 0)
TD_ = [data[data.fold .!= k, :] for k in sort(unique(data.fold))]
VD_ = [data[data.fold .== k, :] for k in sort(unique(data.fold))]

vrbl = [:dx, :dy, :dz], [:x, :y, :z]

Tmse = []
Vmse = []
rankAnalysis = []
d_ = 0:5
for d = d_
    cnfg = cook(last(vrbl), poly = 0:d)
    ΘX = Θ(data[:, last(vrbl)], cnfg)
    push!(rankAnalysis, [rank(ΘX), size(ΘX, 2)])

    f_ = [SINDy(TD, vrbl, cnfg, λ = 0) for TD in TD_]
    push!(Tmse, mean([f.mse for f in f_]))
    push!(Vmse, mean([sum(abs2, stack(residual(f, VD))) / prod(size(VD[:, first(vrbl)])) for (f, VD) in zip(f_, VD_)]))
end
rankAnalysis
plot(xticks = d_, legend = :topleft)
plot!(Tmse, label = "training MSE")
plot!(Vmse, label = "validation MSE")
argmin(Vmse)

λ_ = exp10.(-4:0.1:-1)
λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:3); λ).mse for λ in λ_]
plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:.5:-1))

f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:3), λ = exp10(-2.5))
f |> print

for DAT in groupby(data, :ISO3)[1:1:end]
    tend = (DAT.t[end] - DAT.t[1])
    prdt = DataFrame(solve(f, [DAT[1, last(vrbl)]...], 0:0.01:tend), last(vrbl))
    @info "ISO3: $(DAT.ISO3[1])"
    
    pxy=plot(xlabel = "x", ylabel = "y", lims = [0, 0.5], legend = :topright)
        plot!(DAT.x, DAT.y, color = 1, label = :none); scatter!(DAT.x[[end]], DAT.y[[end]], color = 1, label = "data")
        plot!(prdt.x, prdt.y, color = 2, label = :none); scatter!(prdt.x[[end]], prdt.y[[end]], color = 2, label = "model")
    # pyz=plot(xlabel = "y", ylabel = "z", xlims = [0, 0.5], ylims = [1, 5.5], legend = :none)
    #     plot!(DAT.y, DAT.z, color = 1); scatter!(DAT.y[[end]], DAT.z[[end]], color = 1)
    #     plot!(prdt.y, prdt.z, color = 2); scatter!(prdt.y[[end]], prdt.z[[end]], color = 2)
    # pxz=plot(xlabel = "x", ylabel = "z", xlims = [0, 0.5], ylims = [1, 5.5], legend = :none)
    #     plot!(DAT.x, DAT.z, color = 1); scatter!(DAT.x[[end]], DAT.z[[end]], color = 1)
    #     plot!(prdt.x, prdt.z, color = 2); scatter!(prdt.x[[end]], prdt.z[[end]], color = 2)
    # xyz=plot(xlabel = "x", ylabel = "y", zlabel = "z", xlims = [0, 0.5], ylims = [0, 0.5], zlims = [1, 5.5], legend = :none)
    #     plot!(DAT.x, DAT.y, DAT.z, color = 1); scatter!(DAT.x[[end]], DAT.y[[end]], DAT.z[[end]], color = 1)
    #     plot!(prdt.x, prdt.y, prdt.z, color = 2); scatter!(prdt.x[[end]], prdt.y[[end]], prdt.z[[end]], color = 2)
    plot(pxy, size = (800, 800), plot_title = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]")
    png("figure/$(DAT.ISO3[1]).png")
end
DAT = data[data.ISO3 .== "KOR", :]
plot(plot(DAT.t, DAT.x, xlabel = "t", ylabel = "x", color = :black), 
     plot(DAT.t, DAT.y, xlabel = "t", ylabel = "y", color = :black), 
     plot(DAT.t, DAT.z, xlabel = "t", ylabel = "z", color = :black), 
     layout = (3, 1), size = (500, 500), legend = :none,
     plot_title = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]")
png("temp.png")
extrema(data.t)