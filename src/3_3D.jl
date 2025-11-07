include("../../DataDrivenModel/core/header.jl")
include("0_datacall.jl")

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                    This is 3D analysis

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
vrbl = [:dx, :dy, :dz], [:x, :y, :z]

dropmissing!(data)
add_fold!(data, seed = 0)
TD_ = [data[data.fold .!= k, :] for k in sort(unique(data.fold))]
VD_ = [data[data.fold .== k, :] for k in sort(unique(data.fold))]


# Taic = []
# Vaic = []
# rankAnalysis = []
# d_ = 1:6
# for d = d_
#     cnfg = cook(last(vrbl), poly = 0:d)
#     ΘX = Θ(data[:, last(vrbl)], cnfg)
#     push!(rankAnalysis, [d, rank(ΘX), size(ΘX, 2)])

#     f_ = [SINDy(TD, vrbl, cnfg, λ = 0) for TD in TD_]
#     push!(Taic, mean([f.aic for f in f_]))
#     push!(Vaic, mean([AIC(f, VD) for (f, VD) in zip(f_, VD_)]))
# end
# rankAnalysis
# plot(xticks = d_)
# # plot!(d_, Taic, label = "training AIC")
# plot!(d_, Vaic, label = "AIC")

# SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)])) |> print

# λ_ = exp10.(-4:0.1:-1)
# λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]); λ).mse for λ in λ_]
# plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:.5:-1), legend = :none, color = :black)

# f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]), λ = exp10(-3.0))
f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:2), λ = exp10(-3.0))
f |> print

# pxyargs = (; xlabel = "x", ylabel = "y", lims = [0, 6], legend = :topright)
# pxzargs = (; xlabel = "x", ylabel = "z", lims = [0, 6], legend = :none)
# pyzargs = (; xlabel = "y", ylabel = "z", lims = [0, 6], legend = :none)
# xyzargs = (; xlabel = "x", ylabel = "y", lims = [0, 6], zlabel = "z", legend = :none)
# assessment = []
# for DAT in groupby(data, :ISO3)
#     tend = (DAT.t[end] - DAT.t[1])
#     prdt = DataFrame(solve(f, [DAT[1, last(vrbl)]...], 0:0.01:tend), last(vrbl))
#     push!(assessment, prdt.z[end])
#     pt = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]"
#     @info "$pt"

#     pxy=plot(; pxyargs...)
#         plot!(DAT.x, DAT.y, color = 1, label = :none); scatter!(DAT.x[[end]], DAT.y[[end]], color = 1, label = "data")
#         plot!(prdt.x, prdt.y, color = 2, label = :none); scatter!(prdt.x[[end]], prdt.y[[end]], color = 2, label = "model")
#     pyz=plot(; pyzargs...)
#         plot!(DAT.y, DAT.z, color = 1); scatter!(DAT.y[[end]], DAT.z[[end]], color = 1)
#         plot!(prdt.y, prdt.z, color = 2); scatter!(prdt.y[[end]], prdt.z[[end]], color = 2)
#     pxz=plot(; pxzargs...)
#         plot!(DAT.x, DAT.z, color = 1); scatter!(DAT.x[[end]], DAT.z[[end]], color = 1)
#         plot!(prdt.x, prdt.z, color = 2); scatter!(prdt.x[[end]], prdt.z[[end]], color = 2)
#     xyz=plot(; xyzargs...)
#         plot!(DAT.x, DAT.y, DAT.z, color = 1); scatter!(DAT.x[[end]], DAT.y[[end]], DAT.z[[end]], color = 1)
#         plot!(prdt.x, prdt.y, prdt.z, color = 2); scatter!(prdt.x[[end]], prdt.y[[end]], prdt.z[[end]], color = 2)
#     plot(pxy, pyz, pxz, xyz, size = (800, 800), plot_title = pt)
#     png("figure/$(DAT.ISO3[1]).png")
# end

ic = [[dr...] for dr in eachrow(combine(groupby(data, :ISO3), row -> last(row))[:, [:x, :y, :z]])]
maxk = length(ic)
tspan = 0:0.01:50
trajx = zeros(maxk, length(tspan))
trajy = zeros(maxk, length(tspan))
trajz = zeros(maxk, length(tspan))
bit_valid = zeros(Bool, maxk)
# if isfile("data/trajx3.csv") && isfile("data/trajy3.csv") && isfile("data/trajz3.csv")
#     trajx = Matrix(CSV.read("data/trajx3.csv", DataFrame))'
#     trajy = Matrix(CSV.read("data/trajy3.csv", DataFrame))'
#     trajz = Matrix(CSV.read("data/trajz3.csv", DataFrame))'
# else
    @showprogress for k in 1:maxk
        # if x + 1 < y continue end
        traj = solve(f, ic[k], tspan)
        if all(0 .< traj[end, :])
            trajx[k, :] .= traj[:, 1]
            trajy[k, :] .= traj[:, 2]
            trajz[k, :] .= traj[:, 3]
            bit_valid[k] = true
        end
    end
    trajx = trajx[bit_valid, :]
    trajy = trajy[bit_valid, :]
    trajz = trajz[bit_valid, :]
    CSV.write("data/trajx3.csv", DataFrame(trajx', :auto))
    CSV.write("data/trajy3.csv", DataFrame(trajy', :auto))
    CSV.write("data/trajz3.csv", DataFrame(trajz', :auto))
# end

frame = plot(xlabel = "x", ylabel = "y", zlabel = "z", legend = :none, size = [600, 600], xlims = [0, 6.5], ylims = [0, 6.5])

plt_vf = deepcopy(frame)
for k in axes(trajx, 1)
    plot!(plt_vf, trajx[k, :], trajy[k, :], trajz[k, :], color = :black, label = :none)
    scatter!(plt_vf, trajx[k, [end]], trajy[k, [end]], trajz[k, [end]], color = :black, ms = 1, alpha = .5)
end
# plt_vf
savefig("G:/population/pred50y.html")

Plots.plotly()
Plots.gr()

plot(vec(sum(trajz, dims = 1))/212)
plot(eachrow(trajz))