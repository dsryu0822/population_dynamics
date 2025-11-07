include("../../DataDrivenModel/core/header.jl")
include("0_datacall.jl")

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                    This is 2D analysis

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
vrbl = [:dx, :dy], [:x, :y]

data = data[:, Not([:z, :dz])]
dropmissing!(data)
add_fold!(data, seed = 0)
data_ = Dict(unique(data.ISO3) .=> collect(groupby(data, :ISO3)))
TD_ = [data[data.fold .!= k, :] for k in sort(unique(data.fold))]
VD_ = [data[data.fold .== k, :] for k in sort(unique(data.fold))]


####################### parameter selection #######################
# Taic = []
# Vaic = []
# rankAnalysis = []
# d_ = 1:10
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

# λ_ = exp10.(-4:0.1:0)
# λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]); λ).mse for λ in λ_]
# plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:.5:0), legend = :none, color = :black)

# f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]), λ = exp10(-2.5))
# f |> print

f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:6), λ = exp10(-2.5))
export_tex(f, digits = 3) |> print
f |> print

####################### trajectory for each ISO3 #######################
# pxyargs = (; xlabel = "x", ylabel = "y", lims = [0, 6], legend = :topright)
# for DAT in groupby(data, :ISO3)[1:1:end]
#     tend = (DAT.t[end] - DAT.t[1])
#     prdt = DataFrame(solve(f, [DAT[1, last(vrbl)]...], 0:0.01:tend), last(vrbl))
#     pt = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]"
#     @info "$pt"

#     pxy=plot(; pxyargs...)
#         plot!(DAT.x, DAT.y, color = 1, label = :none); scatter!(DAT.x[[end]], DAT.y[[end]], color = 1, label = "data")
#         plot!(prdt.x, prdt.y, color = 2, label = :none); scatter!(prdt.x[[end]], prdt.y[[end]], color = 2, label = "model")
#     plot(pxy, size = (800, 800), plot_title = pt)
#     png("figure/$(DAT.ISO3[1]).png")
# end
# DAT = data[data.ISO3 .== "KOR", :]
# plot(plot(DAT.t, DAT.x, xlabel = "t", ylabel = "x", color = :black), 
#      plot(DAT.t, DAT.y, xlabel = "t", ylabel = "y", color = :black), 
#      layout = (2, 1), size = (500, 500), legend = :none,
#      plot_title = "$(DAT.ISO3[1]): $(iso3[DAT.ISO3[1]]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]")
# png("temp.png")


####################### visualization #######################
emptyframe = plot(xlabel = L"x", ylabel = L"y", legend = :none, size = [600, 600], lims = [0, 6.5])


####################### real field #######################
plt_rf = deepcopy(emptyframe)
for DAT in groupby(data, :ISO3)
    uc = first(DAT.fold)
    plot!(plt_rf, DAT.x, DAT.y, color = uc, label = :none, alpha = .5)
    scatter!(plt_rf, DAT.x[[end]], DAT.y[[end]], color = uc, label = "data", msw = 0, alpha = .5, text = DAT.ISO3[1])
end
plot(plt_rf, size = [2000, 2000]);
png("rf.png")

####################### vector field #######################
xy = Base.product(0:.1:6, 0:.1:6)
trajx = zeros(length(xy), 2001)
trajy = zeros(length(xy), 2001)
bit_valid = zeros(Bool, length(xy))
if isfile("data/trajx2.csv") && isfile("data/trajy2.csv")
    trajx = Matrix(CSV.read("data/trajx2.csv", DataFrame))'
    trajy = Matrix(CSV.read("data/trajy2.csv", DataFrame))'
else
    @showprogress for (k, (x,y)) in enumerate(xy)
        if x + 1 < y continue end
        traj = solve(f, [x, y], 0:0.01:20)
        if all(0 .< traj[end, :] .< 6)
            trajx[k, :] .= traj[:, 1]
            trajy[k, :] .= traj[:, 2]
            bit_valid[k] = true
        end
    end
    trajx = trajx[bit_valid, :]
    trajy = trajy[bit_valid, :]
    CSV.write("data/trajx2.csv", DataFrame(trajx', :auto))
    CSV.write("data/trajy2.csv", DataFrame(trajy', :auto))
end

plt_vf = deepcopy(emptyframe)
for k in axes(trajx, 1)
    plot!(plt_vf, trajx[k, :], trajy[k, :], color = :black, label = :none)
    scatter!(plt_vf, trajx[k, [end]], trajy[k, [end]], color = :black, ms = 3, alpha = .5)
end
png("vf")

####################### method of nullcline #######################
x_ = 0:.005:7
y_ = 0:.005:7
xy = [[x, y] for (x, y) in Base.product(x_, y_)]
@time fxy = f.(xy) # about 20sec
dx = first.(fxy)
dy = last.(fxy)

boundary = ones(Bool, size(dx)...); boundary[1:end, 1] .= false; boundary[1:end, end] .= false; boundary[1, 1:end] .= false; boundary[end, 1:end] .= false
# nc_x = (((circshift(dx, (1, 0)) .* circshift(dx, (-1, 0))) .< 0) .|| ((circshift(dx, (0, 1)) .* circshift(dx, (0, -1))) .< 0)) .&& boundary
# nc_y = (((circshift(dy, (1, 0)) .* circshift(dy, (-1, 0))) .< 0) .|| ((circshift(dy, (0, 1)) .* circshift(dy, (0, -1))) .< 0)) .&& boundary
nc_x = (((circshift(dx, (1, 0)) .* circshift(dx, (-1, 0))) .< 0) .|| ((circshift(dx, (0, 1)) .* circshift(dx, (0, -1))) .< 0)) .&& boundary
nc_y = (((circshift(dy, (1, 0)) .* circshift(dy, (-1, 0))) .< 0) .|| ((circshift(dy, (0, 1)) .* circshift(dy, (0, -1))) .< 0)) .&& boundary
CSV.write("data/nc_x.csv", DataFrame(Int64.(nc_x), :auto))
CSV.write("data/nc_y.csv", DataFrame(Int64.(nc_y), :auto))

plt_nc = deepcopy(emptyframe)
heatmap!(plt_nc, x_, y_, nc_x', colorbar = :none, transpose = true, color = [:transparent, :red])
heatmap!(plt_nc, x_, y_, nc_y', colorbar = :none, transpose = true, color = [:transparent, :blue])
png("nc")

# ####################### fixed point #######################
# cadny = [xy[ij] for ij in findall(nc_x .&& nc_y)]
# scatter(plt_vf, first.(cadny), last.(cadny), msw = 0)
# @time futurecandy = [solve(f, cdy, 0:0.01:10)[end, :] for cdy in cadny]
# @time futurecandy = [solve(f, cdy, 0:0.01:800)[end, :] for cdy in futurecandy]

# zoom = plot(plt_vf, xlims = [3.5, 4.5], ylims = [4, 5])
# scatter(zoom, first.(cadny), last.(cadny))
# scatter(zoom, first.(futurecandy), last.(futurecandy))

# # fp = futurecandy[argmin(sum.(abs2, f.(futurecandy)))]
# # fp = [4.058707750181835, 4.48902300450415]
# # 6.894396291195874e-15

# Df = jacobian(Function, f)
# Df(fp)
# using LinearAlgebra
# eigen(Df(fp)) |> propertynames
# plot([0, eigen(Df(fp)).vectors[1, 1]], [0, eigen(Df(fp)).vectors[2, 1]])
# plot([0, eigen(Df(fp)).vectors[1, 2]], [0, eigen(Df(fp)).vectors[2, 2]], lims = [0, 6], size = [600, 600], xlabel = "x", ylabel = "y", legend = :none)

# eigen(Df(fp)).vectors[:, 2]

# line(x) = (0.6766074657953228 / 0.7363438987524994)*x + 4.48902300450415 - 4.058707750181835*(0.6766074657953228 / 0.7363438987524994)
# plot(vf, x_, line.(x_), color = :red, lw = 2)

####################### classification #######################
include("0_datacall.jl")
dropmissing!(data)
add_fold!(data, seed = 0)

class = []
for DAT = groupby(data, :ISO3)
    xend = DAT.x[end]
    yend = DAT.y[end]
    dx, dy = f([xend, yend])
    push!(class, dx*dy < 0 ? 2 : yend > xend ? 1 : 3)
end
data_ecnm = combine(groupby(data, :ISO3), row -> last(row))
data_ecnm[!, :class] = class
scatter(emptyframe, data_ecnm.x, data_ecnm.y, group = data_ecnm.class)
groupby(sort(data_ecnm, :z, rev = true), :class)
png("cs")

plot(legend = :none, ylims = [2, 6], size = [400, 400], xlabel = "class", ylabel = "log GDP per capita")
boxplot!(["A"], data_ecnm.z[data_ecnm.class .== 1])
boxplot!(["B"], data_ecnm.z[data_ecnm.class .== 2])
boxplot!(["C"], data_ecnm.z[data_ecnm.class .== 3])
png("bx")

"""
1. We found a universal trend of world-wide population dynamics.
 - Regardless of the scale of the population, all the countries are lie on the band.
2. For advanced countries, the governing equation describes the super aging behavior.
 - A stable manifold exists, and populous countries are distributed around the fixed point.
3. A classification is possible without the knowledge of the economy.
 - Method of nullcline can suggest the concise criterion for the wealthy level.
"""

distance_fixedpoint = sqrt.(abs2.(filter(row -> row.t == 2022, data).x .- fp[1]) + abs2.(filter(row -> row.t == 2022, data).y .- fp[2]))
asdf = DataFrame(
    ISO3 = unique(data.ISO3)[sortperm(distance_fixedpoint)],
    farfromfx = distance_fixedpoint[sortperm(distance_fixedpoint)]
)
qwer = rename(_ISO3[:, [3, 5, 6]], ["ISO3", "LAT", "LON"])
zxcv = sort(innerjoin(asdf, qwer, on = :ISO3), :farfromfx)

scatter(zxcv[:, 2], zxcv[:, 4])

sort!(data_ecnm, :class)
using GLM
data_class_ = groupby(data_ecnm, :class)
model_linear = []
for df in data_class_
    out = lm(@formula(w ~ x + y), df);
    push!(model_linear, out)
    println("Class: ", df.class[1])
    println("cor(x,w) = ", cor(df.x, df.w))
    println("cor(y,w) = ", cor(df.y, df.w))
    println("cor(x,y) = ", cor(df.x, df.y))
    println("rank = ", rank(Matrix(data_class_[3][:, [:x, :y]])))
    println("R^2 = ", r2(out))
    println("="^10)
end
plot(
    scatter(data_class_[1].w, GLM.predict(model_linear[1], data_class_[1]), color = 1, aspect_ratio = 1, size = [400, 400], legend = :none),
    scatter(data_class_[2].w, GLM.predict(model_linear[2], data_class_[2]), color = 2, aspect_ratio = 1, size = [400, 400], legend = :none),
    scatter(data_class_[3].w, GLM.predict(model_linear[3], data_class_[3]), color = 3, aspect_ratio = 1, size = [400, 400], legend = :none),
    size = [1200, 400], layout = (1, 3), msw = 0
)

proj = DataFrame(solve(f, [data_["JPN"].x[end], data_["JPN"].y[end]], 0:0.01:77)[1:100:end, :], [:x, :y])
insertcols!(proj, 1, :t => (1:nrow(proj)) .+ 2022)
proj[!, :w] = GLM.predict(model_linear[2], proj)


plot(ylims = [3, 6], legend = :none)
plot!(data_["JPN"].t, data_["JPN"].x, color = :red, lw = 2)
plot!(data_["JPN"].t, data_["JPN"].y, color = :blue, lw = 2)
plot!(data_["JPN"].t, data_["JPN"].w, color = :black, lw = 2)
plot!(data_["JPN"].t, GLM.predict(model_linear[2], data_["JPN"]), color = :black, ls = :dash)
plot!(proj.t, proj.x, color = :red, ls = :dash)
plot!(proj.t, proj.y, color = :blue, ls = :dash)
plot!(proj.t, proj.w, color = :black, ls = :dash)
f([proj.x[end], proj.y[end]])