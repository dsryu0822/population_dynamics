include("../../DataDrivenModel/core/header.jl")
include("../src/0_datacall.jl")

vrbl = [:dx, :dy], [:x, :y]

data = data[:, Not([:z, :dz])]
dropmissing!(data)
add_fold!(data, seed = 0)
data_ = Dict(unique(data.ISO3) .=> collect(groupby(data, :ISO3)))
TD_ = [data[data.fold .!= k, :] for k in sort(unique(data.fold))]
VD_ = [data[data.fold .== k, :] for k in sort(unique(data.fold))]

pythonplot()
plt_rf = plot(lims = [1, 6])
for DAT in groupby(data, :ISO3)
    uc = first(DAT.fold)
    plot!(plt_rf, DAT.x, DAT.y, color = uc, label = :none, alpha = .5, arrow = Plots.arrow(:open, :head, .6, .6)) #, text = DAT.ISO3[1])
    # scatter!(plt_rf, DAT.x[[end]], DAT.y[[end]], color = uc, label = "data", msw = 0, alpha = .5) #, text = DAT.ISO3[1])
end
plot(plt_rf, size = [400, 400])

trajx = Matrix(CSV.read("data/trajx2.csv", DataFrame))'
trajy = Matrix(CSV.read("data/trajy2.csv", DataFrame))'
plt_vf = plot(lims = [1, 6])
for k in axes(trajx, 1)
    plot!(plt_vf, trajx[k, :], trajy[k, :], color = :black, label = :none, alpha = .5, arrow = Plots.arrow(:head, .6, .6))
    plot!(plt_vf, trajx[k, :], trajy[k, :], color = :black, label = :none, alpha = .5, arrow = Plots.arrow(:head, .6, .6))
    # scatter!(plt_vf, trajx[k, [end]], trajy[k, [end]], color = :black, ms = 3, alpha = .5)
end

x_ = 0:.005:7
y_ = 0:.005:7
nc_x = Bool.(Matrix(CSV.read("data/nc_x.csv", DataFrame)))
nc_y = Bool.(Matrix(CSV.read("data/nc_y.csv", DataFrame)))
pt_nc_x = collect.(getproperty.(findall(nc_x), :I))
xy_ = [[50, 0]]
while pt_nc_x |> !isempty
    temp = pt_nc_x .- Ref(xy_[end])
    push!(xy_, popat!(pt_nc_x, argmin(sum.(abs2, temp))))
end
popfirst!(xy_)
xXY_ = 1000exp10.(stack(xy_) / 200)
pt_nc_y = collect.(getproperty.(findall(nc_y), :I))
xy_ = [[50, 0]]
while pt_nc_y |> !isempty
    temp = pt_nc_y .- Ref(xy_[end])
    push!(xy_, popat!(pt_nc_y, argmin(sum.(abs2, temp))))
end
popfirst!(xy_)
yXY_ = 1000exp10.(stack(xy_) / 200)

plt_nc = plot()
scatter!(plt_nc, xXY_[1, :], xXY_[2, :], label = :none, ms = 1, size = [600, 600], color = :gold, msw = 0)
scatter!(plt_nc, yXY_[1, :], yXY_[2, :], label = :none, ms = 1, size = [600, 600], color = :brown, msw = 0)
plot!(plt_nc, scale = :log10, ticks = exp10.(3:2:9), lims = exp10.([3, 9]))

fp = [4.058707750181835, 4.48902300450415]
scatter!(plt_nc, [1000exp10(fp[1])], [1000exp10(fp[2])], color = :black, label = :none, ms = 10)

plot(plt_rf, plt_vf, plt_nc, ticks = [], layout = (1, 3), size = 2*[900, 300], left_margin = 1cm, bottom_margin = 1cm)
png("temp")
