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
scatter!(plt_nc, xXY_[1, :], xXY_[2, :], label = :none, ms = .5, size = [600, 600], alpha = .5, color = :gold, msw = 0)
scatter!(plt_nc, yXY_[1, :], yXY_[2, :], label = :none, ms = .5, size = [600, 600], alpha = .5, color = :brown, msw = 0)
plot!(plt_nc, scale = :log10, ticks = exp10.(3:2:9), lims = exp10.([3, 9]))

fp = [4.058707750181835, 4.48902300450415]
scatter!(plt_nc, [1000exp10(fp[1])], [1000exp10(fp[2])], color = :black, label = :none, ms = 5)

plot(plt_rf, plt_vf, plt_nc, ticks = [], layout = (1, 3), size = [900, 300], left_margin = 1cm, bottom_margin = 1cm)
png("temp")


plt_compare = deepcopy(plt_nc)
class = []
for DAT = groupby(data, :ISO3)
    xend = DAT.x[end]
    yend = DAT.y[end]
    dx, dy = f([xend, yend])
    push!(class, dx*dy < 0 ? 2 : yend > xend ? 1 : 3)
end
data_ecnm = combine(groupby(data, :ISO3), row -> last(row))
data_ecnm[!, :class] = class
plt_compare = deepcopy(plt_nc)
scatter!(plt_compare, 1000exp10.(data_ecnm.x[data_ecnm.class .== 1]), 1000exp10.(data_ecnm.y[data_ecnm.class .== 1]), color = 1, shape = :rect, msc = :white, ms = 3)
scatter!(plt_compare, 1000exp10.(data_ecnm.x[data_ecnm.class .== 2]), 1000exp10.(data_ecnm.y[data_ecnm.class .== 2]), color = 2, shape = :rect, msc = :white, ms = 3)
scatter!(plt_compare, 1000exp10.(data_ecnm.x[data_ecnm.class .== 3]), 1000exp10.(data_ecnm.y[data_ecnm.class .== 3]), color = 3, shape = :rect, msc = :white, ms = 3)
scatter(plt_compare, 1000exp10.(data_ecnm.x), 1000exp10.(data_ecnm.y), text = data_ecnm.ISO3, size = [2000, 2000])
png("temp")


clr_ = [2, 1, 2, 2, 3, 3]
tgt_ = ["JPN", "PRI", "GUM", "CHN", "NGA", "QAT"]
for k = 1:6
    plot(data_[tgt_[k]].x, data_[tgt_[k]].y, size = [200, 200], color = clr_[k], ticks = [], framestyle = :none)
    scatter!([data_[tgt_[k]].x[end]], [data_[tgt_[k]].y[end]], color = clr_[k], ms = 5, msw = 0)
    png("temp_$(tgt_[k])")
end