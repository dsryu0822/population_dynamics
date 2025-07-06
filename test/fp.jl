using ProgressMeter
packages = [:CSV, :DataFrames, :Plots, :Colors, :ColorSchemes]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm

shape_ = Dict("stbl" => :circle, "sddl" => :square, "unst" => :star)

info_ = CSV.read.(["data/info_fp_$N.csv" for N in 1:12], DataFrame)
for n in 1:12
    info_[n][!,:shape] = replace(info_[n].stability, shape_...)
end

info_[12]

info_[6]
color_ = palette(:darktest, 12)
p1 = plot(lims = [0, .7], size = [800, 800])
for n in 1:12
    p1 = scatter(p1, info_[n].x, info_[n].y, color = color_[n], shape = info_[n].shape, msw = 0, alpha = 0.5, ms = 10)
end
p1

setN = 6
p2 = plot(xlims = [0, .7], ylims = [0, .7], size = [800, 800])
for setN in 1:11
    temp = CSV.read("G:/population/N=$setN/basin.csv", DataFrame)
    temp = temp[(0 .< temp.xend .< 0.5) .&& (0 .< temp.yend .< 0.5), :]
    p3 = scatter(p2, temp.xend, temp.yend, color = color_[setN], msw = 0)
    png(p3, "G:/population/N=$setN/endpoints.csv")
end
p2