include("dataload.jl")
include("Backend_GR.jl")

KLD(p,q) = -sum(p .* log.(q ./ p))

l = @layout [[a; b] c{0.2w}]

가구 = vldn_[2014]

ptemp__ = []
dtemp__ = []
ltemp__ = []
ptemp2 = plot(legend = :none, title = "Kullback-Leibler Divergence")
ptemp3 = plot(legend = :none, title = "L1 Error")
# lagendtrick = plot(ticks = false,  lims = (-1,0))
lagendtrick = plot(ticks = false, showaxis = false, grid = false, lims = (-1,0), legend = :outerleft, right_margin = -10mm, legend_foreground_color = false)
for i = 1:17
    ptemp_ = []
    dtemp_ = []
    ltemp_ = []
    for j = 1:10
        px = trunc.(Int64, real[(1:100) .+ 100(2i-1), j + 15] .+ real[(1:100) .+ 100(2i-2), j + 15])
        qx =               vldn[(1:100) .+ 100(2i-1), j + 15] .+ vldn[(1:100) .+ 100(2i-2), j + 15]
        rx =               가구[(1:100) .+ 100(2i-1), j + 15] .+ 가구[(1:100) .+ 100(2i-2), j + 15]
        ptemp = plot(ticks = false, grid = false, xlims = (0, 99), ylims = (0, Inf), legend = :none)
        plot!(ptemp, 0:99, px, alpha = 0.5, lw = 2, color = :black)
        # plot!(ptemp, 0:99, qx, alpha = 0.5, lw = 2, color = :red)
        if (j % 5) == 3 title!(ptemp, name_location[i]) end
        if j > 2
            plot!(ptemp, 0:99, rx, alpha = 0.5, lw = 2, color = :blue)
        end
        push!(ltemp_, sum(abs, px - qx))
        px = px ./ sum(px)
        qx = qx ./ sum(qx)
        kld = KLD(px, qx)
        # annotate!(ptemp, 30, 0, "$(trunc(kld, digits = 4))")
        push!(dtemp_, kld)
        push!(ptemp_, ptemp)
    end
    push!(ptemp__, ptemp_)
    plot!(ptemp2, 2012:2021, dtemp_, color = 무지개[i], label = name_location[i])
    plot!(ptemp3, 2012:2021, ltemp_, color = 무지개[i], label = name_location[i])
    scatter!(lagendtrick, 1:2, 1:2, color = 무지개[i], label = name_location[i], shape = :circ, msw = 0)
end

pcplt = plot(vcat(ptemp__...)..., layout = (17,10), size = 250 .* (10, 14), dpi = 200);
pcplt2 = plot(ptemp3, ptemp2, lagendtrick, layout = l, size = (800, 400))

savefig(pcplt, "G:/figure/19 00.pdf")
savefig(pcplt2, "G:/figure/19 01.svg")

plot(
  plot(ptemp__[8][1], xlabel = "2012")
, plot(ptemp__[8][3], xlabel = "2014")
, plot(ptemp__[8][6], xlabel = "2018")
, plot(ptemp__[8][10], xlabel = "2021")
, layout = (4,1), size = (400, 800), title = "")

# ---

ts_entropy = zeros(17, 9)
rslt_tnsr = reshape(rslt.y2021, 100, :)
rslt_dstr = Float64.(rslt_tnsr[:, 1:2:end] + rslt_tnsr[:, 2:2:end])
for j in 1:17
    rslt_dstr[:,j] ./= sum(rslt_dstr, dims = 1)[1,j]
end

for t in 1:9
    vldn_tnsr = reshape(vldn_[2011+t].y2021, 100, :)
    vldn_dstr = Float64.(vldn_tnsr[:, 1:2:end] + vldn_tnsr[:, 2:2:end])
    temp = sum(vldn_dstr, dims = 1)
    for j in 1:17
        vldn_dstr[:,j] ./= temp[1,j]
        ts_entropy[j, t] = KLD(rslt_dstr[:,j], vldn_dstr[:,j])
    end
end

using StatsBase
평균 = vec(mean(ts_entropy, dims = 1))
표편 = vec(std(ts_entropy, dims = 1))

pcplt31 = plot(palette = 무지개, ylabel = "KLD", legend = :none, ylims = (0, 0.025))
pcplt32 = plot(palette = 무지개, xlabel = "Initial point", ylabel = "KLD", legend = :none, yscale = :log10)
for i in 1:17
    packet = (; label = name_location[i], msw = 0, shape = :rect, st = :line)
    scatter!(pcplt31, 2012:2020, ts_entropy[i, :]; packet...)
    scatter!(pcplt32, 2012:2020, ts_entropy[i, :]; packet...)
end

# ---

pure_mgrn = mgrn[mgrn.from .!= mgrn.to,:]
birth = rslt[rslt.age .== 0,:]

pcplt33 = plot(palette = 무지개, xlabel = "Year", ylabel = "출생대비유입", legend = :none)
for k in 1:17
    유입 = ts_sum(select(groupby(pure_mgrn, :to)[k], Not(2)))
    유출 = ts_sum(select(groupby(pure_mgrn, :from)[k], Not(2)))
    인구 = ts_sum(groupby(rslt, :location)[k])
    출생 = ts_sum(groupby(birth, :location)[k])
    사망 = ts_sum(groupby(dead, :location)[k])
    인구대비이동 = ((유입 + 유출)  ./ 인구)[1:9]
    성장대비이동 = ((유입)  ./ (출생))[1:9]
    # plot!(pcplt33, 2012:2021, 인구대비이동)
    plot!(pcplt33, 2012:2020, 성장대비이동, lw = 2)
end

pcplt3 = plot(pcplt31, pcplt32, pcplt33, layout = (3,1), size = (800, 800))
savefig(pcplt3, "G:/figure/19 03.svg")

plot(1 ./ rand(10000), st = :scatterhist, scale = :log10, ylims = (0.001, 100000))