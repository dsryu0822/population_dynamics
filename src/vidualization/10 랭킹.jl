include("dataload.jl")
include("Backend_GR.jl")
include("04.jl")

using StatsBase
default()
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
function arr_ranking(arrarr)
    return [mat for mat ∈ eachrow(hcat(ordinalrank.(eachrow(hcat(arrarr...)), rev = true)...))]
end

function myrankingplot(arrarr, titlestr)
    AR = arr_ranking(arrarr)

    p1 = plot(yflip = true, showaxis = false, legend = :none, grid = false, xlims = (2002,yend+2),
    title = titlestr, size = (1000,600))
    for k in 1:17
        plot!(p1, ybgn:yend, AR[k], lw = 4, color = 무지개[AR[k][1]])
        plot!(p1, rectangle(10, 1, 2002, AR[k][1] - 1/2), color = 무지개[AR[k][1]], linealpha = 0)
        annotate!(p1, 2004.5, AR[k][1], text(name_location[k], :white, halign = :left))
        annotate!(p1, 2004, AR[k][1], text(AR[k][1], :white, halign = :right))
        등락 = AR[k][1] - AR[k][end]
        등락_문자 = ""
        if 등락 > 0
            등락_문자 = "▲ $등락"
        elseif 등락 < 0
            등락_문자 = "▼ $(-등락)"
        else
            등락_문자 = "   0"
        end
        annotate!(p1, yend + .5, AR[k][end], text(등락_문자, 무지개[AR[k][1]], halign = :left, pointsize = 10))
    end
    annotate!(p1, 2012, 18, 2012)
    annotate!(p1, 2030, 18, 2030)
    annotate!(p1, 2050, 18, 2050)

    return p1
end

p10_01 = myrankingplot(netgrowth_, "Net Growth Ranking")
p10_02 = myrankingplot(netmigrat_, "Net Migration Ranking")
p10_03 = myrankingplot([ts_sum(gdf) for gdf in groupby(rslt, :location)], "Population Ranking")
savefig(p10_01, "G:/figure/10_01.svg")
savefig(p10_02, "G:/figure/10_02.svg")
savefig(p10_03, "G:/figure/10_03.svg")