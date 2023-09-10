include("dataload.jl")
include("backend_GR.jl")

임시데프 = CSV.read("G:/recent/rslt 1009.csv", DataFrame)

default(plot_titlefontfamily = "맑은 고딕")
for k in 1:17
    aaaa = marginby(groupby(rslt, :location)[k], :gender).y2050
    bbbb = marginby(groupby(임시데프, :location)[k], :gender).y2050

    plot(title = 이름_지역[k])
    plot!(a, aaaa, lw = 2, label = "Control")
    plot!(a, bbbb, lw = 2, label = "x2Gravity")
    png("G:/figure/49 $(이름_지역[k]).png")
end