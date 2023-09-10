include("dataload.jl")
include("Backend_GR.jl")

using Plots
using CSV, DataFrames
using Random
default(titleloc = :right)

DATA = CSV.read("data/KOSIS/인구동태.csv", DataFrame, missingstring = "-")

TFRplot = plot(title = "a", xlims = (1990, 2050), ylabel = "TFR", yticks = [0,1,2], ylims = (0, 2), legend = :none)
plot!(TFRplot, DATA.시점, DATA[:, 2], color = :black, lw = 2)
Random.seed!(1)
for k in 1:5
    plot!(TFRplot, 2021:2050, 0.808 .+ cumsum([0, .02randn(29)...] .+ 0.005(k-3)), lw = 2, color = :black, ls = :dash)
end
annotate!(TFRplot, 2040, 1.5, text("?", :bold, 28, :black))
plot!(TFRplot, left_margin = 4mm)

DEADplot = plot(title = "b")
plot!(DEADplot, 1990:2021, DATA[21:52, 6], color = :black, ylabel = "Life Expectancy", ylims = (70, 85), xlims = (1990, 2020), legend = :none)
plot!(twinx(DEADplot), 1990:2021, DATA[21:52, 5], color = :blue, ylabel = "Death Rate", ylims = (450, 650), xlims = (1990, 2020), legend = :none, y_foreground_color_border=:blue, y_foreground_color_text=:blue, y_guidefontcolor =:blue)
plot!(DEADplot)

MGRNplot = plot(title = "c")
plot!(MGRNplot, 1990:2021, DATA[21:52, 4], color = :black, legend = :none, ylabel = "Domestic migration", xlims = (1990, 2020), yformatter = x -> x/(10^6))
annotate!(MGRNplot, 1990, 9.8 * 10^6, text(L"\times 10^6", :left, 10, :black))

plot(TFRplot, DEADplot, MGRNplot, size = (800, 400), right_margin = 4mm, layout = @layout [a [b ; c]])
png("G:/figure/22_1.png")

# x=rand(10)
# p = plot(x,xlabel="x axis",ylabel="y axis", lc=:blue, lw=3)
# plot!(x_guidefontcolor=:red, y_guidefontcolor=:green)
# plot!(x_foreground_color_axis=:red, y_foreground_color_axis=:green)
# plot!(x_foreground_color_text=:red, y_foreground_color_text=:green)


# MGRN = CSV.read("data/KOSIS/mobility.csv", DataFrame)

# ingdf = groupby(MGRN, :전출지별); names()
# outgdf = groupby(MGRN, :전입지별); names()

# MGRNlot = plot(legend = :none)
# for k in 1:17
#     plot!(MGRNlot, 2012:2021,
#     vec(sum(Matrix(ingdf[k][:, Not(1:16)] .- outgdf[k][:, Not(1:16)]), dims = 1))
#     , color = 무지개[k])
# end
# MGRNlot

# plot(rand(10))
# plot!(twinx(), rand(10), axis = :right)
