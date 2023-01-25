vldn_rslt = CSV.read("G:/recent/vldn/rslt 0001.csv", DataFrame)

pp09_01 = plot(title = "Total population of Korea", ylabel = "Population", xticks = [2012, 2017, 2021])
plot!(pp09_01, 2012:2021, ts_sum(rslt)[1:10], label = "Real Data", color = :black, style = :solid, lw = 2)
plot!(pp09_01, 2012:2021, ts_sum(vldn_rslt), label = "Simulation", color = red   , style = :dash , lw = 2)

pp09_02 = plot(title = L"L_1" * " error", xlabel = "year", xticks = [2012, 2017, 2021])
plot!(pp09_02, 2012:2021, abs.(ts_sum(rslt)[1:10] - ts_sum(vldn_rslt)), legend = :none,
ylabel = "Number (#)", color = red, lw = 2, fillrange = 0, fa = 0.25)
plot!(twinx(), 2012:2021, 100abs.(ts_sum(rslt)[1:10] - ts_sum(vldn_rslt)) ./ ts_sum(rslt)[1:10], legend = :none,
ylabel = "Percentage (%)", color = red, lw = 2, fillrange = 0, fa = 0.25)

plot(pp09_01, pp09_02, layout = (2, 1), size = (800, 400),
margin = 5mm, right_margin = 5mm, bottom_margin = 5mm, dpi = 200)
png("G:/figure/09.png")