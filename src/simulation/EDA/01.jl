include("basic.jl")


df = deepcopy(rslt)
age_l = filter(:age => age ->      age < 15, df) |> ts_sum
age_m = filter(:age => age -> 15 ≤ age < 65, df) |> ts_sum
age_h = filter(:age => age -> 65 ≤ age     , df) |> ts_sum
age_c = cumsum([age_l age_m age_h], dims = 2)
age_r = 100 * cumsum([(age_l ./ age_m) (age_h ./ age_m)], dims = 2)
age_cr = 100 * (age_c ./ age_c[:, 3])
t_super_aging = 2020 + findfirst((age_h ./ age_c[:,3]) .> .21)

pp2 = plot(legend = :topright,
    xlims = (2021, 2070), ylims = (0,Inf), xticks = [2021, (2030:10:yend)...],
    ylabel = "Population")
plot!(pp2, 2021:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
plot!(pp2, 2021:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
plot!(pp2, 2021:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
# annotate!(pp2, [2022], [52], [(L"\times 10^6", 8)])
png(pp2, "D:/figure/subfigure/01 0 korea pp2.png")

pp4 = plot(legend = :topleft,
    xlims = (2021, 2070), ylims = (0,150), xticks = [2021, (2030:10:yend)...],
    ylabel = "Ratio(%)")
plot!(pp4, 2021:yend, age_r[:,2], color =    red, fa = .5, fillrange = age_r[:,1], label = "Elder Dependency Ratio")
plot!(pp4, 2021:yend, age_r[:,1], color = yellow, fa = .5, fillrange = 0         , label = "Child Dependency Ratio")
png(pp4, "D:/figure/subfigure/01 0 korea pp4.png")

pp5 = plot(legend = :none,
    xlims = (2021, 2070), ylims = (0,100), xticks = [2021, (2030:10:yend)...],
    xlabel = "Year", ylabel = "Ratio(%)")
plot!(pp5, 2021:yend, age_cr[:,3], color =    red, fa = .5, fillrange = age_cr[:,2], label = "65-")
plot!(pp5, 2021:yend, age_cr[:,2], color = orange, fa = .5, fillrange = age_cr[:,1], label = "15-64")
plot!(pp5, 2021:yend, age_cr[:,1], color = yellow, fa = .5, fillrange = 0          , label = "-14")
vline!([t_super_aging], style = :dash, color = :black, label = "super aging")
png(pp5, "D:/figure/subfigure/01 0 korea pp5.png")

pp6 = plot(legend = :bottomleft,
    xlims = (2021, 2070), xticks = [2021, (2030:10:yend)...],
    xlabel = "Year", ylabel = "Number")
plot!(pp6, 2021:yend, filter(:age => age -> age == 0, rslt) |> ts_sum, fa = .5, fillrange = 0, color = :green, label = "Birth")
plot!(pp6, 2021:yend, -(dead |> ts_sum),                               fa = .5, fillrange = 0, color = :black, label = "Death")
png(pp6, "D:/figure/subfigure/01 0 korea pp6.png")

plot(pp2, pp4, pp5, pp6, layout = (2,2), size = 60 .* (16, 9), plot_title = "Korea",
    left_margin = 3Plots.mm, right_margin = 3Plots.mm, dpi = 300)
png("D:/figure/01 0 korea.png")
# savefig("D:/figure/svg/01 0 korea.svg")