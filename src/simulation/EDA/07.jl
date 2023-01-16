include("basic.jl")

rslt_cptl = filter(:location => is_capital, rslt)
dead_cptl = filter(:location => is_capital, dead)
age_c = cumsum(ts_age(rslt_cptl), dims = 2)

birth = filter(:age => age -> age == 0, rslt_cptl) |> ts_sum
death = (dead_cptl |> ts_sum)

mgrn_cptl_out = filter(:to => !is_capital, filter(:from => is_capital, mgrn))
rename!(mgrn_cptl_out, :from => :location)
select!(mgrn_cptl_out, Not(:to))
mgrn_cptl_in = filter(:to => is_capital, filter(:from => !is_capital, mgrn))
rename!(mgrn_cptl_in, :to => :location)
select!(mgrn_cptl_in, Not(:from))
age_migration = ts_age(mgrn_cptl_in) - ts_age(mgrn_cptl_out)

netmigrat = sum(age_migration, dims = 2)
netgrowth = birth - death

pp07_01 = plot(xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,14,65,99])
plot!(pp07_01, 0:99, df_age.y2012[1:100], label = "2012", lw = 2, color = :black)
plot!(pp07_01, 0:99, df_age.y2046[1:100], label = "2046", lw = 2, color = :navy)
plot!(pp07_01, 0:99, df_age.y2070[1:100], label = "2070", lw = 2, color = :blue)
png(pp07_01, "G:/figure/subfigure/07 0 capital pp07_01.png")

pp07_02 = plot(legend = :topright,
    xlims = (2012, yend), ylims = (0,Inf), xticks = [2012, (2030:10:yend)...],
    ylabel = "Population")
plot!(pp07_02, 2012:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
plot!(pp07_02, 2012:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
plot!(pp07_02, 2012:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
png(pp07_02, "G:/figure/subfigure/07 0 capital pp07_02.png")

pp07_04 = plot(legend = :topright,
    xlims = (2012, yend), xticks = [2012, (2030:10:yend)...],
    ylabel = "Age-specific Migration")
plot!(pp07_04, 2012:yend, age_migration[:,3], color =    red, msw = 0, shape = :circ, lw = 2, fa = 0.1, fillrange = 0, label = "65-")
plot!(pp07_04, 2012:yend, age_migration[:,2], color = orange, msw = 0, shape = :rect, lw = 2, fa = 0.1, fillrange = 0, label = "15-64")
plot!(pp07_04, 2012:yend, age_migration[:,1], color = yellow, msw = 0, shape = :diamond, lw = 2, fa = 0.1, fillrange = 0, label = "-14")
png(pp07_04, "G:/figure/subfigure/07 0 capital pp07_04.png")

pp07_06 = plot(legend = :bottomleft,
    xlims = (2012, yend), xticks = [2012, (2030:10:yend)...],
    ylabel = "Number")
plot!(pp07_06, 2012:yend, birth, fa = .5, fillrange = 0, color = :green, label = "Birth")
plot!(pp07_06, 2012:yend, -death, fa = .5, fillrange = 0, color = :black, label = "Death")
png(pp07_06, "G:/figure/subfigure/07 0 capital pp07_06.png")
    
pp07_08 = plot(ylabel = "Number", xlabel = "Year", xticks = [2012, (2030:10:yend)...], xlims = (2012, yend))
plot!(pp07_08, 2012:yend, netmigrat, lw = 2, fa = .5, fillrange = 0, color = :olive, label = "net migration")
plot!(pp07_08, 2012:yend, netgrowth, lw = 2, fa = .5, fillrange = 0, color = :darkgreen, label = "net growth")

plot(pp07_02, pp07_06, pp07_04, pp07_08, 
    plot_title = "Capital Area", layout = (2,2), size = (1200, 800),
    leftmargin = 5Plots.mm, dpi = 200)
png("G:/figure/07 0 capital.png")