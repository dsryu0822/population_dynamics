diridx = [1,4,9,2,3,5,6,7,8,(10:17)...]
invidx = [1,4,5,2,6,7,8,9,3,(10:17)...]

rslt_gdf = groupby(rslt, :location)
dead_gdf = groupby(dead, :location)

pp04_04_ = []
pp04_08_ = []
for k ∈ 1:17
    rslt_k = rslt_gdf[k]
    dead_k = dead_gdf[k]
    age_c = cumsum(ts_age(rslt_k), dims = 2)
    
    birth = filter(:age => age -> age == 0, rslt_k) |> ts_sum
    death = (dead_k |> ts_sum)
    
    mgrn_k_out = filter(:to => col -> col != 이름_지역[k], filter(:from => col -> col == 이름_지역[k], mgrn))
    rename!(mgrn_k_out, :from => :location)
    select!(mgrn_k_out, Not(:to))
    mgrn_k_in = filter(:to => col -> col == 이름_지역[k], filter(:from => col -> col != 이름_지역[k], mgrn))
    rename!(mgrn_k_in, :to => :location)
    select!(mgrn_k_in, Not(:from))
    age_migration = ts_age(mgrn_k_in) - ts_age(mgrn_k_out)
    
    netmigrat = sum(age_migration, dims = 2)
    netgrowth = birth - death
    
    pp04_01 = plot(xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,14,65,99])
    plot!(pp04_01, 0:99, df_age.y2020[1:100], label = "2020", lw = 2, color = :black)
    plot!(pp04_01, 0:99, df_age.y2045[1:100], label = "2045", lw = 2, color = :navy)
    plot!(pp04_01, 0:99, df_age.y2070[1:100], label = "2070", lw = 2, color = :blue)
    png(pp04_01, "G:/figure/subfigure/04 $k $(name_location[k]) pp04_01.png")
    
    pp04_02 = plot(legend = :topright,
        xlims = (2012, yend), ylims = (0,Inf), xticks = [2012, (2020:10:yend)...],
        ylabel = "Population")
    plot!(pp04_02, 2012:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
    plot!(pp04_02, 2012:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
    plot!(pp04_02, 2012:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
    png(pp04_02, "G:/figure/subfigure/04 $k $(name_location[k]) pp04_02.png")
    
    pp04_04 = plot(
        xlims = (2012, yend), xticks = [2012, (2020:10:yend)...],
        xlabel = "Year", ylabel = "Age-specific Migration")
    plot!(pp04_04, 2012:yend, age_migration[:,3], color =    red, msw = 0, shape = :circ, lw = 2, fa = 0.1, fillrange = 0, label = "65-")
    plot!(pp04_04, 2012:yend, age_migration[:,2], color = orange, msw = 0, shape = :rect, lw = 2, fa = 0.1, fillrange = 0, label = "15-64")
    plot!(pp04_04, 2012:yend, age_migration[:,1], color = yellow, msw = 0, shape = :diamond, lw = 2, fa = 0.1, fillrange = 0, label = "-14")
    png(pp04_04, "G:/figure/subfigure/04 $k $(name_location[k]) pp04_04.png")
    push!(pp04_04_, pp04_04)
    
    pp04_06 = plot(legend = :bottomleft,
        xlims = (2012, yend), xticks = [2012, (2020:10:yend)...],
        ylabel = "Number", bgcolor_legend = colorant"#DDDDDD")
    plot!(pp04_06, 2012:yend, birth, fa = .5, fillrange = 0, color = :green, label = "Birth")
    plot!(pp04_06, 2012:yend, -death, fa = .5, fillrange = 0, color = :black, label = "Death")
    plot!(pp04_06, 2012:yend, birth - death, color = :white, label = "Net growth", lw = 2)
    png(pp04_06, "G:/figure/subfigure/01 0 korea pp01_06.png")

    pp04_08 = plot(ylabel = "Number", xlabel = "Year", xticks = [2012, (2020:10:yend)...], xlims = (2012, yend))
    plot!(pp04_08, 2012:yend, netmigrat, lw = 2, fa = .5, fillrange = 0, color = :olive, label = "net migration")
    plot!(pp04_08, 2012:yend, netgrowth, lw = 2, fa = .5, fillrange = 0, color = :darkgreen, label = "net growth")
    push!(pp04_08_, pp04_08)
    
    plot(pp04_02, pp04_06, pp04_04, pp04_08, 
        plot_title = "$k $(name_location[k])", layout = (2,2), size = (1200, 800),
        leftmargin = 5Plots.mm, dpi = 200)
    png("G:/figure/04 $k $(name_location[k]).png")
end

plot(
    plot(pp04_08_[1],  title = "Seoul"     , legend = :none),
    plot(pp04_08_[2],  title = "Busan"     , legend = :none),
    plot(pp04_08_[4],  title = "Incheon"   , legend = :none),
    plot(pp04_08_[16], title = "Gyoungnam" , legend = :none),
    plot(pp04_08_[9],  title = "Gyounggi"  , legend = :none),
    plot(pp04_08_[10], title = "Gangwon"   , legend = :none),
    plot(pp04_08_[8],  title = "Sejong"    , legend = :none),
    plot(pp04_08_[17], title = "Jeju")     ,
    layout = (4,2), size = (1200, 1200), dpi = 200
)
png("G:/figure/04.png")