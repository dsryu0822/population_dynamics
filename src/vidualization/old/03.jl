rslt_gdf = groupby(rslt, :location)
dead_gdf = groupby(dead, :location)
t_super_aging_ = Int64[]
y2021_ = Int64[]
y2070_ = Int64[]
for k ∈ 1:17
    df = rslt_gdf[k]
    age_l = filter(:age => age ->      age < 15, df) |> ts_sum
    age_m = filter(:age => age -> 15 ≤ age < 65, df) |> ts_sum
    age_h = filter(:age => age -> 65 ≤ age     , df) |> ts_sum
    age_c = cumsum([age_l age_m age_h], dims = 2)
    age_r = 100 * cumsum([(age_l ./ age_m) (age_h ./ age_m)], dims = 2)
    age_cr = 100 * (age_c ./ age_c[:, 3])
    t_super_aging = 2011 + findfirst((age_h ./ age_c[:,3]) .> .21)
    push!(t_super_aging_, t_super_aging)
    push!(y2021_, sum(df.y2021))
    push!(y2070_, sum(df.y2070))

    아기 = filter(:age => age -> age == 0, df) |> ts_sum    
    ddf = dead_gdf[k]
    시체 = ddf |> ts_sum
    
    mtrx_mgrn_ = Dict()
    유입 = Int64[]
    유출 = Int64[]
    for t ∈ 2012:yend
        vctr_mgrn = mgrn[:, "y$t"]
        mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
        mtrx_mgrn .*= (1 .- I(17))
        push!(mtrx_mgrn_, t => mtrx_mgrn)
        push!(유입, sum(mtrx_mgrn[k, :]))
        push!(유출, sum(mtrx_mgrn[:, k]))
        if t == 2012
            heatmap(mtrx_mgrn, size = (800,800), ticks = [1,4,9], xlabel = "from", ylabel = "to", title = "Migration Matrix at 2012")
            png("G:/figure/03 00 Migration Matrix")
        end
    end
    
    pp2 = plot(legend = :topright,
    xlims = (2012, yend), ylims = (0,Inf), xticks = [2012, (2020:10:yend)...],
    ylabel = "Population")
    plot!(pp2, 2012:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
    plot!(pp2, 2012:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
    plot!(pp2, 2012:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
    png(pp2, "G:/figure/subfigure/03 $k $(name_location[k]) pp2.png")

    pp4 = plot(legend = :topleft,
    xlims = (2012, yend), ylims = (0,150), xticks = [2012, (2020:10:yend)...],
    ylabel = "Ratio(%)")
    plot!(pp4, 2012:yend, age_r[:,2], color =    red, fa = .5, fillrange = age_r[:,1], label = "Old-age Dependency Ratio")
    plot!(pp4, 2012:yend, age_r[:,1], color = yellow, fa = .5, fillrange = 0         , label = "Youth Dependency Ratio")
    png(pp4, "G:/figure/subfigure/03 $k $(name_location[k]) pp4.png")

    pp5 = plot(legend = :none,
    xlims = (2012, yend), ylims = (0,100), xticks = [2012, (2020:10:yend)...],
    xlabel = "Year", ylabel = "Ratio(%)")
    plot!(pp5, 2012:yend, age_cr[:,3], color =    red, fa = .5, fillrange = age_cr[:,2], label = "65-")
    plot!(pp5, 2012:yend, age_cr[:,2], color = orange, fa = .5, fillrange = age_cr[:,1], label = "15-64")
    plot!(pp5, 2012:yend, age_cr[:,1], color = yellow, fa = .5, fillrange = 0          , label = "-14")
    vline!([t_super_aging], style = :dash, color = :white, label = "super aging", lw = 2)
    png(pp5, "G:/figure/subfigure/03 $k $(name_location[k]) pp5.png")

    pp6 = plot(legend = :bottomleft,
    xlims = (2012, yend), xticks = [2012, (2020:10:yend)...],
    xlabel = "Year", ylabel = "Number")
    plot!(pp6, 2012:yend, filter(:age => age -> age == 0, df) |> ts_sum, fa = .5, fillrange = 0, color = :green, label = "Birth")
    plot!(pp6, 2012:yend, -(ddf |> ts_sum),                               fa = .5, fillrange = 0, color = :black, label = "Death")
    png(pp6, "G:/figure/subfigure/03 $k $(name_location[k]) pp6.png")

    plot(pp2, pp4, pp5, pp6, layout = (2,2), size = 60 .* (16, 9), plot_title = "$k $(name_location[k])",
    left_margin = 3Plots.mm, right_margin = 3Plots.mm, dpi = 200)
    png("G:/figure/02 $k $(name_location[k]).png")
end

poploss = trunc.(1 .- (y2070_ ./ y2021_), digits = 2)
DataFrame(;
name_location,
t_super_aging_,
y2021_,
y2070_,
poploss
)
scatter(t_super_aging_, y2021_, text = name_location,
    legend = :none, 
    ms = 0, xlims = (2018, 2042), xticks = t_super_aging_, yaxis = :log10,)
png("G:/figure/temp1.png")


#  1 & Seoul*      & 2027 &  9,452,735 & 4,397,762 &  53 \% \\
#  2 & Busan       & 2023 &  3,342,975 & 1,341,650 &  59 \% \\
#  3 & Daegu       & 2025 &  2,387,193 &   910,889 &  61 \% \\
#  4 & Incheon*    & 2028 &  2,921,486 & 1,795,707 &  38 \% \\
#  5 & Gwangju     & 2029 &  1,438,662 &   663,556 &  53 \% \\
#  6 & Daejeon     & 2028 &  1,449,428 &   719,868 &  50 \% \\
#  7 & Ulsan       & 2028 &  1,123,933 &   418,347 &  62 \% \\
#  8 & Sejong      & 2040 &    362,853 &   260,817 &  28 \% \\
#  9 & Gyeonggi*   & 2030 & 13,397,444 & 8,825,604 &  34 \% \\
# 10 & Gangwon     & 2021 &  1,530,154 &   946,321 &  38 \% \\
# 11 & Chungbuk    & 2024 &  1,590,962 &   942,687 &  40 \% \\
# 12 & Chungnam    & 2023 &  2,108,672 & 1,314,459 &  37 \% \\
# 13 & Jeonbuk     & 2021 &  1,785,904 &   870,279 &  51 \% \\
# 14 & Jeonnam     & 2017 &  1,832,551 &   922,135 &  49 \% \\
# 15 & Gyeongbuk   & 2020 &  2,619,286 & 1,280,885 &  51 \% \\
# 16 & Gyeongnam   & 2024 &  3,310,829 & 1,452,166 &  56 \% \\
# 17 & Jeju        & 2027 &    670,767 &   410,916 &  38 \% \\