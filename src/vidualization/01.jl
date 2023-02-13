begin
    age_l = filter(:age => age ->      age < 15, rslt) |> ts_sum
    age_m = filter(:age => age -> 15 ≤ age < 65, rslt) |> ts_sum
    age_h = filter(:age => age -> 65 ≤ age     , rslt) |> ts_sum
    age_c = cumsum([age_l age_m age_h], dims = 2)
    age_r = 100 * cumsum([(age_l ./ age_m) (age_h ./ age_m)], dims = 2)
    age_cr = 100 * (age_c ./ age_c[:, 3])
    t_super_aging = 2011 + findfirst((age_h ./ age_c[:,3]) .> .21)

    birth = filter(:age => age -> age == 0, rslt) |> ts_sum
    death = dead |> ts_sum
    pop = rslt |> ts_sum
    argmaxpop = argmax(pop)

    df_age = marginal(rslt, :age)

    오십미만 = filter(:age => age -> (age < 50), rslt) |> ts_sum
    오십이상 = filter(:age => age -> (age ≥ 50), rslt) |> ts_sum

    pp01_02 = plot(legend = :topright,
        xlims = (2012, yend), ylims = (0,1.05maximum(pop)), xticks = [2012, (2020:10:yend)...],
        ylabel = "Population")
    plot!(pp01_02, 2012:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
    plot!(pp01_02, 2012:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
    plot!(pp01_02, 2012:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
    scatter!(pp01_02, (2012:yend)[[argmaxpop]], pop[[argmaxpop]], color = :black, shape = :+, label = :none)
    png(pp01_02, "G:/figure/subfigure/01 0 korea pp01_02.png")

    pp01_04 = plot(legend = :topleft,
        xlims = (2012, yend), ylims = (0,150), xticks = [2012, (2020:10:yend)...],
        xlabel = "Year", ylabel = "Ratio(%)")
    plot!(pp01_04, 2012:yend, age_r[:,2], color =    red, fa = .5, fillrange = age_r[:,1], label = "Old-age Dependency Ratio")
    plot!(pp01_04, 2012:yend, age_r[:,1], color = yellow, fa = .5, fillrange = 0         , label = "Youth Dependency Ratio")
    png(pp01_04, "G:/figure/subfigure/01 0 korea pp01_04.png")

    pp01_05 = plot(legend = :none,
        xlims = (2012, yend), ylims = (0,100), xticks = [2012, (2020:10:yend)...],
        xlabel = "Year", ylabel = "Ratio(%)")
    plot!(pp01_05, 2012:yend, age_cr[:,3], color =    red, fa = .5, fillrange = age_cr[:,2], label = "65-")
    plot!(pp01_05, 2012:yend, age_cr[:,2], color = orange, fa = .5, fillrange = age_cr[:,1], label = "15-64")
    plot!(pp01_05, 2012:yend, age_cr[:,1], color = yellow, fa = .5, fillrange = 0          , label = "-14")
    vline!([t_super_aging], style = :dash, color = :white, label = "super aging", lw = 2)
    hline!([79],                           color = :white, label = :none)
    png(pp01_05, "G:/figure/subfigure/01 0 korea pp01_05.png")

    pp01_06 = plot(legend = :bottomleft,
        xlims = (2012, yend), xticks = [2012, (2020:10:yend)...],
        ylabel = "Number", bgcolor_legend = colorant"#DDDDDD")
    plot!(pp01_06, 2012:yend, birth, fa = .5, fillrange = 0, color = :green, label = "Birth")
    plot!(pp01_06, 2012:yend, -death, fa = .5, fillrange = 0, color = :black, label = "Death")
    plot!(pp01_06, 2012:yend, birth - death, color = :white, label = "Natural growth", lw = 2)
    png(pp01_06, "G:/figure/subfigure/01 0 korea pp01_06.png")

    pp01_08 = plot(xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,14,65,99])
    plot!(pp01_08, 0:99, df_age.y2020[1:100], label = "2020", lw = 2, color = :black)
    plot!(pp01_08, 0:99, df_age.y2045[1:100], label = "2045", lw = 2, color = :navy)
    plot!(pp01_08, 0:99, df_age.y2070[1:100], label = "2070", lw = 2, color = :blue)
    png(pp01_08, "G:/figure/subfigure/01 0 korea pp01_08.png")

    pp01_09 = plot(xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,25,50,75,99])
    plot!(pp01_09, 0:99, df_age.y2020[1:100], label = "2020", lw = 2, color = :black)
    plot!(pp01_09, 0:99, df_age.y2045[1:100], label = "2045", lw = 2, color = :navy)
    plot!(pp01_09, 25:99, df_age.y2020[1:75], label = :none, lw = 2, color = :black, ls = :dash, alpha = 0.5)
    png(pp01_09, "G:/figure/subfigure/01 0 korea pp01_09.png")

    plot(pp01_02, pp01_06, pp01_05, pp01_04, pp01_08, pp01_09,
        layout = (3,2), size = 60 .* (16, 16), plot_title = "Korea",
        left_margin = 3Plots.mm, right_margin = 3Plots.mm, dpi = 200)
    png("G:/figure/01 0 korea.png")
    # savefig("G:/figure/eps/01_0_korea.svg")

    수도 = filter(:location => is_capital, rslt)
    서울 = 수도[수도.location .== "서울특별시", :]
    인천 = 수도[수도.location .== "인천광역시", :]
    경기 = 수도[수도.location .== "경기도", :]
    지방 = filter(:location => !is_capital, rslt)
    합계 = cumsum([ts_sum(서울) ts_sum(인천) ts_sum(경기) ts_sum(지방)], dims = 2)
    비율 = 100 * (합계 ./ 합계[:, end])
    pp01_11 = plot(legend = :topright,
        xlims = (2012, yend), ylims = (0,100), xticks = [2012, (2020:10:yend)...], yticks = 0:50:100,
        xlabel = "Year", ylabel = "Ratio(%)", dpi = 200)
    plot!(pp01_11, 2012:yend, 비율[:,end], color = :black, fa = 0.2, fillrange = 비율[:,end - 1], label = "Other")
    plot!(pp01_11, 2012:yend, 비율[:,end - 1], color = :black, fa = 0.4, fillrange = 비율[:,end - 2], label = "Gyounggi")
    plot!(pp01_11, 2012:yend, 비율[:,end - 2], color = :black, fa = 0.6, fillrange = 비율[:,end - 3], label = "Incheon")
    plot!(pp01_11, 2012:yend, 비율[:,1], color = :black, fa = 0.8, fillrange = 0, label = "Seoul")
    png("G:/figure/01 0 korea pp01_11.png")

    pp01_13 = plot(legend = :topright,
    xlims = (2012, yend), ylims = (0,1.05maximum(pop)), xticks = [2012, (2020:10:yend)...],
    ylabel = "Population")
    plot!(pp01_13, 2012:yend, 오십미만 + 오십이상, color = :navy, fa = .5, fillrange = 오십미만, label = "≥ 50")
    plot!(pp01_13, 2012:yend, 오십미만, color = :blue, fa = .5, fillrange = 0, label = "< 50")
    pp01_14 = plot(legend = :topright,
        xlims = (2012, yend), ylims = (0,1), xticks = [2012, (2020:10:yend)...],
        ylabel = "Ratio")
    plot!(pp01_14, 2012:yend, (오십미만 + 오십이상) ./ (오십미만 + 오십이상), color = :navy, fa = .5, fillrange = 오십미만 ./ (오십미만 + 오십이상), label = "≥ 50")
    plot!(pp01_14, 2012:yend, 오십미만 ./ (오십미만 + 오십이상), color = :blue, fa = .5, fillrange = 0, label = "< 50")
    plot(pp01_13, pp01_14, layout = (2,1), dpi = 200)
    png("G:/figure/01 0 korea 50yo.png")
end

# 평균연령l50 = zeros(59)
# 평균연령g50 = zeros(59)
# for gdf = groupby(rslt, :age)
#     a = gdf.age[1]
#     if a < 50
#         평균연령l50 += gdf.age[1] * ts_sum(gdf)
#     else
#         평균연령g50 += gdf.age[1] * ts_sum(gdf)
#     end
# end
# plot(ylims = (0,99), yticks = 0:10:100, ylabel = "Avg. of age", xlabel = "year")
# plot!(2012:2070, 평균연령l50 ./ 오십미만, color = :navy, label = "< 50")
# plot!(2012:2070, 평균연령g50 ./ 오십이상, color = :blue, label = "≥ 50")
# png("50대이상이하평균연령.png")