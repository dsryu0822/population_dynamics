using CairoMakie, LaTeXStrings
using StatsBase

symtemp = ' '

seola = Figure()
time_series1 = Axis(seola[1,1], title = "People", ylabel = L"Population Size $S$", yscale = log10)
ranking1 = Axis(seola[2,1], yreversed = true, xlabel = "Year", ylabel = L"Ranking $R$", yticks = [1,5,10,17])
temp1 = zeros(Int64, 59, 17)
for (idf, df) ∈ enumerate(groupby(rslt, :location))
    temp1[:, idf] = ts_sum(df)
    lines!(time_series1, 2012:2070, temp1[:, idf])
end
current_figure()

ts_ranking = hcat(ordinalrank.(eachrow(temp1), rev = true)...)
for (i, mat) in enumerate(eachrow(ts_ranking))
    if i < 8
        symtemp = :vline
    elseif i == 8
        symtemp = :circle
    else
        symtemp = ' '
    end
    scatterlines!(ranking1, 2012:2070, mat, linewidth = 2, label = 이름_지역[i], marker = symtemp)
end
seola[1:2,2] = axislegend(ranking1)
current_figure()

vec(sum(.!iszero.(diff(ts_ranking, dims = 2)), dims = 1))

save("G:/figure/ranking_population.svg", seola)


area = [605243961,770170750,883698218,1066465186,501113181,539503587,1062327676,10196731502,464918346,16829671465,7406988694,8246960405,8072147335,12358940174,19034795138,10541895013,1850278723]

exy = Figure()
time_series2 = Axis(exy[1,1], title = L"People / km$^2$", ylabel = L"Population Size $S$", yscale = log10)
ranking2 = Axis(exy[2,1], yreversed = true, xlabel = "Year", ylabel = L"Ranking $R$", yticks = [1,5,10,17])
temp2 = zeros(Float64, 59, 17)
for (idf, df) ∈ enumerate(groupby(rslt, :location))
    temp2[:, idf] = ts_sum(df) ./ area[idf]
    lines!(time_series2, 2012:2070, temp2[:, idf])
end
current_figure()

ts_ranking = hcat(ordinalrank.(eachrow(temp2), rev = true)...)
for (i, mat) in enumerate(eachrow(ts_ranking))
    if i < 8
        symtemp = :vline
    elseif i == 8
        symtemp = :circle
    else
        symtemp = ' '
    end
    scatterlines!(ranking2, 2012:2070, mat, linewidth = 2, label = 이름_지역[i], marker = symtemp)
end
exy[1:2,2] = axislegend(ranking2)
current_figure()

save("G:/figure/ranking_density.svg", exy)

normalized_df = hcat([(df |> ts_sum) ./ (rslt |> ts_sum) for df in groupby(rslt, :location)]...)
using StatsBase
yeonjung = scatterlines(2012:2070, entropy.(eachrow(normalized_df), 17), axis = (; title = "entropy"))
save("G:/figure/entropy.svg", yeonjung)


linr_netmigrat = hcat(netmigrat_...)
linr_netgrowth = hcat(netgrowth_...)
rank_netmigrat = hcat(ordinalrank.(eachrow(linr_netmigrat), rev = true)...)'
rank_netgrowth = hcat(ordinalrank.(eachrow(linr_netgrowth), rev = true)...)'

bona =  Figure()
time_series3 = Axis(bona[1,1], title = "Net migration", ylabel = L"Number")
ranking3 = Axis(bona[2,1], yreversed = true, xlabel = "Year", ylabel = L"Ranking $R$", yticks = [1,5,10,17])
subin = Figure()
time_series4 = Axis(subin[1,1], title = "Net growth", ylabel = L"Number")
ranking4 = Axis(subin[2,1], yreversed = true, xlabel = "Year", ylabel = L"Ranking $R$", yticks = [1,5,10,17])
for j in 1:17
    if j < 8
        symtemp = :vline
    elseif j == 8
        symtemp = :circle
    else
        symtemp = ' '
    end
    lines!(time_series3, 2012:2070, linr_netmigrat[:, j])
    lines!(time_series4, 2012:2070, linr_netgrowth[:, j])
    scatterlines!(ranking3, 2012:2070, rank_netmigrat[:, j], linewidth = 2, label = 이름_지역[j], marker = symtemp)
    scatterlines!(ranking4, 2012:2070, rank_netgrowth[:, j], linewidth = 2, label = 이름_지역[j], marker = symtemp)
end
bona[1:2,2] = axislegend(ranking3)
subin[1:2,2] = axislegend(ranking4)
save("G:/figure/ranking_netmigrat.svg", bona)
save("G:/figure/ranking_netgrowth.svg", subin)
