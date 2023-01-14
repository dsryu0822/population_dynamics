include("basic.jl")

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
    t_super_aging = 2020 + findfirst((age_h ./ age_c[:,3]) .> .21)
    push!(t_super_aging_, t_super_aging)
    push!(y2021_, sum(df.y2021))
    push!(y2070_, sum(df.y2070))

    아기 = filter(:age => age -> age == 0, df) |> ts_sum    
    ddf = dead_gdf[k]
    시체 = ddf |> ts_sum
    
    mtrx_mgrn_ = Dict()
    유입 = Int64[]
    유출 = Int64[]
    for t ∈ 2021:yend
        vctr_mgrn = mgrn[:, "y$t"]
        mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
        mtrx_mgrn .*= (1 .- I(17))
        push!(mtrx_mgrn_, t => mtrx_mgrn)
        push!(유입, sum(mtrx_mgrn[k, :]))
        push!(유출, sum(mtrx_mgrn[:, k]))
        if t == 2021
            heatmap(mtrx_mgrn, size = (800,800), ticks = [1,4,9], xlabel = "from", ylabel = "to", title = "Migration Matrix at 2021")
            png("G:/figure/03 00 Migration Matrix")
        end
    end
    
    pp2 = plot(legend = :topright,
    xlims = (2021, yend), ylims = (0,Inf), xticks = [2021, (2030:10:yend)...],
    ylabel = "Population")
    plot!(pp2, 2021:yend, age_c[:,3], color =    red, fa = .5, fillrange = age_c[:,2], label = "65-")
    plot!(pp2, 2021:yend, age_c[:,2], color = orange, fa = .5, fillrange = age_c[:,1], label = "15-64")
    plot!(pp2, 2021:yend, age_c[:,1], color = yellow, fa = .5, fillrange = 0         , label = "-14")
    # annotate!(pp2, [2022], [52], [(L"\times 10^6", 8)])
    png(pp2, "G:/figure/subfigure/01 0 korea pp2.png")

    pp4 = plot(legend = :topleft,
    xlims = (2021, yend), ylims = (0,150), xticks = [2021, (2030:10:yend)...],
    ylabel = "Ratio(%)")
    plot!(pp4, 2021:yend, age_r[:,2], color =    red, fa = .5, fillrange = age_r[:,1], label = "Elder Dependency Ratio")
    plot!(pp4, 2021:yend, age_r[:,1], color = yellow, fa = .5, fillrange = 0         , label = "Child Dependency Ratio")
    png(pp4, "G:/figure/subfigure/01 0 korea pp4.png")

    pp5 = plot(legend = :none,
    xlims = (2021, yend), ylims = (0,100), xticks = [2021, (2030:10:yend)...],
    xlabel = "Year", ylabel = "Ratio(%)")
    plot!(pp5, 2021:yend, age_cr[:,3], color =    red, fa = .5, fillrange = age_cr[:,2], label = "65-")
    plot!(pp5, 2021:yend, age_cr[:,2], color = orange, fa = .5, fillrange = age_cr[:,1], label = "15-64")
    plot!(pp5, 2021:yend, age_cr[:,1], color = yellow, fa = .5, fillrange = 0          , label = "-14")
    vline!([t_super_aging], style = :dash, color = :white, label = "super aging", lw = 2)
    png(pp5, "G:/figure/subfigure/01 0 korea pp5.png")

    pp6 = plot(legend = :bottomleft,
    xlims = (2021, yend), xticks = [2021, (2030:10:yend)...],
    xlabel = "Year", ylabel = "Number")
    plot!(pp6, 2021:yend, filter(:age => age -> age == 0, df) |> ts_sum, fa = .5, fillrange = 0, color = :green, label = "Birth")
    plot!(pp6, 2021:yend, -(ddf |> ts_sum),                               fa = .5, fillrange = 0, color = :black, label = "Death")
    png(pp6, "G:/figure/subfigure/01 0 korea pp6.png")

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