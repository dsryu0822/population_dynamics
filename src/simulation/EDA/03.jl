include("basic.jl")

rslt = rslt_[1]
dead = dead_[1]
mgrn = mgrn_[1]

pp0 = plot(legend = :bottomright, xlabel = "year", ylabel = "population")
df = deepcopy(rslt)
age_l = filter(:age => age ->      age < 15, df) |> ts_sum
age_m = filter(:age => age -> 15 ≤ age < 65, df) |> ts_sum
age_h = filter(:age => age -> 65 ≤ age     , df) |> ts_sum
age_c = cumsum([age_l age_m age_h], dims = 2)
t_super_aging = 2020 + findfirst((age_h ./ age_c[:,3]) .> .21)
plot!(pp0, 2021:yend, age_c[:,3], color =    :red, fa = .5, fill = age_c[:,2], label = "65~")
plot!(pp0, 2021:yend, age_c[:,2], color = :orange, fa = .5, fill = age_c[:,1], label = "15~64")
plot!(pp0, 2021:yend, age_c[:,1], color = :yellow, fa = .5, fill = 0         , label = "~14")
vline!([t_super_aging], style = :dash, color = :black, label = "super aging")
png(pp0, "D:/figure/03 00 korea.png")

rslt_gdf = groupby(rslt, :location)
dead_gdf = groupby(dead, :location)
pyramid_ = Plots.Plot[]
t_super_aging_ = Int64[]
for k ∈ 1:17
    df = rslt_gdf[k]
    age_l = filter(:age => age ->      age < 15, df) |> ts_sum
    age_m = filter(:age => age -> 15 ≤ age < 65, df) |> ts_sum
    age_h = filter(:age => age -> 65 ≤ age     , df) |> ts_sum
    age_c = cumsum([age_l age_m age_h], dims = 2)
    age_r = cumsum([(age_l ./ age_m) (age_h ./ age_m)], dims = 2)
    t_super_aging = 2020 + findfirst((age_h ./ age_c[:,3]) .> .21)
    push!(t_super_aging_, t_super_aging)

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
            png("D:/figure/03 00 Migration Matrix")
        end
    end
    
    pp1 = plot(xticks = 0:30:99, title = "$k $(name_location[k])", xlabel = "age", ylabel = "population", legend = :outertopright)
    plot!(pp1, 0:99, combine(groupby(df, :age), :y2021 => sum => :pyramid).pyramid, label = "2021", lw = 2, color = :black)
    plot!(pp1, 0:99, combine(groupby(df, :age), :y2050 => sum => :pyramid).pyramid, label = "2050", lw = 2, color =   :red)
    push!(pyramid_, pp1)
    png(pp1, "D:/figure/03 00 pp1-$k")

    pp2 = plot(legend = :outertopright, xlabel = "year", ylabel = "population")
    plot!(pp2, 2021:yend, age_c[:,3], color =    :red, fa = .5, fill = age_c[:,2], label = "65~")
    plot!(pp2, 2021:yend, age_c[:,2], color = :orange, fa = .5, fill = age_c[:,1], label = "15~64")
    plot!(pp2, 2021:yend, age_c[:,1], color = :yellow, fa = .5, fill = 0         , label = "~14")
    vline!([t_super_aging], style = :dash, color = :black, label = "super aging")
    png(pp2, "D:/figure/03 00 pp2-$k")

    pp3 = plot(legend = :outertopright)
    plot!(pp3, 2021:yend, 아기 - 시체, lw = 2, fa = .5, fill = 0, label = "net growth")
    plot!(pp3, 2021:yend, 유입 - 유출, lw = 2, fa = .5, fill = 0, label = "net migration")
    png(pp3, "D:/figure/03 00 pp3-$k")

    pp4 = plot(legend = :outertopright, xlabel = "year", ylabel = "population")
    plot!(pp4, 2021:yend, age_r[:,2], color =    :red, fa = .5, fillrange = age_r[:,1], label = "older ratio")
    plot!(pp4, 2021:yend, age_r[:,1], color = :yellow, fa = .5, fillrange = 0         , label = "young ratio")
    png(pp4, "D:/figure/03 00 pp4-$k")
    
    plot(pp1, pp2, pp3, pp4, layout = (2,2), size = (1200, 900))
    png("D:/figure/03 00 $k $(name_location[k]).png")
end