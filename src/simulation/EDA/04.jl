include("basic.jl")

diridx = [1,4,9,2,3,5,6,7,8,(10:17)...]
invidx = [1,4,5,2,6,7,8,9,3,(10:17)...]

rslt_gdf = groupby(rslt, :location)
dead_gdf = groupby(dead, :location)

pp04_5_ = []
for k ∈ 1:17
    df = rslt_gdf[k]

    아기 = filter(:age => age -> age == 0, df) |> ts_sum    
    ddf = dead_gdf[k]
    시체 = ddf |> ts_sum
    
    mtrx_mgrn_ = Dict()
    유입시도 = Int64[]
    유출시도 = Int64[]
    유입수도지방2 = Int64[]
    유출수도지방2 = Int64[]
    유입수도지방17 = Int64[]
    유출수도지방17 = Int64[]
    for t ∈ 2021:yend
        vctr_mgrn = mgrn[:, "y$t"]
        mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
        mtrx_mgrn .*= (1 .- I(17))
        push!(유입시도, sum(mtrx_mgrn[k, :]))
        push!(유출시도, sum(mtrx_mgrn[:, k]))
        
        temp = mtrx_mgrn[diridx, diridx]
        temp[1:3,1:3] .= 0
        temp[4:17,4:17] .= 0
        push!(유입수도지방2, sum(temp[1:3, :]))
        push!(유출수도지방2, sum(temp[:, 1:3]))

        mtrx_mgrn = temp[invidx, invidx]
        push!(mtrx_mgrn_, t => mtrx_mgrn)
        push!(유입수도지방17, sum(mtrx_mgrn[k, :]))
        push!(유출수도지방17, sum(mtrx_mgrn[:, k]))
       
        # vctr_mgrn_age = vec(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 2))
        # mtrx_mgrn_child = reshape(sum(reshape(vctr_mgrn_age, 17, 17, 17)[1:10,:,:], dims = 1), 17, 17)
        # mtrx_mgrn_elder = reshape(sum(reshape(vctr_mgrn_age, 17, 17, 17)[11:end,:,:], dims = 1), 17, 17)
    end

    pp1 = plot(title = "$k $(name_location[k])", xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,25,50,75,99])
    plot!(pp1, 0:99, df.y2021[1:100] + df.y2021[101:200], label = "2021", lw = 2, color = :black)
    plot!(pp1, 0:99, df.y2046[1:100] + df.y2046[101:200], label = "2046", lw = 2, color = :navy)
    # plot!(pp1, 0:99, df.y2070[1:100] + df.y2070[101:200], label = "2070", lw = 2, color = :blue)
    plot!(pp1, 25:99, df.y2021[1:75] + df.y2021[101:175], label = :none, lw = 2, color = :black, ls = :dash, alpha = 0.5)
    # plot!(pp1, 25:99, df.y2046[1:75] + df.y2046[101:175], label = :none, lw = 2, color = :navy, ls = :dash, alpha = 0.5)

    pp2 = plot(title = "$k $(name_location[k])", xlabel = "age", ylabel = "Population", xlims = (0,99), ylims = (0, Inf), xticks = [0,25,50,75,99])
    plot!(pp2, 0:99, df.y2021[1:100] + df.y2021[101:200], label = "2021", lw = 2, color = :black)
    plot!(pp2, 0:99, df.y2046[1:100] + df.y2046[101:200], label = "2046", lw = 2, color = :navy)
    plot!(pp2, 0:99, df.y2070[1:100] + df.y2070[101:200], label = "2070", lw = 2, color = :blue)
    # plot!(pp2, 25:99, df.y2021[1:75] + df.y2021[101:175], label = :none, lw = 2, color = :black, ls = :dash, alpha = 0.5)
    # plot!(pp2, 25:99, df.y2046[1:75] + df.y2046[101:175], label = :none, lw = 2, color = :navy, ls = :dash, alpha = 0.5)
    
    pp3 = plot(ylabel = "Number", xlabel = "Year", xticks = [2021, (2030:10:yend)...], xlims = (2021, yend), legend = :none)
    plot!(pp3, 2021:yend, 유입시도 - 유출시도, lw = 2, fa = .5, fillrange = 0, color = :olive, label = "net migration")
    plot!(pp3, 2021:yend, 아기  -  시체, lw = 2, fa = .5, fillrange = 0, color = :darkgreen, label = "net growth")
    
    pp4 = plot(ylabel = "Number(bipartite)", xlabel = "Year", xticks = [2021, (2030:10:yend)...], xlims = (2021, yend))
    plot!(pp4, 2021:yend, 유입수도지방17 - 유출수도지방17, lw = 2, fa = .5, fillrange = 0, color = :olive, label = "net migration")
    plot!(pp4, 2021:yend, 아기  -  시체, lw = 2, fa = .5, fillrange = 0, color = :darkgreen, label = "net growth")
    
    pp04_5 = plot(ylabel = "Number(bipartite)", xlabel = "Year", xticks = [2021, (2030:10:yend)...], xlims = (2021, yend))
    plot!(pp04_5, 2021:yend, 유입수도지방2 - 유출수도지방2, lw = 2, fa = .5, fillrange = 0, color = :olive, label = "net migration")
    plot!(pp04_5, 2021:yend, 아기  -  시체, lw = 2, fa = .5, fillrange = 0, color = :darkgreen, label = "net growth")
    push!(pp04_5_, pp04_5)

    plot(pp1, pp2, pp3, pp4, layout = (2,2), size = 60 .* (16, 9), leftmargin = 6Plots.mm, rightmargin = 6Plots.mm, dpi = 200)
    png("G:/figure/04 $k $(name_location[k]).png")
end
pp04_5_[12]