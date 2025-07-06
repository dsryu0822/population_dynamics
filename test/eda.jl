@time include("datacall.jl")

using StatsPlots

_edata = dropmissing(data)
_edata[:, 3:end] .= 100Matrix(_edata[:, 3:end])
edata = sort(_edata, [:ISO3, :Time])
rename!(edata, :ecnm => :e, :decnm => :de, :yng => :y, :old => :o, :dyng => :dy, :dold => :do)

default(size = [400, 400])
scatter(edata.e, edata.de, color = :black, msw = 0, ms = 1)
scatter(edata.y, edata.dy, color = :black, msw = 0, ms = 1)
scatter(edata.o, edata.do, color = :black, msw = 0, ms = 1)

scatter(edata.y, edata.o, edata.e, color = get(colorschemes[:diverging_linear_bjr_30_55_c53_n256], edata.e ./ 100), msw = 0, ms = 1, label = :none, xlabel = "y", ylabel = "o", zlabel = "e")

@df edata corrplot(cols(3:8), size = [1200, 1200])


gdf = groupby(edata, :ISO3)
for k = 1:length(gdf)
    df = gdf[k]
    p1 = plot(size = [800, 800], xlims = [0, 50], ylims = [0, 30], zlims = [0,100], title = "$(k)_$(df.ISO3[1])", xlabel = "y", ylabel = "o", zlabel = "e")
    if all(isone.(diff(df.Time)))
        plot!(p1, df.y, df.o, df.e, color = :red, label = "data")
        scatter!(p1, [df.y[end]], [df.o[end]], [df.e[end]], color = :red, msw = 0, label = :none)
        plot!(p1, df.y, df.o, zeros(nrow(df)), color = :black, label = "shadow")
        plot!(p1, df.y, fill(30, nrow(df)), df.e, color = :black, label = :none)
        plot!(p1, zeros(nrow(df)), df.o, df.e, color = :black, label = :none)

        p2 = plot(
            plot(df.Time, df.y, color = :black, ylabel = "y", ylims = [0, 50]),
            plot(df.Time, df.o, color = :black, ylabel = "o", ylims = [0, 30]),
            plot(df.Time, df.e, color = :black, ylabel = "e", ylims = [0, 100]),
            layout = (3, 1), size = [800, 800]
        )

        png(
            plot(p1, p2, size = [1600, 800]),
            "G:/population/eda/data3d/$(k)_$(df.ISO3[1]).png"
        )
        # png(p1, "G:/population/eda/data3d/$(k)_$(df.ISO3[1]).png")
    end
end