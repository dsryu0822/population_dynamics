using ProgressBars
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering]
for package in ProgressBar(packages)
    @eval using $(package)
end
cm = Plots.cm
default(); default(legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))

function cosine_matrix(v)
    n = length(v)
    cosM = zeros(n, n)
    for i in 1:n, j in 1:n
        if i < j
            cosM[i,j] = cosine_similarity(v[i], v[j])
        end
    end
    return cosM
end

cd("//155.230.155.221/ty/Population")

DATA_summary = CSV.read("summary.csv", DataFrame)
@time DATA = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale, :PopTotal]);
dropmissing!(DATA)

filter!(:LocTypeName => x -> x == "Country/Area", DATA)
data = filter(:Time => x -> x == 2021, DATA)
data = data[:, Not([:Time, :LocTypeName])]
rename!(data, :ISO3_code => :ISO3)

# ISO3 = unique(data.ISO3)
gdf = groupby(data, :ISO3) # 236 countries
popcut = combine(gdf, :PopTotal => sum => :Pop)
data = reduce(vcat, gdf[popcut.Pop .> 1000])
pop = popcut.Pop[popcut.Pop .> 1000]

areaname = unique(data.Location)
ISO3 = unique(data.ISO3)
gdf = groupby(data, :ISO3) # 161 countries
CSV.write("ISO3.csv", sort(DataFrame(; ISO3, areaname, pop), :ISO3), bom = true)

n = length(gdf)
pop = [[df.PopTotal for df in gdf],
       [df.PopMale for df in gdf],
       [df.PopFemale for df in gdf]]
cosM = cosine_matrix.(pop)[1]
weightss = filter(!iszero, cosM)
cosM = Symmetric(cosM)

idx_top10_ = getindex.(reverse.(sortperm.(eachrow(cosM))), Ref(1:10))
wgt_top10_ = getindex.(eachrow(cosM), idx_top10_)
# CSV.write("CStop10.csv", DataFrame(stack(getindex.(Ref(areaname), idx_top10_)), areaname), bom = true)
# CSV.write("CStop10w.csv", DataFrame(stack(wgt_top10_), areaname), bom = true)

# for i = 1:n
#     p1 = [plot(gdf[i].PopTotal, color = :black, title = areaname[i])]
#     for j in 1:10
#         k = idx_top10_[i][j]
#         push!(p1,
#             plot(gdf[k].PopTotal, color = :black, title = areaname[k], xlabel = round(wgt_top10_[i][j], digits = 4))
#         )
#     end
#     p1_ = plot(p1..., size = 80*[16, 9])
#     png(join([lpad(i, 3, '0'), ISO3[i], areaname[i]], "_"))
# end

# histogram(weightss, bins = 0.5:0.001:1, lw = 0, color = :black, title = "Distribution of cosine similarity")
# heatmap(Matrix(cosM + 0.78*I), color = :RdBu, size = [700, 600], legend = :best)

# θ_ = 0.8:0.001:0.999
# degree_ = [count(cosM .> θ) for θ = θ_]
# plot(θ_, degree_, xlabel = "Threshold θ", ylabel = "sum of degree", yscale = :log10, color = :black, ylims = [50,50000])

# plot(θ_[1:(end-1)], diff(log10.(degree_)))
# plot(θ_[2:(end-1)], diff(diff(log10.(degree_))))

# θ_ = [0.8, 0.9, 0.95, 0.99, 0.999]
# degree_ = [sum.(eachrow(cosM .> θ)) for θ = θ_]
# for k in 1:5
#     histogram(degree_[k], color = k, bins = 1:maximum(degree_[k]), xlabel = :dgree, title = "θ = $(θ_[k])", size = [400, 400], lw = 0, xlims = (0, 170), ylims = (0, 13))
#     png("$k")
# end

oldrate = combine(gdf, :PopTotal => (x -> sum(x[65:end]) / (sum(x[1:15]) + sum(x[65:end]))) => :ratio).ratio
oldrate .-= minimum(oldrate)
oldrate ./= maximum(oldrate)
color1 = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), oldrate)
colorV = color1

# for θ = 0.98:0.001:0.999
for θ = [0.995]
    Random.seed!(0)
    graphplot(cosM .> θ, method = :stress, fontsize = 10, nodecolor = colorV, msw = 0, size = [3000, 3000], nodeshape = :rect, la = 0.1, names = ISO3, curves = false, dpi = 300); png("network cut $(rpad(θ, 5, '0'))")
end

# graphplot(cosM .> θ, method = :stress, fontsize = 10, nodecolor = oldrate, msw = 0, size = [3000, 3000], nodeshape = :rect, la = 0.1, names = ISO3, curves = false, dpi = 300); png("temp")

screeplot = plot(1:10, [kmeans(cosM, k).totalcost for k in 1:10], xticks = 1:10, color = :black, lw = 2, title = "Scree plot", xlabel = L"K", ylabel = "Total cost")

kmeansed = kmeans(cosM, 3)
ij = mod.(kmeansed.assignments, 3) |> sortperm
heatmap(Matrix(cosM[ij, ij] + 0.78*I), color = :RdBu, size = [700, 600], legend = :best)
color2 = kmeansed.assignments
colorV = color2

scatter(oldrate, color2 + 0.1randn(n), yticks = 1:3, size = [400, 400], xlabel = L"O_{p}", ylabel = "Cluster", color = color2, msw = 0)


areaname[findall((color2 .== 2) .&& (oldrate .< 0.5))]
println()
# # hclusted = hclust(D)
# # assigned = cutree(hclusted, k = 6)
# kmedoidsed = kmedoids(D, 6)
# assigned = kmedoidsed.assignments
# kmedoidsed.counts

# areaname[assigned .== 1]
# areaname[assigned .== 2]
# areaname[assigned .== 3]
# areaname[assigned .== 4]
# areaname[assigned .== 5]
# areaname[assigned .== 6]
# areaname[assigned .== assigned[findfirst(ISO3 .== "KOR")]]

# bit_adjM = adjM .> 0.99
# fadj = findall.(eachrow(bit_adjM))

# x_ = 100rand(n) .+ 100mod.(assigned, 3)
# y_ = 100rand(n) .+ 100mod.(assigned, 2)
# p3 = plot(size = (16*200, 9*200), framestyle=:none);
# for i in 1:n
#     for j in fadj[i]
#         if i < j
#             plot!(p3, [x_[i], x_[j]], [y_[i], y_[j]], color = :black, alpha = 0.1)
#         end
#     end
# end
# scatter!(p3, x_, y_, msw = 0, ms = 20, color = assigned, text = areaname);
# png("p3.png")