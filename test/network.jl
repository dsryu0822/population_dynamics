using CSV, DataFrames, LinearAlgebra, Plots
default(color = :black, legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))
@time DATA = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale])
dropmissing!(DATA)

filter!(:Time => x -> x == 2022, DATA)
filter!(:LocTypeName => x -> x == "Country/Area", DATA)
DATA = DATA[:, Not(:Time)]
DATA = DATA[:, Not(:LocTypeName)]

# groupby(DATA, :Location) |> propertynames
# only.(keys(groupby(DATA, :Location).keymap))


areaname = unique(DATA.Location)
ISO3 = unique(DATA.ISO3_code)

gdf = groupby(DATA, :Location)
vec_pop = [select(df, [:PopMale, :PopFemale]) |> Matrix |> vec for df in gdf]

n = length(vec_pop)
adjU = zeros(n, n)
for i in 1:n, j in 1:n
    if i < j
        adjU[i,j] = cosine_similarity(vec_pop[i], vec_pop[j])
    end
end
weights = filter(!iszero, adjU)
histogram(weights, title = "distribution of cosine similarity", legend = :none)
vline!([0.99], color = :red, lw = 2)
histogram(weights, title = "distribution of cosine similarity", legend = :none, yscale = :log10)
adjM = Symmetric(adjU)

D = 1 ./ (adjM .- minimum(weights))
D[diagind(D)] .= 0
@assert all(D .≥ 0) "distance matrix must be non-negative"

using Clustering, StatsPlots
# hclusted = hclust(D)
# assigned = cutree(hclusted, k = 6)
kmedoidsed = kmedoids(D, 6)
assigned = kmedoidsed.assignments
kmedoidsed.counts

areaname[assigned .== 1]
areaname[assigned .== 2]
areaname[assigned .== 3]
areaname[assigned .== 4]
areaname[assigned .== 5]
areaname[assigned .== 6]
areaname[assigned .== assigned[findfirst(ISO3 .== "KOR")]]

bit_adjM = adjM .> 0.99
fadj = findall.(eachrow(bit_adjM))

x_ = 100rand(n) .+ 100mod.(assigned, 3)
y_ = 100rand(n) .+ 100mod.(assigned, 2)
p1 = plot(size = (16*200, 9*200), framestyle=:none);
for i in 1:n
    for j in fadj[i]
        plot!(p1, [x_[i], x_[j]], [y_[i], y_[j]], color = :black, alpha = 0.1)
    end
end
scatter!(p1, x_, y_, msw = 0, ms = 20, color = assigned, text = areaname);
png("p1.png")