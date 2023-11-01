using CSV, DataFrames, LinearAlgebra, Plots
default(color = :black, legend = :none)

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

@time DATA = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale])
dropmissing!(DATA)

filter!(:Time => x -> x == 2022, DATA)
filter!(:LocTypeName => x -> x == "Country/Area", DATA)
DATA = DATA[:, Not([:Time, :LocTypeName])]

# groupby(DATA, :Location) |> propertynames
# only.(keys(groupby(DATA, :Location).keymap))


areaname = unique(DATA.Location)
ISO3 = unique(DATA.ISO3_code)

gdf = groupby(DATA, :Location)
# vec_pop = [select(df, [:PopMale, :PopFemale]) |> Matrix |> vec for df in gdf]
pop = [[df.PopMale for df in gdf] .+ [df.PopFemale for df in gdf],
       [df.PopMale for df in gdf],
       [df.PopFemale for df in gdf]]
cosM_ = cosine_matrix.(pop)
weights = filter.(!iszero, cosM_)

plot(legend = :best)
histogram!(weights[1], label = "male", alpha = 0.5, color = :red)
histogram!(weights[2], label = "female", alpha = 0.5, color = :blue)
histogram!(weights[3], label = "total", alpha = 0.5, color = :green)

histogram(weights, title = "distribution of cosine similarity", legend = :none)
vline!([0.99], color = :red, lw = 2)
histogram(weights, title = "distribution of cosine similarity", legend = :none, yscale = :log10)
adjM = Symmetric(adjU)

D = 1 ./ (adjM .- minimum(weights))
D[diagind(D)] .= 0
@assert all(D .≥ 0) "distance matrix must be non-negative"

bit_adjM = adjM .> 0.99
histogram(count.(eachrow(bit_adjM)), bin = 10)

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
p3 = plot(size = (16*200, 9*200), framestyle=:none);
for i in 1:n
    for j in fadj[i]
        if i < j
            plot!(p3, [x_[i], x_[j]], [y_[i], y_[j]], color = :black, alpha = 0.1)
        end
    end
end
scatter!(p3, x_, y_, msw = 0, ms = 20, color = assigned, text = areaname);
png("p3.png")