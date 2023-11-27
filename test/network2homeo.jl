using ProgressBars
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering, :Distances]
for package in ProgressBar(packages)
    @eval using $(package)
end
cm = Plots.cm
default(); default(legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))
normalize(v) = (v .- minimum(v)) / (maximum(v) - minimum(v))
getrratio(v) = sum(v[65:end]) / sum(v[[1:15; 65:end]])
getaratio(v) = sum(v[65:end]) / sum(v[1:15])

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

cd(@__DIR__)
DATA_summary = CSV.read("G:/summary.csv", DataFrame)
ISO3 = DATA_summary.ISO3; areaname = DATA_summary.areaname

@time DATA = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale, :PopTotal]);
dropmissing!(DATA)
filter!(:LocTypeName => x -> x == "Country/Area", DATA)

rDATA = filter(:ISO3_code => x -> (x ∈ ISO3), DATA)
ISO3 = unique(rDATA.ISO3_code); areaname = unique(rDATA.Location)

DATAy = collect(groupby(rDATA, :Time));
DATAy = Dict(1950:2021 .=> groupby.(DATAy, :ISO3_code))

vec_Pop_ = getproperty.(collect(DATAy[1950]), :PopTotal);
Op = getrratio.(vec_Pop_) |> normalize

cosM = vec_Pop_ |> cosine_matrix |> Symmetric

# using Graphs, GraphRecipes, Plots, LinearAlgebra, NetworkLayout
# n = 100
# adjM = Symmetric(rand(n, n) .< 0.02)
# adjM[diagind(adjM)] .= 0
# G = Graph(adjM)

# CG = []
# _G = deepcopy(G)
# deg = degree(_G)
# while any(deg .> 2)
#     k = argmax(deg)
#     deg[k] = -1

#     selected = setdiff(_G.fadjlist[k], findall(deg .≤ 0))
#     deg[selected] .= -2
#     for i in selected
#         for j in _G.fadjlist[i]
#             rem_edge!(_G, i, j)
#             if j ≠ k
#                 add_edge!(_G, j, k)
#             end
#         end
#     end
#     [rem_edge!(_G, i, k) for i in _G.fadjlist[k]]
#     push!(CG, k => selected)
# end
# CG
# _G.fadjlist[21]
# _G.fadjlist[12]
# deg[12]
# G.fadjlist[12]


# deg[first.(CG)] .= length.(last.(CG))

# x_, y_ = eachrow(stack(spring(G)))
# # x_[(deg .== -2)] .= 20; y_[(deg .== -2)] .= 20
# x_[(deg .==  0)] .= -20; y_[(deg .==  0)] .= -20
# graphplot(_G, x = x_, y = y_, names = deg, nodesize = 0.1, nodeweights = deg,
# axis_buffer = 0, lims = [-12, 12], size = [800, 800])

# graphplot(G, x = x_, y = y_, names = deg, nodesize = 0.1, nodeweights = deg,
# axis_buffer = 0, lims = [-12, 12], size = [800, 800])


# H_ = []
# _G = deepcopy(G)
# deg = degree(_G)
# while any(deg .> 0)
#     k = argmax(deg)
#     deg[k] = -1
#     push!(H_, k => setdiff(_G.fadjlist[k], findall(deg .< 0)))
#     deg[_G.fadjlist[k]] .= -2
# end
# @assert all(reduce(vcat, last.(H_)) .== union(last.(H_)...))
# @assert 100 == (length(H_) + length(union(last.(H_)...)) + count(deg .≥ 0))

# _adjM = zeros(Int64, n,n)

# ball_ = last.(H_)

# # newnodes = [first.(H_); findall(deg .> 0)]
# _G.fadjlist[findall(deg .== 2)]
# _G.fadjlist[findall(deg .== 1)]
# _G.fadjlist[findall(deg .== 0)]



H_ = []
_G = deepcopy(G)
deg = degree(_G)
while any(deg .> 0)
    k = argmax(deg)
    deg[k] = -1

    push!(H_, k => k)
    push!(H_, (setdiff(_G.fadjlist[k], findall(deg .< 0)) .=> k)...)
    deg[_G.fadjlist[k]] .= -2
end
isolated = findall(iszero.(deg))
push!(H_, (isolated .=> isolated)...)
H_ = Dict(H_)

_adjM = zeros(Int64, n,n)
for i in 1:n, j in 1:n
    _adjM[H_[i],H_[j]] += adjM[i,j]
end
_adjM[diagind(_adjM)] .= 0
_adjM = Symmetric(_adjM)
findall(deg .== -1)

represent = sort(unique(values(H_)))
nw = [count(values(H_) .== k) for k in 1:n]

_layout = spring(adjM)
__layout = stack(_layout)

__layout[:, deg .== -2] .= -20
__layout[:, isolated] = randn(2, length(isolated)) .+ [5,5]

x_, y_ = eachrow(__layout)
graphplot(_adjM, x = x_, y = y_, lims = [-8, 8], names = nw, size = [400,400],
nodesize = 3, nodeweights = nw, ma = 0, shape = :rect)
# graphplot( adjM, x = x_, y = y_, lims = [-8, 8], names = 1:n, size = [400,400])

