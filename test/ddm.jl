using ProgressMeter
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :SparseArrays, :Colors, :ColorSchemes, :JLD2, 
            :Clustering, :Symbolics]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm
default(); default(legend = :none)

grid(tick::Float64, x_ = [0, 1], y_ = [0, 1]) = Base.product(x_[1]:tick:x_[2], y_[1]:tick:y_[2]) |> collect |> vec .|> collect
grid(tick::Int64, x_ = [0, 1], y_ = [0, 1]) = Base.product(LinRange(x_[1],x_[2],tick), LinRange(y_[1],y_[2],tick)) |> collect |> vec .|> collect
square(X) = reshape(X, isqrt(length(X)), isqrt(length(X)))
function solve(f, x; T = 10, h = 1e-2, shorten = true)
    x_ = [x]
    for t in 1:h:T
        if !(-1000 ≤ x_[end][1] ≤ 1000) || !(-1000 ≤ x_[end][2] ≤ 1000)
            break
        elseif rand() < 0.1
            if norm(f(x_[end])) < 1e-6
                break
            end
        end
        push!(x_, x_[end] + h*f(x_[end]))
    end    
    return shorten ? x_[1:max(1, (length(x_) ÷ 100)):end] : x_
end

# @time include("datacall.jl")
data = CSV.read("cached.csv", DataFrame)
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

const xmax = 0.7; const ymax = 0.7; resolution = 100
qargs = (; xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [800, 800], framestyle = :box, formatter = x -> "$(round(100x; digits = 2))%")
layer_ = [plot(; qargs...)]

### SINDy
# data[data.ISO3 .== "KOR", :]
odata = dropmissing(select(data, Not(:ecnm)))
sort(sort(odata, :Time), :ISO3)
f = SINDy(odata, [:dyng, :dold], [:yng, :old], N = 6)
print(f, ["y", "o"])

### Find fixed points
pos = grid(resolution, [0, xmax], [0, ymax])
dir = f.(pos)
fixedcandy = pos[norm.(dir) .< 1e-4]
dbscaned = dbscan(stack(fixedcandy), 0.01)
fixedcandy = fixedcandy[rand.(getproperty.(dbscaned.clusters, :core_indices))]
fixedpoint = [findfixed(f, fc; maxitr = 50000, atol = 1e-8) for fc in fixedcandy]
idx_unfixed = findall(norm.(f.(fixedpoint)) .> 1e-8)
unfixed = fixedpoint[idx_unfixed]
gridsearch = [pairwiseadd(uf, grid(1000, [-0.01, 0.01], [-0.01, 0.01])) for uf in unfixed]
for (i, j) in enumerate(idx_unfixed)
    better = norm.(f.(gridsearch[i]))
    better[argmin(better)]
    fixedpoint[j] .= gridsearch[i][argmin(better)]
end
fixedpoint .=> norm.(f.(fixedpoint))
fp5 = [0.19891399160390585, 0.22150148422443738]

### Stability analysis
@variables t x y
∂x = Differential(x)
∂y = Differential(y)
Θxy = Θ([x, y], N = 6)
∂xΘ = ∂x.(Θxy) .|> expand_derivatives
∂yΘ = ∂y.(Θxy) .|> expand_derivatives
eigenvalues = []
for j ∈ eachindex(fixedpoint)
    fp = (Dict(x => fixedpoint[j][1], y => fixedpoint[j][2]),)
    J = Float64.([substitute.(∂xΘ, fp) .|> Symbolics.value
                substitute.(∂yΘ, fp) .|> Symbolics.value]) * Matrix(f.matrix)
    push!(eigenvalues, eigen(J).values)
end
encode_stability = sum.(real.(eigenvalues) .> 0)
stabiliary_colors = [:black, :gray, :red][encode_stability .+ 1]
push!(layer_, deepcopy(fisrt(layer_)))
scatter!(layer_[2], first.(fixedpoint), last.(fixedpoint), mc = stabiliary_colors, ms = 10, msw = 0, txt = "\t" .* string.(1:7))

### Separatrix
push!(layer_, deepcopy(fisrt(layer_)))
@showprogress for ic = pos
    sol = solve(f, ic, T = 10000, shorten = false)
    endpoint = last(sol); push!(Z_, ic => endpoint)
    plot!(layer_[3], first.(sol), last.(sol), lw = 1, color = :black, alpha = 0.1)
end
basin = DataFrame(x0 = first.(first.(Z_)), y0 = last.(first.(Z_)), xend = first.(last.(Z_)), yend = last.(last.(Z_)))
CSV.write("basin.csv", basin); basin = CSV.read("basin.csv", DataFrame)

### Basin of attraction
temp = square((0.19.< basin.xend .< 0.222) .&& (0.19 .< basin.yend .< 0.222) .&& (basin.y0 .< 0.65))'
# temp = (circshift(temp, [1,0]) .&& circshift(temp, [-1,0])) .|| ((circshift(temp, [0,1]) .&& circshift(temp, [0,-1])))
temp = (circshift(temp, [1,0]) + circshift(temp, [-1,0]) + circshift(temp, [0,1]) + circshift(temp, [0,-1]))
points = pos[vec((temp .== 3)')]
points = points[sortperm([atan((pt - fp5)...) for pt in points])]
push!(points, points[1])
plgn = DataFrame(stack(points, dims = 1), :auto)
CSV.write("plgn.csv", plgn); plgn = CSV.read("plgn.csv", DataFrame)
push!(layer_, deepcopy(fisrt(layer_)))
plot!(layer_[4], first.(points), last.(points), lw = 2, color = :black)

@info "end of code"

### SINDy with external data
edata = dropmissing(data)
color_ecnm = get.(Ref(ColorSchemes.plasma), edata.ecnm)

push!(layer_, deepcopy(fisrt(layer_)))
for dat in eachrow(edata)
    plot!(layer_[5], [dat.yng[1], dat.yng+dat.dyng[1]], [dat.old[1], dat.old+dat.dold[1]]
    , arrow = arrow(:closed), color = get.(Ref(ColorSchemes.plasma), dat.ecnm[1]), α = 0.5)
end

g = SINDy(edata, [:dyng, :dold], [:yng, :old, :ecnm], N = 2)
print(g, ["y", "o", "e"])

dir_e = g.(pos .⊕ [0.66])
heatmap(reshape(log10.(norm.(dir_e)), length(0:.0005:xmax), length(0:.0005:ymax))', size = [600, 600], legend = :best)