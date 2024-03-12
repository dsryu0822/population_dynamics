using ProgressMeter
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :SparseArrays, :Colors, :ColorSchemes, :JLD2, 
            :Clustering, :Symbolics]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm; cm = Plots.cm;
default(); default(legend = :none)

grid(tick::Float64, x_ = [0, 1], y_ = [0, 1]) = Base.product(x_[1]:tick:x_[2], y_[1]:tick:y_[2]) |> collect |> vec .|> collect
grid(tick::Int64, x_ = [0, 1], y_ = [0, 1]) = Base.product(LinRange(x_[1],x_[2],tick), LinRange(y_[1],y_[2],tick)) |> collect |> vec .|> collect
square(X) = reshape(X, isqrt(length(X)), isqrt(length(X)))
function solve(f, x; T = 10, h = 1e-2, shorten = true)
    x_ = [x]
    for t in 1:h:T
        if !(0 ≤ x_[end][1] ≤ xmax) || !(0 ≤ x_[end][2] ≤ ymax)
            break
        elseif rand() < 0.01
            if norm(f(x_[end])) < 1e-6
                break
            end
        end
        push!(x_, x_[end] + h*f(x_[end]))
    end    
    return shorten ? x_[1:max(1, (length(x_) ÷ 100)):end] : x_
end

@time include("datacall.jl")
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

### SINDy
odata = dropmissing(select(data, Not(:ecnm)))
# data[data.ISO3 .== "KOR", :]
sort(sort(odata, :Time), :ISO3)
f = SINDy(odata, [:dyng, :dold], [:yng, :old], N = 5)
print(f, ["y", "o"])

const xmax = 0.7; const ymax = 0.7; resolution = 200
qargs = (; xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [800, 800], framestyle = :box, formatter = x -> "$(round(100x; digits = 2))%")
layer_ = [plot(; qargs...) for _ in 1:5]

### Find fixed points
pos = grid(resolution, [0, xmax], [0, ymax])
dir = f.(pos)
fixedcandy = pos[norm.(dir) .< 1e-3]
dbscaned = dbscan(stack(fixedcandy), 0.01)
fixedcandy = fixedcandy[rand.(getproperty.(dbscaned.clusters, :core_indices))]
fixedpoint = [findfixed(f, fc; maxitr = 50000, atol = 1e-8) for fc in fixedcandy]
fixedpoint = [findfixed(f, fc; maxitr = 50000, atol = 1e-8, h = 1e-6) for fc in fixedpoint]
idx_unfixed = findall(norm.(f.(fixedpoint)) .> 1e-8)
unfixed = fixedpoint[idx_unfixed]
gridsearch = [pairwiseadd(uf, grid(100, [-0.01, 0.01], [-0.01, 0.01])) for uf in unfixed]
# gridsearch = [pairwiseadd(uf, [randn(2) for _ in 1:1000]) for uf in unfixed]
for (i, j) in enumerate(idx_unfixed)
    better = norm.(f.(gridsearch[i]))
    fixedpoint[j] .= gridsearch[i][argmin(better)]
end
fixedpoint .=> norm.(f.(fixedpoint))
fp0 = [0.19888149603646155, 0.22155592768409718]
# heatmap(square(log10.(norm.(dir)))')

### Stability analysis
@variables t x y
∂x(f) = Differential(x).(f) .|> expand_derivatives
∂y(f) = Differential(y).(f) .|> expand_derivatives
Θxy = Θ([x, y], N = 6)
∂xΘ = ∂x(Θxy)
∂yΘ = ∂y(Θxy)
eigenvalues = []
for j ∈ eachindex(fixedpoint)
    fp = (Dict(x => fixedpoint[j][1], y => fixedpoint[j][2]),)
    J = Float64.([substitute.(∂xΘ, fp) .|> Symbolics.value
                substitute.(∂yΘ, fp) .|> Symbolics.value]) * Matrix(f.matrix)
    push!(eigenvalues, eigen(J).values)
end
encode_stability = sum.(map(x -> x .> 0, real.(eigenvalues))) .+ 1
stabiliary_colors = [:black, :gray, :red][encode_stability]
scatter!(layer_[2], first.(fixedpoint), last.(fixedpoint), mc = stabiliary_colors, ms = 10, msw = 0, txt = "        " .* string.(1:7))

### Separatrix
traj_ = []
layer_[3] = deepcopy(layer_[1])
@showprogress for ic = grid(20, [0, .60], [0, .35])
    sol = solve(f, ic, T = 10000, shorten = false)
    push!(traj_, ic => last(sol))
    if any(sol[2] .< 0) continue end
    plot!(layer_[3], first.(sol), last.(sol), lw = 1, color = :black, alpha = 0.1)
end
png(layer_[3], "layer_[3].png")
basin = DataFrame(x0 = first.(first.(traj_)), y0 = last.(first.(traj_)), xend = first.(last.(traj_)), yend = last.(last.(traj_)))
# CSV.write("basin.csv", basin); basin = CSV.read("basin.csv", DataFrame)

### Some trajectories
layer_3 = deepcopy(layer_[1]);
sol = solve(f, [.20, .55], T = 1.1, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.40, .50], T = 1.1, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.12, .34], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.25, .34], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.13, .34], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.16, .34], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.50, .20], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.20, .11], T = 57, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.30, .10], T = 75, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.20, .10], T = 75, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.06, .10], T = 10, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.50, .10], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.45, .08], T = 300, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.40, .05], T = 50, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.60, .03], T = 5, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
sol = solve(f, [.30, .03], T = 50, shorten = false); plot!(layer_3, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
png(layer_3, "layer_3.png")

### Basin of attraction
temp = square((0.19.< basin.xend .< 0.222) .&& (0.19 .< basin.yend .< 0.222) .&& (basin.y0 .< 0.65))'
temp = (circshift(temp, [1,0]) + circshift(temp, [-1,0]) + circshift(temp, [0,1]) + circshift(temp, [0,-1]))
points = pos[vec((temp .== 3)')]
points = points[sortperm([atan((pt - fp0)...) for pt in points])]
push!(points, points[1])
plgn = DataFrame(stack(points, dims = 1), :auto) # CSV.write("plgn.csv", plgn);
plot!(layer_[4], first.(points), last.(points), lw = 2, color = :black)
# plgn = CSV.read("plgn.csv", DataFrame)
# plot!(layer_[4], plgn.x1, plgn.x2, lw = 2, color = :black)


### SINDy with external data
edata = dropmissing(data)

color_ecnm = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, edata.ecnm)
for dat in eachrow(edata)
    plot!(layer_[5], [dat.yng[1], dat.yng+dat.dyng[1]], [dat.old[1], dat.old+dat.dold[1]]
    , arrow = arrow(:closed), color = get.(Ref(ColorSchemes.diverging_linear_bjr_30_55_c53_n256), dat.ecnm[1]), α = 0.5)
end
g = SINDy(edata, [:dyng, :dold], [:yng, :old, :ecnm], N = 6)
print(g, ["y", "o", "e"])

candy_ = []
@showprogress for ee = 0:0.01:1.0
    dir_e = g.(pos .⊕ [ee])
    candy = pos[log10.(norm.(dir_e)) .< -4]
    append!(candy_, candy .⊕ [ee])
end
scatter(first.(candy_), getindex.(candy_, 2), last.(candy_); qargs..., color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)), xlims = [0.1, 0.5], ylims = [0.0, 0.2], zlabel = "economy", msw = 0)
scatter(first.(candy_), getindex.(candy_, 2), color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)); qargs..., xlims = [0.1, 0.5], ylims = [0.0, 0.2], msw = 0)
