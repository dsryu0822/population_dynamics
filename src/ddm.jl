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

@variables t x y
∂x(f) = Differential(x).(f) .|> expand_derivatives
∂y(f) = Differential(y).(f) .|> expand_derivatives

grid(tick::Float64, x_ = [0, 1], y_ = [0, 1]) = Base.product(x_[1]:tick:x_[2], y_[1]:tick:y_[2]) |> collect |> vec .|> collect
grid(tick::Int64, x_ = [0, 1], y_ = [0, 1]) = Base.product(LinRange(x_[1],x_[2],tick), LinRange(y_[1],y_[2],tick)) |> collect |> vec .|> collect
square(X) = reshape(X, isqrt(length(X)), isqrt(length(X)))
function solve(f, ic; tend = 10, h = 1e-2, dense = true)
    idx_end = length(1:h:tend)
    traj = zeros(2, idx_end); traj[:,1] = ic
    for tk in 2:(idx_end-1)
        if (!(0 ≤ traj[1,tk-1] ≤ xmax) || !(0 ≤ traj[2,tk-1] ≤ ymax)) ||
            (mod(tk, 100) == 0 && (norm(f(traj[:,tk-1])) < 1e-8))
            traj = traj[:,1:tk]
            break
        end
        traj[:,tk] = traj[:,tk-1] + h*f(traj[:,tk-1])
    end
    traj = collect(eachcol(traj))[1:(end-1)]
    return dense ? traj : traj[1:(idx_end ÷ 100):end]
end

@time include("datacall.jl")
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

const xmax = 0.7; const ymax = 0.7; resolution = 200
qargs = (; xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [800, 800], framestyle = :box, formatter = x -> "$(round(100x; digits = 2))%")
layer_ = [plot(; qargs...) for _ in 1:6]

### SINDy
setN = 6
odata = dropmissing(select(data, Not(:ecnm)))
# data[data.ISO3 .== "KOR", :]
sort(sort(odata, :Time), :ISO3)
f = SINDy(odata, [:dyng, :dold], [:yng, :old], N = setN)
print(f, ["y", "o"])
Θxy = Θ([x, y], N = setN)
∂xΘ = ∂x(Θxy)
∂yΘ = ∂y(Θxy)

### Find fixed points
pos = grid(resolution, [0, xmax], [0, ymax])
function eachfixed(f, pos)
    # slowpoints = pos[norm.(f.(pos)) .< 1e-3]
    slowpoints = [findfixed(f, maxitr = 1000, sp) for sp in pos[norm.(f.(pos)) .< 1e-3]]
    # scatter(first.(slowpoints), last.(slowpoints))
    dbscaned = dbscan(stack(slowpoints), 1e-2)
    fixedcandy = slowpoints[rand.(getproperty.(dbscaned.clusters, :core_indices))]
    @info "$(length(fixedcandy)) fixed points will be found"
    fixedpoint = [findfixed(f, fc; maxitr = 100000, atol = 1e-8, h = 1e-6) for fc in fixedcandy]
    idx_unfixed = findall(norm.(f.(fixedpoint)) .> 1e-8)
    unfixed = fixedpoint[idx_unfixed]
    gridsearch = [pairwiseadd(uf, grid(100, [-0.01, 0.01], [-0.01, 0.01])) for uf in unfixed]
    for (i, j) in enumerate(idx_unfixed)
        better = norm.(f.(gridsearch[i]))
        fixedpoint[j] .= gridsearch[i][argmin(better)]
    end
    return fixedpoint
end
fixedpoint = eachfixed(f, pos)
fixedpoint .=> norm.(f.(fixedpoint))
fp0 = [0.19888149603646155, 0.22155592768409718]

### Stability analysis
eigenvalues = []
for fp in fixedpoint
    xy = Dict((x, y) .=> fp)
    J = Symbolics.value.(substitute.([∂xΘ; ∂yΘ], Ref(xy)) * f.matrix)
    push!(eigenvalues, eigen(J).values)
end
encode_stability = sum.(map(x -> x .> 0, real.(eigenvalues))) .+ 1
stabiliary_colors = [:black, :gray, :red][encode_stability]
scatter!(layer_[2], first.(fixedpoint), last.(fixedpoint), mc = stabiliary_colors, ms = 10, msw = 0, txt = "        " .* string.(1:7))

### Separatrix
traj_ = []
layer_[3] = deepcopy(layer_[1])
# pos_basin = grid(10, [0, .60], [0, .35])
pos_basin = grid(100, [0, .7], [0, .7])
@showprogress for ic = pos_basin
    sol = solve(f, ic, tend = 10000, dense = true)
    push!(traj_, ic => last(sol))
    if (ic[1] > .5) || (ic[2] > .4) continue end
    plot!(layer_[3], first.(sol), last.(sol), lw = 1, color = :black, alpha = 0.1)
end
basin = DataFrame(x0 = first.(first.(traj_)), y0 = last.(first.(traj_)), xend = first.(last.(traj_)), yend = last.(last.(traj_)))
# CSV.write("basin.csv", basin); basin = CSV.read("basin.csv", DataFrame)

### Basin of attraction
layer_[4] = deepcopy(layer_[1])
temp = square((0.198.< basin.xend .< 0.199) .&& (0.220 .< basin.yend .< 0.222) .&& (basin.y0 .< 0.65))'
# scatter(first.(pos_basin[vec(temp)]), last.(pos_basin[vec(temp)]))
temp = (circshift(temp, [1,0]) + circshift(temp, [-1,0]) + circshift(temp, [0,1]) + circshift(temp, [0,-1]))
points = pos_basin[vec((temp .== 3)')]
# scatter(first.(points), last.(points))
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

png.(layer_, ["layer_[$(k)].png" for k in 1:6])

# candy_ = []
# @showprogress for ee = 0:0.01:1.0
#     dir_e = g.(pos .⊕ [ee])
#     candy = pos[log10.(norm.(dir_e)) .< -4]
#     append!(candy_, candy .⊕ [ee])
# end
# scatter(first.(candy_), getindex.(candy_, 2), last.(candy_); qargs..., color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)), xlims = [0.1, 0.5], ylims = [0.0, 0.2], zlabel = "economy", msw = 0)
# scatter(first.(candy_), getindex.(candy_, 2), color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)); qargs..., xlims = [0.1, 0.5], ylims = [0.0, 0.2], msw = 0)

### Some trajectories
layer_[6] = deepcopy(layer_[1]);
function trajectory!(inputplot, ic, tend)
    sol = solve(f, ic, tend = tend, dense = true)
    plot!(inputplot, first.(sol), last.(sol), lw = 1, color = :black, arrow = true);
end
trajectory!(layer_[6], [.20, .55], 1.1)
trajectory!(layer_[6], [.40, .50], 1.1)
trajectory!(layer_[6], [.12, .34], 5)
trajectory!(layer_[6], [.25, .34], 5)
trajectory!(layer_[6], [.13, .34], 5)
trajectory!(layer_[6], [.16, .34], 5)
trajectory!(layer_[6], [.50, .20], 5)
trajectory!(layer_[6], [.20, .11], 57)
trajectory!(layer_[6], [.30, .10], 75)
trajectory!(layer_[6], [.20, .10], 75)
trajectory!(layer_[6], [.06, .10], 10)
trajectory!(layer_[6], [.50, .10], 5)
trajectory!(layer_[6], [.45, .08], 300)
trajectory!(layer_[6], [.40, .05], 50)
trajectory!(layer_[6], [.60, .03], 5)
trajectory!(layer_[6], [.30, .03], 50)