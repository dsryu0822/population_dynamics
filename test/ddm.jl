using ProgressMeter
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :SparseArrays, :Colors, :ColorSchemes]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm
default(); default(legend = :none)

grid(tick, x_ = [0, 1], y_ = [0, 1]) = Base.product(x_[1]:tick:x_[2], y_[1]:tick:y_[2]) |> collect |> vec .|> collect

function solve(f, x; T = 10, h = 1e-2, shorten = true)
    x_ = [x]
    for t in 1:h:T
        if !(0 ≤ x_[end][1] ≤ xmax) || !(0 ≤ x_[end][2] ≤ ymax)
            break
        elseif rand() < 0.1
            if norm(f(x_[end])) < 1e-6
                break
            end
        end
        push!(x_, x_[end] + h*f(x_[end]))
    end    
    return shorten ? x_[1:(length(x_) ÷ 100):end] : x_
end

# @time include("datacall.jl")
data = CSV.read("cached.csv", DataFrame)
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

n = 6
const xmax = 0.6
const ymax = 0.35
resolution = 1000
qargs = (; title = "N = $n", xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [600, 600], formatter = x -> "$(round(100x; digits = 2))%")
# cxcy = cat(eachslice(head(tnsr_yo), dims = 3)..., dims = 2)
# dxdy = cat(eachslice(diff(tnsr_yo, dims = 3), dims = 3)..., dims = 2)
# quiver(cxcy[1, 1:2:end], cxcy[2, 1:2:end], quiver = (dxdy[1, 1:2:end], dxdy[2, 1:2:end]), α = 0.1, c = :black; qargs...)

# data[data.ISO3 .== "KOR", :]
odata = dropmissing(select(data, Not(:ecnm)))
sort(sort(odata, :Time), :ISO3)
f = SINDy(odata, [:dyng, :dold], [:yng, :old], N = 6)
# f = SINDy(cxcy', dxdy', K = k)
print(f, ["y", "o"])
# print(f, ["y", "o", "e"])

pos = grid(0.0005, [0, xmax], [0, ymax])
dir = f.(pos)

arg0x = findall(abs.(getindex.(dir, 1)) .< 0.0001)
arg0y = findall(abs.(getindex.(dir, 2)) .< 0.0001)

ncln = plot(; qargs...)
scatter!(ncln, first.(pos[arg0x]), last.(pos[arg0x]), color = :red, alpha = 0.1, ms = 1, msw = 0;)
scatter!(ncln, first.(pos[arg0y]), last.(pos[arg0y]), color = :blue, alpha = 0.1, ms = 1, msw = 0;)
png("23")

using Clustering
using StatsBase
fixedcandy = pos[arg0x ∩ arg0y]
scatter(ncln, first.(fixedcandy), last.(fixedcandy), color = :black, shape = :x)
png("24")
dbscaned = dbscan(stack(fixedcandy), 0.01)
fixedcandy = fixedcandy[rand.(getproperty.(dbscaned.clusters, :core_indices))]
# norm.(f.(fixedcandy))
fixedpoint = [findfixed(f, fc; maxitr = 50000, atol = 1e-8) for fc in fixedcandy]

idx_unfixed = findall(norm.(f.(fixedpoint)) .> 1e-8)
unfixed = fixedpoint[idx_unfixed]
gridsearch = [pairwiseadd(uf, grid(0.0001, [-0.01, 0.01], [-0.01, 0.01])) for uf in unfixed]

for (i, j) in enumerate(idx_unfixed)
    better = norm.(f.(gridsearch[i]))
    better[argmin(better)]
    fixedpoint[j] .= gridsearch[i][argmin(better)]
end
fixedpoint .=> norm.(f.(fixedpoint))

# fixedpoint = vec.(mean.([pos[:, indices] for indices in getproperty.(dbscaned.clusters, :core_indices)], dims = 2))

using Symbolics
@variables t x y
∂x = Differential(x)
∂y = Differential(y)
Θxy = Θ([x, y], N = 6)
∂xΘ = ∂x.(Θxy) .|> expand_derivatives
∂yΘ = ∂y.(Θxy) .|> expand_derivatives

encode_stability = []
for j ∈ eachindex(fixedpoint)
fp = (Dict(x => fixedpoint[j][1], y => fixedpoint[j][2]),)
J = Float64.([substitute.(∂xΘ, fp) .|> Symbolics.value
              substitute.(∂yΘ, fp) .|> Symbolics.value]) * Matrix(f.matrix)
push!(encode_stability, sum(real.(eigen(J).values) .> 0))
end
stabiliary_colors = [:black, :gray, :red][encode_stability .+ 1]
# stabiliary_styles = [:solid, :dot, :dash][encode_stability .+ 1]

θ = 0:1:359
r = .0001

ring = hcat([fp .+ [(r .* cosd.(θ)) (r .* sind.(θ))]' for fp in fixedpoint]...)
# ring = Base.product(0:.01:xmax, 0:.01:ymax) |> collect .|> collect |> vec |> stack
_ncln = deepcopy(ncln)
@showprogress for k in axes(ring, 2)
    traj_ = solve(f, ring[:, k], T = 10)
    plot!(_ncln, first.(traj_), last.(traj_), arrow = true, lw = 1, color = :black, alpha = 0.1)
    # scatter!(_ncln, first.(traj_)[[end]], last.(traj_)[[end]], shape = :x, color = :black, alpha = 0.1)
    # for kk in 1:10
    #     traj_ = solve(f, traj_[end], T = 10000)
    #     plot!(_ncln, first.(traj_), last.(traj_), arrow = true, alpha = 0.2, color = :black)
    # end
end
_ncln
scatter!(_ncln, first.(fixedpoint), last.(fixedpoint), mc = stabiliary_colors, ms = 10, msw = 0, txt = "       " .* string.(1:7))
png("22")

_jplot = []
for j = 1:7
xlim, ylim = collect(zip(fixedpoint[j] .- 1.5r, fixedpoint[j] .+ 1.5r))
push!(_jplot, plot(_ncln, xlims = xlim, ylims = ylim, title = "j = $j"))
end
jplot = plot(_jplot..., layout = (2, 4), size = 200*[16, 9])
png("23")
run(`explorer $(pwd())`)

heatmap(rand([:a, :b, :c], 512,512))

tick = 0.001
ic_ = grid(tick, [0, xmax], [0, ymax])
Z_ = []
# ic = rand(ic_)
@showprogress for ic = ic_
    code = 0
    endpoint = last(solve(f, ic, T = 1000, shorten = false))
    if endpoint[1] ≤ 0
        code = 1
    elseif endpoint[2] ≤ 0
        code = 2
    elseif endpoint[1] > xmax
        code = 3
    elseif endpoint[2] > ymax
        code = 4
    elseif norm(f(endpoint)) < 1e-6
        code = 5
    else
        code = 6
    end
    push!(Z_, ic => code)
    # push!(Z_, ic => endpoint)
end
CSV.write("basin-3.csv", DataFrame(x0 = first.(first.(Z_)), y0 = last.(first.(Z_)), code = last.(Z_)))
heatmap(reshape(last.(Z_), length(0:tick:xmax), length(0:tick:ymax))',
color = [:red, :blue, :pink, :skyblue, :black, :white], size = [600, 600])
png("basin-3")

# reshape(first.(Z_), length(0:tick:xmax), length(0:tick:ymax))
# cp = palette([:red, :blue, :pink, :skyblue, :black, :white])
# scatter(first.(first.(Z_)), last.(first.(Z_)), color = cp[last.(Z_)], shape = :sq, msw = 0, ms = 5, size = [600, 600])
# heatmap(0:tick:xmax, 0:tick:ymax, reshape(last.(Z_), length(0:tick:xmax), length(0:tick:ymax)))
# heatmap(0:.01:xmax, 0:.01:ymax, reshape(last.(Z_), length(0:.01:xmax), length(0:.01:ymax)))

@info "end of code"
# length(65:101)
# POPy[2020][63]
# plot(xlabel = "Age", ylabel = "Population (k)", title = "Natural age distribution", xlims = [0, 101], ylims = [0, Inf], xticks = [0, 15, 65])
# plot!(0:15, POPy[2020][67].PopTotal[1:16], lw = 1, color = 1, fill = 0)
# plot!(64:100, POPy[2020][67].PopTotal[65:end], lw = 1, color = 2, fill = 0)
# plot!(0:100, POPy[2020][67].PopTotal, lw = 2, color = :black)
# png("17")
# tnsr_yo[:, 63, 71]
# tnsr_yo[:, 63, 72]

edata = dropmissing(data)
color_ecnm = get.(Ref(ColorSchemes.plasma), edata.ecnm)

quiver(edata.yng, edata.old, quiver = (edata.dyng, edata.dold), α = 0.1, lw = 2, color = :black; qargs...)
scatter!(edata.yng, edata.old, edata.ecnm, α = 0.1, color = color_ecnm, msw = 0, ms = 5)

p1 = plot(; qargs...)
for dat in eachrow(edata)
    plot!(p1, [dat.yng[1], dat.yng+dat.dyng[1]], [dat.old[1], dat.old+dat.dold[1]]
    , arrow = arrow(:closed), color = get.(Ref(ColorSchemes.plasma), dat.ecnm[1]), α = 0.5)
end
p1