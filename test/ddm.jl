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

function solve(f, x; T = 1000)
    x_ = [x]
    for t in 1:T
        push!(x_, x_[end] + 0.01f(x_[end]))
    end
    return x_[1:100:end]
end

@time include("datacall.jl")
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")

cxcy = cat(eachslice(head(tnsr_yo), dims = 3)..., dims = 2)
dxdy = cat(eachslice(diff(tnsr_yo, dims = 3), dims = 3)..., dims = 2)

k = 2
qargs = (; title = "k = $k", xlims = [0, 0.5], ylims = [0, 0.3], xlabel = "young", ylabel = "old", size = [600, 600])
f = SINDy(cxcy', dxdy', λ = 0.01, K = k)

pos = Base.product(0:.01:.55, 0:.01:0.35) |> collect .|> collect |> vec |> stack
dir = f.(eachcol(pos)) |> stack
quiver(pos[1, :], pos[2, :], quiver = (dir[1, :], dir[2, :]), color = :black; qargs...)
png("240124_k=$k vectorfield.png")

arg0x = findall(abs.(dir[1,:]) .< 0.00005)
arg0y = findall(abs.(dir[2,:]) .< 0.00005)

ncln = plot(; qargs...)
scatter!(ncln, pos[1,arg0x], pos[2,arg0x], color = :red, alpha = 0.5, ms = 1, msw = 0;)
scatter!(ncln, pos[1,arg0y], pos[2,arg0y], color = :blue, alpha = 0.5, ms = 1, msw = 0;)

_ncln = deepcopy(ncln)
fixedcandy = pos[:, abs.(dir[1,:]) .< 0.00001 .&& abs.(dir[2,:]) .< 0.00001]
fixedpoint = fixedcandy |> eachcol .|> collect

# traj_ = solve(f, fp, T = 10000)
# plot!(_ncln, first.(traj_), last.(traj_), arrow = true, color = :black)

fp = fixedpoint[3]
θ = 0:10:359
ring = fp .+ [0.01cosd.(θ) 0.01sind.(θ)]'
@showprogress for k in eachindex(θ)
    traj_ = solve(f, ring[:, k], T = 10000)
    plot!(_ncln, first.(traj_), last.(traj_), arrow = true, alpha = 0.2, color = :black)
    for kk in 1:10
        traj_ = solve(f, traj_[end], T = 10000)
        plot!(_ncln, first.(traj_), last.(traj_), arrow = true, alpha = 0.2, color = :black)
    end
end
_ncln
png(_ncln, "240124_k=$k.png")