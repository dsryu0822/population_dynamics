using ProgressMeter
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :SparseArrays, :Colors, :ColorSchemes, :JLD2, 
            :Clustering, :Symbolics]
@showprogress for package in packages
    @eval using $(package)
end; @info "packages $(packages) are loaded"
mm = Plots.mm; cm = Plots.cm;
default(); default(legend = :none)
using Base.Threads: @threads


function grid(tick::Int64, xyz...)
    return Base.product(LinRange.(first.(xyz), last.(xyz), tick)...) |> collect |> vec .|> collect
end
function solve(f, ic; tend = 10, h = 1e-2, dense = true)
    idx_end = length(1:h:tend)
    traj = zeros(length(ic), idx_end); traj[:,1] .= ic
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

# function eachfixed(f, pos)
#     # slowpoints = pos[norm.(f.(pos)) .< 1e-3]
#     slowpoints = [findfixed(f, maxitr = 1000, sp) for sp in pos[norm.(f.(pos)) .< 1e-3]]
#     # scatter(first.(slowpoints), last.(slowpoints))
#     dbscaned = dbscan(stack(slowpoints), 1e-2)
#     fixedcandy = slowpoints[rand.(getproperty.(dbscaned.clusters, :core_indices))]
#     @info "$(length(fixedcandy)) fixed points will be found"
#     fixedpoint = [findfixed(f, fc; maxitr = 100000, atol = 1e-8, h = 1e-6) for fc in fixedcandy]
#     idx_unfixed = findall(norm.(f.(fixedpoint)) .> 1e-8)
#     unfixed = fixedpoint[idx_unfixed]
#     gridsearch = [pairwiseadd(uf, grid(100, [-0.01, 0.01], [-0.01, 0.01])) for uf in unfixed]
#     for (i, j) in enumerate(idx_unfixed)
#         better = norm.(f.(gridsearch[i]))
#         fixedpoint[j] .= gridsearch[i][argmin(better)]
#     end
#     return fixedpoint
# end

@time include("datacall.jl")
include("../../DataDrivenModel/core/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

const xmax = 50; const ymax = 50; const zmax = 100; resolution = 200
qargs = (; xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [800, 800], framestyle = :box, formatter = x -> "$(round(x; digits = 2))%")
maxN = 20

### SINDy
_odata = dropmissing(select(data, Not(:ecnm)))
_odata[:, 3:end] .= 100Matrix(_odata[:, 3:end])
odata = sort(sort(_odata, :Time), :ISO3)
# data[data.ISO3 .== "KOR", :]

@time f_ = [SINDy(odata, [:dyng, :dold], [:yng, :old], N = setN) for setN in 1:maxN];
plot(log10.(maximum.(getproperty.(f_, :matrix))), xticks = 1:maxN, color = :black)
plot(size.(getproperty.(f_, :matrix), 1), xticks = 1:maxN, color = :black)
print(f_[7], ["y", "o"])

function eachfixed(f, rsln, pos = grid(rsln, [0, xmax], [0, ymax], [0, 100]))
    df = f.(pos)
    dim = length(first(df))
    onechain = [rsln for _ in 1:dim]; onechain[end] = 1

    null_ = []
    for k in 1:dim
        _df = reshape(getindex.(df, k), [rsln for _ in 1:dim]...) .< 0
    
        null = []
        for _k = 1:dim
            zero_pad = circshift(onechain, _k)
            push!(null, cat(zeros(Bool, zero_pad...), .!iszero.(diff(_df, dims = _k)), dims = _k))
        end
        push!(null_, reduce(.|, null))
    end
    return pos[vec(reduce(.&, null_))]
end
@time temp = eachfixed(e8, 501)
CSV.write("e8_fixed.csv", DataFrame(stack(temp)', [:y, :o, :e]))

scatter(eachrow(stack(temp))..., msw = 0, ms = 1)
pos = grid(10, [0, xmax], [0, ymax], [0, 100])

try mkdir("G:/population/ddm/_nullcline") catch end
for setN = 1:maxN
    @time df = f_[setN].(pos);
    null_x = [zeros(Bool, 1, resolution); .!iszero.(diff(reshape(first.(df), resolution, resolution) .< 0, dims = 1))] .|| [zeros(Bool,    resolution)  .!iszero.(diff(reshape(first.(df), resolution, resolution) .< 0, dims = 2))];
    null_y = [zeros(Bool,    resolution)  .!iszero.(diff(reshape( last.(df), resolution, resolution) .< 0, dims = 2))] .|| [zeros(Bool, 1, resolution); .!iszero.(diff(reshape( last.(df), resolution, resolution) .< 0, dims = 1))];
    null__ = null_x .&& null_y;

    plot(title = "Nullcline for N=$setN", legend = :topright; qargs...)
    scatter!(first.(pos[vec(null_x)]), last.(pos[vec(null_x)]), ms = .5, msw = 0, label = "dx = 0")
    scatter!(first.(pos[vec(null_y)]), last.(pos[vec(null_y)]), ms = .5, msw = 0, label = "dy = 0")
    scatter!(first.(pos[vec(null__)]), last.(pos[vec(null__)]), color = :black, shape = :x, label = "dy = dx = 0")
    png("G:/population/ddm/_nullcline/N=$setN.png")
end

# try mkdir("G:/population/ddm/_quiver") catch end
# for setN in 1:maxN
#     cf = zip(odata[:, 3:4].yng, odata[:, 3:4].old)
#     df = (f_[setN]).(eachrow(Matrix(odata[:, 3:4])))
#     quiver(first.(cf), last.(cf), quiver = (first.(df), last.(df)), size = [800, 800], color = :black)
#     png("G:/population/ddm/_quiver/N=$setN.png")
# end

# MSE_loo = []
# for setN = 1:20
#     MSE_ = []
#     for k in axes(odata, 1)
#         f = SINDy(odata[Not(k), :], [:dyng, :dold], [:yng, :old], N = setN)
#         mse = norm(collect(odata[k, 5:6]) .- f(collect(odata[k, 3:4])))
#         push!(MSE_, mse)
#         iszero(mod(k, 1000)) && print("█")
#     end
#     push!(MSE_loo, sum(MSE_)/nrow(odata))
#     println(last(MSE_loo))
# end
# plot(MSE_loo, xticks = 1:20, color = :black)
# scatter!((1:20)[[argmin(MSE_loo)]], [minimum(MSE_loo)])

try mkdir("G:/population/ddm/_separatrix") catch end
@threads for setN = 1:20
    ### Separatrix
    
    traj_ = []
    pos_basin = grid(50, [0, xmax], [0, ymax])
    plt_sptr = plot(; qargs...)
    @showprogress for ic = pos_basin
        sol = solve(f_[setN], ic, tend = 1000, dense = true)
        push!(traj_, ic => last(sol))
        if !(0 < ic[1] < xmax) || !(0 < ic[2] < ymax) continue end
        plt_sptr = plot(plt_sptr, first.(sol), last.(sol), lw = 1, color = :black, alpha = 0.1)
    end
    zero2end = DataFrame(x0 = first.(first.(traj_)), y0 = last.(first.(traj_)), xend = first.(last.(traj_)), yend = last.(last.(traj_)))
    CSV.write("G:/population/ddm/_separatrix/N=$(setN) zero2end.csv", zero2end); # basin = CSV.read("basin.csv", DataFrame)

    zero2end = zero2end[(0 .< zero2end.xend .< xmax) .&& (0 .< zero2end.yend .< ymax), :]
    scatter!(plt_sptr, zero2end.xend, zero2end.yend, color = :black, msw = 0, shape = :x)
    png(plt_sptr, "G:/population/ddm/_separatrix/N=$(setN).png")
end



_edata = dropmissing(data)
_edata[:, 3:end] .= 100Matrix(_edata[:, 3:end])
edata = sort(_edata, [:ISO3, :Time])
rename!(edata, :ecnm => :e, :decnm => :de, :yng => :y, :old => :o, :dyng => :dy, :dold => :do)

p1 = plot(; qargs...)
color_ecnm = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, edata.e)
for dat in eachrow(edata[edata.ISO3 .== "KOR", :])
    plot!(p1, [dat.y[1], dat.y+dat.dy[1]], [dat.o[1], dat.o+dat.do[1]]
    , arrow = arrow(:closed), color = get.(Ref(ColorSchemes.diverging_linear_bjr_30_55_c53_n256), dat.e[1]))
end
plot(p1, size = [350, 350], formatter = x -> "", xlabel = "", ylabel = "")
png("temp")

e5 = SINDy(edata, [:dy, :do, :de], [:y, :o, :e], N = 5); print(e5)
e8 = SINDy(edata, [:dy, :do, :de], [:y, :o, :e], N = 8); print(e8)

qargs = (; xlabel = "y", ylabel = "o", zlabel = "e", xlims = [0, xmax], ylims = [0, ymax], zlims = [0, 100])
p1 = plot(; qargs...)
@showprogress for ic = Base.product(0:10:50, 0:10:50, 0:10:100)
    ic = [50, 5, 10]
    temp = stack(solve(e5, ic, tend = 100))
    plot!(p1, eachrow(temp)..., color = :red)
    plot!(p1, eachrow(temp[[1,2],:])..., zeros(size(temp, 2)), color = :black)
end
p1
png("temp1")

resolution = 1000
pos = grid(resolution, [0, xmax], [0, ymax]);
candidates = DataFrame(y = [], o = [], e = [])
try mkdir("G:/population/ddm/ncln_e5") catch end
@showprogress for _e = 0:100
    filename = "G:/population/ddm/ncln_e5/e=$(_e).png"
    if isfile(filename) continue end

    df = e5.(pos .⊕ _e);
    null_x = [zeros(Bool, 1, resolution); .!iszero.(diff(reshape(first.(df), resolution, resolution) .< 0, dims = 1))] .|| [zeros(Bool,    resolution)  .!iszero.(diff(reshape(first.(df), resolution, resolution) .< 0, dims = 2))];
    null_y = [zeros(Bool,    resolution)  .!iszero.(diff(reshape( last.(df), resolution, resolution) .< 0, dims = 2))] .|| [zeros(Bool, 1, resolution); .!iszero.(diff(reshape( last.(df), resolution, resolution) .< 0, dims = 1))];
    null__ = null_x .&& null_y;

    plot(title = "Nullcline with e = $(_e)", legend = :topright; qargs...)
    scatter!(first.(pos[vec(null_x)]), last.(pos[vec(null_x)]), ms = .5, msw = 0, label = "dx = 0")
    scatter!(first.(pos[vec(null_y)]), last.(pos[vec(null_y)]), ms = .5, msw = 0, label = "dy = 0")
    scatter!(first.(pos[vec(null__)]), last.(pos[vec(null__)]), color = :black, shape = :x, label = "dy = dx = 0")

    append!(candidates, DataFrame(y = first.(pos[vec(null__)]), o = last.(pos[vec(null__)]), e = fill(_e, count(null__))))
    CSV.write("G:/population/ddm/ncln_e5/candidates.csv", candidates)
    png("G:/population/ddm/ncln_e5/e=$(_e).png")
end

# --------------------------------------------
setN = 11

concolor = get(colorschemes[:rainbow], odata[odata.Time .== 1950, :old] / maximum(odata[odata.Time .== 1950, :old]))
background = deepcopy(plt_sptr);
for T = 1950:2022
# anime = @showprogress @animate for T = 1950:2022
    temp = odata[odata.Time .== T, :]
    # scatter(background, temp.yng, temp.old, color =:white);
    scatter(background, temp.yng, temp.old, color = concolor, msw = 0)
    png("$T.png")
end
# mp4(anime, "7.mp4", fps = 10)

10