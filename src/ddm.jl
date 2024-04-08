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

@time include("datacall.jl")
include("../../DataDrivenModel/src/DDM.jl")
include("../../setup/Portables.jl")
include("../../setup/Dynamics.jl")

const xmax = 0.5; const ymax = 0.5; resolution = 200
qargs = (; xlims = [0, xmax], ylims = [0, ymax], xlabel = "young", ylabel = "old", size = [800, 800], framestyle = :box, formatter = x -> "$(round(100x; digits = 2))%")

### SINDy
odata = dropmissing(select(data, Not(:ecnm)))
sort(sort(odata, :Time), :ISO3)
# data[data.ISO3 .== "KOR", :]



f_ = [SINDy(odata, [:dyng, :dold], [:yng, :old], N = setN) for setN in 1:20]
plot(log10.(maximum.(getproperty.(f_, :matrix))), xticks = 1:20, color = :black)
plot(size.(getproperty.(f_, :matrix), 1), xticks = 1:20, color = :black)
print(f_[20], ["y", "o"])
try mkdir("G:/population/ddm/quiver") catch end
for setN in 1:20
    cf = zip(odata[:, 3:4].yng, odata[:, 3:4].old)
    df = (f_[setN]).(eachrow(Matrix(odata[:, 3:4])))
    quiver(first.(cf), last.(cf), quiver = (first.(df), last.(df)), size = [800, 800])
    png("G:/population/ddm/quiver/N=$setN.png")
end

# setN = 6
fold10 = rand(1:10, nrow(odata))
MSE_fold10 = []
for setN = 1:20
    MSE_ = []
    for k in 1:10
        f = SINDy(odata[fold10 .!= k, :], [:dyng, :dold], [:yng, :old], N = setN)
        mse = sum(norm, eachrow(Matrix(odata[fold10 .== k, 5:6])) .- f.(eachrow(Matrix(odata[fold10 .== k, 3:4]))))
        push!(MSE_, mse)
        print("█")
    end
    MSE_ = push!(MSE_fold10, sum(MSE_)/10)
    println(last(MSE_fold10))
end
plot(MSE_fold10, xticks = 1:20, color = :black)


layer_ = [plot(; qargs...) for _ in 1:3]
for setN = 11:20
    mkdir("G:/population/ddm/N=$setN")
    # f = SINDy(odata, [:dyng, :dold], [:yng, :old], N = setN)

    ### Separatrix
    traj_ = []
    layer_[2] = deepcopy(layer_[1])
    # pos_basin = grid(10, [0, .60], [0, .35])
    pos_basin = grid(100, [0, .5], [0, .5])
    @showprogress for ic = pos_basin
        sol = solve(f_[setN], ic, tend = 1000, dense = true)
        push!(traj_, ic => last(sol))
        if !(0 < ic[1] < .5) || !(0 < ic[2] < .5) continue end
        plot!(layer_[2], first.(sol), last.(sol), lw = 1, color = :black, alpha = 0.1)
    end
    zero2end = DataFrame(x0 = first.(first.(traj_)), y0 = last.(first.(traj_)), xend = first.(last.(traj_)), yend = last.(last.(traj_)))
    CSV.write("G:/population/ddm/N=$(setN)/zero2end.csv", zero2end); # basin = CSV.read("basin.csv", DataFrame)

    layer_[3] = deepcopy(layer_[1])
    zero2end = zero2end[(0 .< zero2end.xend .< 0.5) .&& (0 .< zero2end.yend .< 0.5), :]
    scatter!(layer_[3], zero2end.xend, zero2end.yend, color = :black, msw = 0)
    # png(layer_[3], "layer_[3].png")

    ### SINDy with external data
    # edata = dropmissing(data)

    # layer_[5] = deepcopy(layer_[1])
    # color_ecnm = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, edata.ecnm)
    # for dat in eachrow(edata)
    #     plot!(layer_[5], [dat.yng[1], dat.yng+dat.dyng[1]], [dat.old[1], dat.old+dat.dold[1]]
    #     , arrow = arrow(:closed), color = get.(Ref(ColorSchemes.diverging_linear_bjr_30_55_c53_n256), dat.ecnm[1]), α = 0.5)
    # end
    # g = SINDy(edata, [:dyng, :dold], [:yng, :old, :ecnm], N = 6)
    # print(g, ["y", "o", "e"])

    # candy_ = []
    # @showprogress for ee = 0:0.01:1.0
    #     dir_e = g.(pos .⊕ [ee])
    #     candy = pos[log10.(norm.(dir_e)) .< -4]
    #     append!(candy_, candy .⊕ [ee])
    # end
    # scatter(first.(candy_), getindex.(candy_, 2), last.(candy_); qargs..., color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)), xlims = [0.1, 0.5], ylims = [0.0, 0.2], zlabel = "economy", msw = 0)
    # scatter(first.(candy_), getindex.(candy_, 2), color = get(ColorSchemes.diverging_linear_bjr_30_55_c53_n256, last.(candy_)); qargs..., xlims = [0.1, 0.5], ylims = [0.0, 0.2], msw = 0)

    png.(layer_, ["G:/population/ddm/N=$(setN)/layer_[$(k)].png" for k in 1:length(layer_)])
    rm("G:/population/ddm/N=$(setN)/layer_[1].png")
end