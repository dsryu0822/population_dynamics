include("../../DataDrivenModel/core/header.jl")
include("0_datacall.jl")

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                    This is 3D analysis

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
vrbl = [:dx, :dy, :dw], [:x, :y, :w]

data = data[:, [[:ISO3, :t]; [vrbl...;]]]
dropmissing!(data)
add_fold!(data, seed = 0)
TD_ = [data[data.fold .!= k, :] for k in sort(unique(data.fold))]
VD_ = [data[data.fold .== k, :] for k in sort(unique(data.fold))]
overM = unique(data.ISO3)[[sum(exp10.(Matrix(combine(groupby(data, :ISO3), last)[:, last(vrbl)])), dims = 2)...;] .> 1e+3]
data = data[data.ISO3 .∈ Ref(overM), :]

# Taic = []
# Vaic = []
# rankAnalysis = []
# d_ = 1:7
# for d = d_
#     cnfg = cook(last(vrbl), poly = 0:d)
#     ΘX = Θ(data[:, last(vrbl)], cnfg)
#     push!(rankAnalysis, [d, rank(ΘX), size(ΘX, 2)])

#     f_ = [SINDy(TD, vrbl, cnfg, λ = 0) for TD in TD_]
#     push!(Taic, mean([f.aic for f in f_]))
#     push!(Vaic, mean([AIC(f, VD) for (f, VD) in zip(f_, VD_)]))
# end

# plot(xticks = d_)
# # plot!(d_, Taic, label = "training AIC")
# plot!(d_, Vaic, label = "AIC")

# SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)])) |> print

# λ_ = exp10.(-4:0.1:4)
# λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]); λ).mse for λ in λ_]
# plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:4), legend = :none, color = :black)

# # f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:d_[argmin(Vaic)]), λ = exp10(-3.0))
# # f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:6), λ = exp10(-0.25))
# # f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:6), λ = exp10(-0.5))
# f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:6), λ = exp10(-0.5))
# f |> print

# data = deepcopy(_data)


rankAnalysis = []
d_ = 1:7
for d = d_
    cnfg = cook(last(vrbl), poly = 0:d)
    ΘX = Θ(data[:, last(vrbl)], cnfg)
    push!(rankAnalysis, [d, rank(ΘX), size(ΘX, 2)])
end

pxyargs = (; xlabel = "x", ylabel = "y", lims = [0, 6], legend = :topright)
pxzargs = (; xlabel = "x", ylabel = "w", lims = [0, 6], legend = :none)
pyzargs = (; xlabel = "y", ylabel = "w", lims = [0, 6], legend = :none)
xyzargs = (; xlabel = "x", ylabel = "y", lims = [0, 6], zlabel = "w", legend = :none)
vldn = DataFrame(zeros(length(unique(data.ISO3)), 7), ["d$d" for d = 1:7])
for d = 1:7
    try mkdir("figure/$d") catch end
    f = SINDy(data[data.t .≤ 2000, :], vrbl, cook(last(vrbl), poly = 0:d))
    lasterror = []
    @info "Validating d = $d"
    for DAT = groupby(data, :ISO3)
        # tend = (DAT.t[end] - DAT.t[1])
        # prdt = DataFrame(solve(f, [DAT[1, last(vrbl)]...], 0:0.01:tend), last(vrbl))
        tbegin = findfirst(DAT.t .== 2000)
        tend = DAT.t[end] - DAT.t[tbegin]
        prdt = DataFrame(solve(f, [DAT[tbegin, last(vrbl)]...], 0:0.01:tend), last(vrbl))
        push!(lasterror, sum(abs2, [prdt[end, :]...;] - [DAT[end, last(vrbl)]...;]))
        print(DAT.ISO3[1])
        # pt = "$(DAT.ISO3[1]): $(DAT.ISO3[1]), t ∈ [$(DAT.t[1]), $(DAT.t[end])]"
        # @info "$pt"

        pxy=plot(; pxyargs...)
            plot!(DAT.x, DAT.y, color = 1, label = :none); scatter!(DAT.x[[end]], DAT.y[[end]], color = 1, label = "data")
            plot!(prdt.x, prdt.y, color = 2, label = :none); scatter!(prdt.x[[end]], prdt.y[[end]], color = 2, label = "model")
        pyz=plot(; pyzargs...)
            plot!(DAT.y, DAT.w, color = 1); scatter!(DAT.y[[end]], DAT.w[[end]], color = 1)
            plot!(prdt.y, prdt.w, color = 2); scatter!(prdt.y[[end]], prdt.w[[end]], color = 2)
        pxz=plot(; pxzargs...)
            plot!(DAT.x, DAT.w, color = 1); scatter!(DAT.x[[end]], DAT.w[[end]], color = 1)
            plot!(prdt.x, prdt.w, color = 2); scatter!(prdt.x[[end]], prdt.w[[end]], color = 2)
        xyz=plot(; xyzargs...)
            plot!(DAT.x, DAT.y, DAT.w, color = 1); scatter!(DAT.x[[end]], DAT.y[[end]], DAT.w[[end]], color = 1)
            plot!(prdt.x, prdt.y, prdt.w, color = 2); scatter!(prdt.x[[end]], prdt.y[[end]], prdt.w[[end]], color = 2)
        plot(pxy, pyz, pxz, xyz, size = (800, 800), plot_title = "$(DAT.ISO3[1]) with (d=$d)")
        png("figure/$d/$(DAT.ISO3[1]).png")
    end
    println("$d done!")
    vldn[:, d] .= lasterror
end
CSV.write("validation_3D.csv", vldn)
vldn = CSV.read("validation_3D.csv", DataFrame)
pfmc = [mean(skipmissing(replace(col, NaN => missing))) for col in eachcol(vldn)]
plot(pfmc, yscale = :log10)

λ_ = exp10.(-4:0.1:4)
λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:2); λ).mse for λ in λ_]
plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:4), legend = :none, color = :black)

f = SINDy(data[data.t .< 2010, :], vrbl, cook(last(vrbl), poly = 0:2), λ = exp10(-2.5))
f |> print
CSV.write("final_model.csv", DataFrame(ISO3 = unique(data.ISO3), error = vldn[:, 2]))

λmse = [SINDy(data, vrbl, cook(last(vrbl), poly = 0:2); λ).mse for λ in λ_]
plot(λ_, λmse, xscale = :log10, xticks = exp10.(-4:4), legend = :none, color = :black)

f = SINDy(data, vrbl, cook(last(vrbl), poly = 0:2), λ = exp10(-2.5))
f |> print


# ic = [[dr...] for dr in eachrow(combine(groupby(data, :ISO3), row -> last(row))[:, [:x, :w, :y]])]
ic = vec([[dr...] for dr in Base.product(0:0.2:6, 0:0.2:6, 0:0.2:6)])

maxk = length(ic)
tspan = 0:0.01:1
trajx = zeros(maxk, length(tspan))
trajy = zeros(maxk, length(tspan))
trajw = zeros(maxk, length(tspan))
bit_valid = zeros(Bool, maxk)
# if isfile("data/trajx3.csv") && isfile("data/trajy3.csv") && isfile("data/trajz3.csv")
#     trajx = Matrix(CSV.read("data/trajx3.csv", DataFrame))'
#     trajy = Matrix(CSV.read("data/trajy3.csv", DataFrame))'
#     trajz = Matrix(CSV.read("data/trajz3.csv", DataFrame))'
# else
    @showprogress @threads for k in 1:maxk
        # if x + 1 < y continue end
        traj = solve(f, ic[k], tspan)
        if all(0 .< traj[end, :])
            trajx[k, :] .= traj[:, 1]
            trajy[k, :] .= traj[:, 2]
            trajw[k, :] .= traj[:, 3]
            bit_valid[k] = true
        end
    end
    trajx = trajx[bit_valid, :]
    trajy = trajy[bit_valid, :]
    trajw = trajw[bit_valid, :]
    CSV.write("data/trajx3.csv", DataFrame(trajx', :auto)[:, 1:2:end])
    CSV.write("data/trajy3.csv", DataFrame(trajy', :auto)[:, 1:2:end])
    CSV.write("data/trajw3.csv", DataFrame(trajw', :auto)[:, 1:2:end])
# end

frame = plot(xlabel = "x", ylabel = "y", zlabel = "w", legend = :none, size = [600, 600])

plt_vf = deepcopy(frame)
for k in axes(trajx, 1)
    plot!(plt_vf, trajx[k, :], trajy[k, :], trajw[k, :], color = :black, label = :none)
    scatter!(plt_vf, trajx[k, [end]], trajy[k, [end]], trajw[k, [end]], color = :black, ms = 1, alpha = .5)
end
# plt_vf
savefig("G:/population/pred50y.html")

Plots.plotly()
Plots.gr()

plot(vec(sum(trajz, dims = 1))/212)
plot(eachrow(trajz))
cumprod(1:7)

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                        ball trajectories

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
sargs = (; lims = [0, 8], ms = 1, msw = 0)

@async for k in 1:5
    tmp = [randn(3) for _ in 1:500]
    ic = (tmp ./ sqrt.(sum.(abs2, tmp))) .+ Ref([k, k, k])

    @time after10 = [last(eachrow(solve(f, i, 0:0.01:20))) for i in ic]
    @time after20 = [last(eachrow(solve(f, i, 0:0.01:100))) for i in ic]
    S00 = scatter(eachrow(stack(ic))..., label = "t = 0", color = 1; sargs...)
    S01 = scatter(eachrow(stack(after10))..., label = "t = 20", color = 2; sargs...)
    S02 = scatter(eachrow(stack(after20))..., label = "t = 100", color = 3; sargs...)
    S03 = scatter(eachrow(stack(ic))..., label = "t = 0"; sargs...)
    scatter!(eachrow(stack(after10))..., label = "t = 20"; sargs...)
    scatter!(eachrow(stack(after20))..., label = "t = 100"; sargs...)
    plot(S00, S01, S02, S03, size = (500, 500), plot_title = "IC = [$k, $k, $k] + noise", dpi = 180)
    png("dough_$k")
end

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                        fixed point

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""

fic = [f(p) for p in ic]
sum(abs2, f(ic[argmin(sum.(abs2, fic))]))

solve(f, ic[argmin(sum.(abs2, fic))], 0:0.01:100)
solve(f, [3, 3, 3], 0:0.01:100)

"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                        nullclines

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
ics = [[dr...] for dr in Base.product(0:0.05:6, 0:0.05:6, 0:0.05:6)]
tic = f.(ics)

xic = getindex.(tic, 1) .> 0
yic = getindex.(tic, 2) .> 0
wic = getindex.(tic, 3) .> 0

ncx = ics[(xic .!= circshift(xic, 1)) .&& (xic .!= circshift(xic, 2))]
ncy = ics[(yic .!= circshift(yic, 1)) .&& (yic .!= circshift(yic, 2))]
ncz = ics[(wic .!= circshift(wic, 1)) .&& (wic .!= circshift(wic, 2))]
candy = intersect(ncx, ncy, ncz)
sum.(abs2, f.(candy))

scatter(getindex.(ncx, 1), getindex.(ncx, 2), getindex.(ncx, 3), color = :red, ma = 0.1, label = "dx=0")
scatter!(getindex.(ncy, 1), getindex.(ncy, 2), getindex.(ncy, 3), color = :blue, ma = 0.1, label = "dy=0")
scatter!(getindex.(ncz, 1), getindex.(ncz, 2), getindex.(ncz, 3), color = :green, ma = 0.1, label = "dw=0")

after100 = [solve(f, ic, 0:0.01:1000) for ic in candy]


[last(eachrow(af100)) for af100 in after100]
plot(eachcol(after100[2])...)




"""'''''''''''''''''''''''''''''''''''''''''''''''''''''

                    reservoir computing

'''''''''''''''''''''''''''''''''''''''''''''''''''''"""
struct Reservoir
    M::Int64
    P::Int64
    N::Int64
    α::Float64
    β::Float64
    ρ::Float64
    A::AbstractMatrix
    Win::AbstractMatrix
    Wout::AbstractMatrix
    c::AbstractVector
    r::AbstractVector
end
function Base.show(io::IO, s::Reservoir)
    print(io, "f: R^$(s.M) → R^$(s.P)")
end
function print(s::Reservoir)
    println("""
    f: R^$(s.M) → R^$(s.P)
      reservoir size N = $(s.N)
       leakage reate α = $(s.α)
    sparsity of Wout β = $(s.β)
     spectral radius ρ = $(s.ρ)
    """)
    return nothing
end
function (s::Reservoir)(U)    
    encode(r, u) = (1-s.α)*r + tanh.(s.A * r + s.Win * u)
    decode(r, u) = s.Wout * encode(r, u) .+ s.c

    r_ = [s.r]
    [push!(r_, encode(r_[end], u)) for u in eachcol(U)]
    popfirst!(r_); R = stack(r_)
    S_pred = s.Wout * R .+ s.c
    return S_pred    
end
function reservoir_computing(U::AbstractMatrix, S::AbstractMatrix;
    seed = -1, N = 500, D = 2, α = 0.5, β = 1e-2, ρ = 1.0)

    if seed ≥ 0 Random.seed!(seed) end
    M = size(U, 1)
    P = size(S, 1)
    RN = erdos_renyi(N, (D*N ÷ 2))
    A = adjacency_matrix(RN) .* 2(rand(N, N) .- .5)
    A = ρ*A / maximum(sqrt.(sum.(abs2, eigen(Matrix(A)).values)))
    Win = sparse(stack([shuffle([1; zeros(M-1)]) for _ in 1:N])') .* 2(rand(N, M) .- .5)
    encode(r, u) = (1-α)*r + tanh.(A * r + Win * u)
    r_ = [randn(N)]
    [push!(r_, encode(r_[end], u)) for u in eachcol(U)]
    popfirst!(r_); R = stack(r_)
    Wout = S*R'inv(R*R' + (β*LinearAlgebra.I))
    c = -vec(Wout * mean(R, dims = 2) - mean(S, dims = 2))
    return Reservoir(M, P, N, α, β, ρ, A, Win, Wout, c, r_[end])
end

actl2022 = []
prdt2022 = []
for dat = groupby(data, :ISO3)
    iso3 = dat.ISO3[1]
    # U = Matrix(dat[:, [:x, :y]])'
    U = Matrix(dat[:, [:x]])'
    S = Matrix(dat[:, [:w]])'
    g = reservoir_computing(U[:, 1:60], S[:, 1:60]; seed = 0, ρ = 1, α = 1)
    prdS = g(U[:, 61:end])
    push!(actl2022, last(S))
    push!(prdt2022, last(prdS))
    # sqrt(mean(abs2, prdS - S[:, 61:end]))
    # plot(dat.t, vec(S), color = :black, label = "actual working population", title = iso3)
    # plot!(1950:2009, vec(g(U[:, 1:60])), color = :red, label = "fitted in training process")
    # plot!(2010:2022, vec(prdS), color = :red, ls = :dash, label = "predicted")
    # png("RC_figure/$(iso3).png")
end
actl2022
prdt2022

scatter(actl2022, prdt2022, lims = [2, 7], size = [400, 400], title = "(x,y) ↦ w at t = 2022, R²=$(round(100Rsq(prdt2022, actl2022), digits = 2)) %", xlabel = "actual w", ylabel = "predicted w", legend = :none)
