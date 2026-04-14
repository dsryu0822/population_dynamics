include("../../DataDrivenModel/core/header.jl")
include("0_datacall.jl")

temp = data[data.t .== 2023, [:ISO3, :x, :w, :y]]
ISO3161 = [temp.ISO3[sum(exp10.(Matrix(temp[:, [:x, :w, :y]])), dims = 2) .> 1e3]...]
_data = dropmissing(data0[data0.ISO3 .∈ Ref(ISO3161), [:ISO3, :t, :X, :Y, :dX, :dY]])

vrbl = ([:dX, :dY], [:X, :Y])
XY_2020 = _data[_data.t .== 2020, :]
XY_2022 = _data[_data.t .== 2022, :]

groupby(_data, :ISO3)[1]

results = DataFrame(d = [], sse = [])
for d in 1:4 # NaN error for d = 9
    cnfg = cook(last(vrbl), poly = 0:d)
    f = SINDy(_data[_data.t .≤ 2020, :], vrbl, cnfg)
    error_abs2 = zeros(161)
    @time for k in 1:161
        prdt = solve(f, [XY_2020[k, [:X, :Y]]...;], 0:1e-2:2)[end, :]
        actl = [XY_2022[k, [:X, :Y]]...;]
        error_abs2[k] = sum(abs2, actl - prdt)
    end
    push!(results, (d, sum(error_abs2)))
end
results
plot(results.d, results.sse, legend = :none, color = :black, xlabel = "d")
cnfg = cook(last(vrbl), poly = 0:4)
f = SINDy(_data, vrbl, cnfg)
λ_ = logrange(1e-4, 1e+0, 10)
plot(λ_, [SINDy(_data, vrbl, cnfg; λ).mse for λ in λ_], xscale = :log10)
f = SINDy(_data, vrbl, cnfg, λ = 1e-2)
f |> print

ic_ = vec([[x...] for x in Base.product(0:0.001:0.5, 0:0.001:0.5)])

fic_ = f.(ic_)
_X = reshape(first.(fic_), 501, :) .> 0
_Y = reshape(last.(fic_), 501, :) .> 0

null_X = (_X .!= circshift(_X, (1, 0)) .|| (_X .!= circshift(_X, (-1, 0)))); null_X[1, :] .= false; null_X[end, :] .= false; null_X[:, 1] .= false; null_X[:, end] .= false; 
null_Y = (_Y .!= circshift(_Y, (0, 1)) .|| (_Y .!= circshift(_Y, (0, -1)))); null_Y[1, :] .= false; null_Y[end, :] .= false; null_Y[:, 1] .= false; null_Y[:, end] .= false; 
plt_nc = plot(legend = :none, lims = [0.01, .5], size = [800, 800])
heatmap!(plt_nc, 0:0.001:0.5, 0:0.001:0.5, null_X', color = [:transparent, :blue])
heatmap!(plt_nc, 0:0.001:0.5, 0:0.001:0.5, null_Y', color = [:transparent, :red])
# scatter!(XY_2020.X, XY_2020.Y)
scatter(plt_nc, XY_2022.X, XY_2022.Y, text = ISO3161, ms = 0, size = [1600, 1600])
candy = ic_[findall(vec(null_X) .&& vec(null_Y))]
f.(candy)

f.(eachrow([XY_2022.X XY_2022.Y]))
# quiver(100reshape(first.(ic_), 501, :)[1:10:end, 1:10:end], 100reshape(last.(ic_), 501, :)[1:10:end, 1:10:end], quiver = (reshape(first.(fic_), 501, :)[1:10:end, 1:10:end]*2, reshape(last.(fic_), 501, :)[1:10:end, 1:10:end]*2), color = :black, size = [800, 800])

for onecountry in  groupby(_data, :ISO3)
    plot!(plt_nc, onecountry.X, onecountry.Y, arrow = true, color = :black, alpha = 0.1)
end

eigen(jacobian(Function, f)(fp)).values
f(fp)
fp = [mean(first.(candy)), mean(last.(candy))]
for _ in 1:100
    randirec = randn(2)
    randirec /= norm(randirec)

    pX, pY = eachcol(solve(f, fp + 0.0001randirec, 0:1e-2:5))
    plot!(plt_nc, pX, pY, arrow = true, color = :black)
end
plot(plt_nc, lims = [])


# k_ = []
for _ in 1:300
    k = rand(eachindex(ic_))
    pX, pY = eachcol(solve(f, ic_[k], 0:1e-2:10))
    plot!(plt_nc, pX, pY, arrow = true, color = :black)
end
plt_nc
png(plt_nc, "quiver.png")
push!(k_, k)
pwd()