include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2]])
include(joinpath(@__DIR__, "06_analyze.jl"))
using Flux
@load joinpath(@__DIR__, "../data.jld2") ISO3 iso3

vrbl = data_config.vars
diso3 = population_data(iso3; data_config...)
data3 = deepcopy(diso3)
X, Y = map_data(diso3; vars = vrbl)
μX, σX = vec(mean(X, dims = 2)), max.(vec(std(X, dims = 2)), 1e-8)
μY, σY = vec(mean(Y, dims = 2)), max.(vec(std(Y, dims = 2)), 1e-8)
Xn, Yn = Float32.((X .- μX) ./ σX), Float32.((Y .- μY) ./ σY)

Random.seed!(0)
model = Chain(Dense(size(X, 1) => 32, tanh), Dense(32 => 32, tanh), Dense(32 => size(Y, 1)))
optimizer = Flux.setup(Adam(1e-3), model)
EPOCH = 5000
L3_ = [mlp_loss(model, Xn, Yn)]
best = deepcopy(model)
blas_threads = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    @showprogress for epoch in 1:EPOCH
        Flux.train!(mlp_loss, model, [(Xn, Yn)], optimizer)
        L = mlp_loss(model, Xn, Yn)
        if L < minimum(L3_)
            global best = deepcopy(model)
        end
        push!(L3_, L)
    end
finally
    BLAS.set_num_threads(blas_threads)
end

# Freeze scaling and use Float64 for finite-difference stability calculations.
model3 = Flux.f64(best)
f3 = let μX = copy(μX), σX = copy(σX), μY = copy(μY), σY = copy(σY)
    Chain(x -> (x .- μX) ./ σX, model3, y -> y .* σY .+ μY)
end
pfmc3 = DataFrame(variable = vrbl,
    training_rmse = vec(sqrt.(mean(abs2, f3(X) - Y, dims = 2))))
display(pfmc3)
plt_loss3 = plot(0:EPOCH, L3_, yscale = :log10,
    xlabel = "epoch", ylabel = "normalized training MSE", label = :none)

pred3 = Dict(k => forecast(f3, Matrix(v[1:2, vrbl]), nrow(v)-2; discrete = true, lags = 1)
    for (k, v) in diso3)
k = haskey(diso3, "KOR") ? "KOR" : first(sort(collect(keys(diso3))))
plt_prediction3 = prediction_plot(diso3[k][3:end, :], ["MLP" => pred3[k]]; vars = vrbl, title = k)
display(plt_prediction3)

data = vcat([diso3[k] for k in sort(collect(keys(diso3)))]...)
fixedpoints3 = find_fixedpoints(f3, data; vars = vrbl, discrete = true, lags = 1)
display(fixedpoints3)

output_dir = joinpath(@__DIR__, "../archive")
mkpath(output_dir)
JLD2.@save joinpath(output_dir, "MLP.jld2") model3 μX σX μY σY L3_ data3 fixedpoints3 pfmc3
savefig(plt_prediction3, joinpath(output_dir, "MLP_prediction.png"))
savefig(plt_loss3, joinpath(output_dir, "MLP_loss.png"))
