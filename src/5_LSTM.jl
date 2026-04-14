using Flux
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

X_ = []
Y_ = []
for gdf = groupby(data, :ISO3)
    df = gdf[:, [:x, :y, :w]]

    tnsr = Matrix(df)
    x_ = [tnsr[t0 .+ (0:9), :] for t0 in 1:50]
    y_ = [tnsr[[t0 + 10], :] for t0 in 1:50]
    push!(X_, cat(x_..., dims = 3))
    push!(Y_, cat(y_..., dims = 3))
end
X = cat(X_..., dims = 3) .|> Float32 # |> gpu
Y = cat(Y_..., dims = 3) .|> Float32 # |> gpu

loader = Flux.DataLoader((X, Y), batchsize = 128, shuffle=true);
Loss(f, x, y) = Flux.mse(f(x), y)

model = Chain(
    LSTM(size(X, 1) => 128),
    Dense(128 => size(Y, 1))
) # |> gpu
optimizer = Flux.setup(Flux.Adam(), model); # will store optimiser momentum, etc.

losses = [Inf]
model_ = [deepcopy(model)]
@showprogress for epoch in 1:100
    Flux.train!(Loss, model, loader, optimizer)
    loss = Loss(model, X, Y)
    # println("Epoch $epoch: Loss = $loss")
    # print(".")
    if loss < minimum(losses)
        push!(model_, deepcopy(model))
        print("New best model found with loss: $(Flux.mse(model_[end](X), Y))\n")
    end
    push!(losses, loss)  # logging, outside gradient context
end
best_model = deepcopy(model_[end])
# JLD2.save("LSTM.jld2", "best_model", best_model)

function prdt_lead(inputX, lead = 1)
    outputY = deepcopy(inputX)
    for ld in 1:lead
        outputY = [outputY; best_model(outputY[ld:end, :, :])]
    end
    return outputY[(size(inputX, 1)+1):end, :, :]
end

t_ = 1950:2022
_t_ = 2008:2022

result = DataFrame(x1 = [], y1 = [], w1 = [], x2 = [], y2 = [], w2 = [], x3 = [], y3 = [], w3 = [])
for k in 1:161
    _X = X[:,:,50k]
    ic = [_X[end, :]...]
    df = groupby(data, :ISO3)[k]
    lstmx, lstmy, lstmw = eachcol([ic'; prdt_lead(_X, 14)[:,:]])
    sindyx, sindyy, sindyw = eachcol(solve(f, ic, 0:0.01:14)[1:100:end, :])

    px = plot(t_, df.x, color = :black, label = "actual", ylabel = "x")
    plot!(px, _t_, lstmx, color = :red, label = "LSTM")
    plot!(px, _t_, sindyx, color = :blue, label = "SINDy")
    py = plot(t_, df.y, color = :black, label = "actual", ylabel = "y")
    plot!(py, _t_, lstmy, color = :red, label = "LSTM")
    plot!(py, _t_, sindyy, color = :blue, label = "SINDy")
    pw = plot(t_, df.w, color = :black, label = "actual", ylabel = "w", xlabel = "t")
    plot!(pw, _t_, lstmw, color = :red, label = "LSTM")
    plot!(pw, _t_, sindyw, color = :blue, label = "SINDy")

    plot(px, pw, py, layout = (3, 1), size = [600, 600], xticks = [1950, 2008, 2022], plot_title = "$(k)_$(only(unique(df.ISO3)))", dpi = 180, legend = :topleft)
    push!(result, last.([df.x, df.y, df.w, lstmx, lstmy, lstmw, sindyx, sindyy, sindyw]))
    png("G:/population/LSTMvsSINDy/$(k)_$(only(unique(df.ISO3))).png")
end


pxlstm  = scatter(result.x1, result.x2, color = 1, msw = 0, label = "LSTM", xlabel = "actual x(2022)", ylabel = "LSTM x(2022)", title = "R2 = $(round(100Rsq(result.x1, result.x2), digits = 2))%")
pylstm  = scatter(result.y1, result.y2, color = 2, msw = 0, label = "LSTM", xlabel = "actual y(2022)", ylabel = "LSTM y(2022)", title = "R2 = $(round(100Rsq(result.y1, result.y2), digits = 2))%")
pwlstm  = scatter(result.w1, result.w2, color = 3, msw = 0, label = "LSTM", xlabel = "actual w(2022)", ylabel = "LSTM w(2022)", title = "R2 = $(round(100Rsq(result.w1, result.w2), digits = 2))%")
pxsindy = scatter(result.x1, result.x3, color = 1, msw = 0, label = "SINDy", xlabel = "actual x(2022)", ylabel = "SINDy x(2022)", title = "R2 = $(round(100Rsq(result.x1, result.x3), digits = 2))%")
pysindy = scatter(result.y1, result.y3, color = 2, msw = 0, label = "SINDy", xlabel = "actual y(2022)", ylabel = "SINDy y(2022)", title = "R2 = $(round(100Rsq(result.y1, result.y3), digits = 2))%")
pwsindy = scatter(result.w1, result.w3, color = 3, msw = 0, label = "SINDy", xlabel = "actual w(2022)", ylabel = "SINDy w(2022)", title = "R2 = $(round(100Rsq(result.w1, result.w3), digits = 2))%")
plot(pxlstm, pwlstm, pylstm, pxsindy, pwsindy, pysindy, layout = (2, 3), size = [800, 600], legend = :none, dpi = 180)
png("G:/population/LSTMvsSINDy/_summary.png")
# (1:161)[sortperm(abs.(result.w1 - result.w3))]

