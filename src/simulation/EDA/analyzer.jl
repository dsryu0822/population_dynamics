using CSV, DataFrames, Plots
default(lw = 2, size = (900, 400), xticks = 0:30:100)

rslt_ = [CSV.read("D:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:10]

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 1], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 1], dims = 2), label = "2051")

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 3], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 3], dims = 2), label = "2051")

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2051 seed 1")

tensor_rslt = reshape(rslt_[2][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2051 seed 2", color = 2)
tensor_rslt = reshape(rslt_[3][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2051 seed 3", color = 2)
tensor_rslt = reshape(rslt_[4][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2051 seed 4", color = 2)
tensor_rslt = reshape(rslt_[4][:, 4 + 30], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2051 seed 5", color = 2)

p1 = plot()
for k ∈ 1:10
    tensor_rslt = reshape(rslt_[k][:, end], 100, 2, 17);
    plot!(p1, 2021:2100, vec(sum(Matrix(rslt_[k][:,4:end]), dims = 1)),
     xticks = 2020:20:2100, label = "seed $k", title = "total population")
end
p1


tensor_rslt = reshape(Matrix(rslt_[1][:, 4:end]), 100, 2, 17, 80)

수도권 = vec(sum(tensor_rslt[:, :, [1,4,9], :], dims = 1:3));
대구 = vec(sum(tensor_rslt[:, :, [3], :], dims = 1:3));
경북 = vec(sum(tensor_rslt[:, :, [15], :], dims = 1:3));
지방 = vec(sum(tensor_rslt[:, :, Not([1,4,9]), :], dims = 1:3));

plot(
  plot(2021:2100, 대구, xticks = 2020:20:2100, )
, plot(2021:2100, 경북, xticks = 2020:20:2100, )
)
plot(2021:2100, 수도권 ./ (수도권 + 지방), xticks = 2020:20:2100)

생산인구 = vec(sum(tensor_rslt[(15:64), :, :, :], dims = 1:3))
부양인구 = vec(sum(tensor_rslt[Not(15:64), :, :, :], dims = 1:3))
plot(2021:2100, 부양인구 ./ (생산인구 + 부양인구), xticks = 2020:20:2100)

# ----------------

rslt_ = [CSV.read("D:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:10]

p2 = plot()
for k ∈ 1:10
    tensor_rslt = reshape(rslt_[k][:, end], 100, 2, 17);
    plot!(p2, 2021:2050, vec(sum(Matrix(rslt_[k][:,4:end]), dims = 1)),
     xticks = 2020:10:2050, label = "seed $k", title = "total population")
end
p2