using CSV, DataFrames, Plots
default(lw = 2, size = (900, 400), xticks = 0:20:100)


rslt_ = [CSV.read("D:/tlabfqpddjq/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:10]
rslt = rslt_[1]

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 1], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 1], dims = 2), label = "2041")

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 4], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 4], dims = 2), label = "2041")

tensor_rslt = reshape(rslt_[1][:, 4], 100, 2, 17);
plot(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2021")
tensor_rslt = reshape(rslt_[1][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2041 seed 1")

tensor_rslt = reshape(rslt_[2][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2041 seed 2", color = 2)
tensor_rslt = reshape(rslt_[3][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2041 seed 3", color = 2)
tensor_rslt = reshape(rslt_[4][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2041 seed 4", color = 2)
tensor_rslt = reshape(rslt_[4][:, end], 100, 2, 17);
plot!(sum(tensor_rslt[:, 1:2, 14], dims = 2), label = "2041 seed 5", color = 2)
