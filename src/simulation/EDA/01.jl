using CSV, DataFrames, Plots

n_seed = 10
rslt_ = [CSV.read("D:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]

p1 = plot()
esbl = zeros(100, 2, 17)
for k in 1:n_seed
    select!(rslt_[k], 1:80)
    rslt = Matrix(rslt_[k][:, 4:end])
    tnsr_rslt = reshape(rslt[:, 4:end], 100, 2, 17, :)
    ts_total_pop = vec(sum(tnsr_rslt, dims = 1:3))
    yend = 2021 + (length(ts_total_pop)-1)
    plot!(p1, 2021:yend, ts_total_pop)
    println(sum(ts_total_pop))
end

p1