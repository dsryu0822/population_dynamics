using CSV, DataFrames, Plots

name_location = ["Seoul","Busan","Daegu","Incheon","Gwangju","Daejeon","Ulsan","Sejong","Gyeonggi","Gangwon","Chungbuk","Chungnam","Jeonbuk","Jeonnam","Gyeongbuk","Gyeongnam","Jeju"]

n_seed = 10
yend = 2100

rslt_ = [CSV.read("D:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
dead_ = [CSV.read("D:/recent/dead $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
for k in 1:n_seed
    select!(rslt_[k], 1:(yend - 2017))
    select!(dead_[k], 1:(yend - 2017))
end

p1 = plot()
esbl = zeros(100, 2, 17)
for k in 1:n_seed
    # select!(rslt_[k], 1:80)
    rslt = Matrix(rslt_[k][:, 4:end])
    tnsr_rslt = reshape(rslt, 100, 2, 17, :)
    ts_total_pop = vec(sum(tnsr_rslt, dims = 1:3))
    plot!(p1, 2021:yend, ts_total_pop)
    println(sum(ts_total_pop))
end
png("D:/figure/stochastic.png")

for loc ∈ 1:17
    tnsr_rslt = reshape(Matrix(rslt_[1][:, 4:end]), 100, 2, 17, :)
    ts_local_pop = vec(sum(tnsr_rslt[:, :, loc, :], dims = [1,2,4]))
    plot(2021:yend, ts_local_pop, label = name_location[loc], lw = 2, color = :black, legend = :topright)
    extinction50 = 2021 + findfirst(ts_local_pop .< (ts_local_pop[1]/2))
    vline!([extinction50], label  = extinction50)
    png("D:/figure/01 $loc $(name_location[loc]).png")
end
plot(2021:yend, ts_local_pop, label = name_location[loc], lw = 2, color = :black, legend = :topright)

# ----------------------------
