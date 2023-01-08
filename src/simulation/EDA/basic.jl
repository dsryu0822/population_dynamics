using CSV, DataFrames, Plots
using LinearAlgebra

function ts_sum(df)
    return sum.(eachcol(select(df, Not([:location, :gender, :age]))))
end

name_location = ["Seoul","Busan","Daegu","Incheon","Gwangju","Daejeon","Ulsan","Sejong","Gyeonggi","Gangwon","Chungbuk","Chungnam","Jeonbuk","Jeonnam","Gyeongbuk","Gyeongnam","Jeju"]

n_seed = 10
yend = 2050

rslt_ = [CSV.read("D:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
dead_ = [CSV.read("D:/recent/dead $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
mgrn_ = [CSV.read("D:/recent/mgrn $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
for k in 1:n_seed
    select!(rslt_[k], 1:(yend - 2017))
    select!(dead_[k], 1:(yend - 2017))
    select!(mgrn_[k], 1:(yend - 2017 + 1))
    rslt_[k].age = parse.(Int, replace.(rslt_[k].age, '세' => ""))
    dead_[k].age = parse.(Int, replace.(dead_[k].age, '세' => ""))
    mgrn_[k].age = parse.(Int, first.(mgrn_[k].age, 2))
end

