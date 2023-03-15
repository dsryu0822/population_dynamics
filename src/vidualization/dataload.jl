using CSV, DataFrames
using LinearAlgebra

function ts_sum(df)
    return sum.(eachcol(select(df, Not(1:3))))
end

function marginal(df, col)
    return combine(groupby(df, col), ["y$t" => sum => "y$t" for t = ybgn:yend])[:, Not(col)]
end

is_capital = (loc -> loc ∈ ["서울특별시", "인천광역시", "경기도"])

function ts_age(df)
    age_l = ts_sum(df[df.age .< 15,:])
    age_m = ts_sum(df[15 .≤ df.age .< 65,:])
    age_h = ts_sum(df[65 .≤ df.age,:])
    return [age_l age_m age_h]
end


name_location = ["Seoul","Busan","Daegu","Incheon","Gwangju","Daejeon","Ulsan","Sejong","Gyeonggi","Gangwon","Chungbuk","Chungnam","Jeonbuk","Jeonnam","Gyeongbuk","Gyeongnam","Jeju"]

n_seed = 10
ybgn = 2012
yend = 2050


rslt_ = [CSV.read("G:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
dead_ = [CSV.read("G:/recent/dead $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
mgrn_ = [CSV.read("G:/recent/mgrn $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
vldn_ = Dict(2000:2020 .=> [CSV.read("G:/recent/vldn/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 2000:2020])
for k in 1:n_seed
    select!(rslt_[k], ["location", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    select!(dead_[k], ["location", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    select!(mgrn_[k], ["from", "to", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    rslt_[k].age = parse.(Int, replace.(rslt_[k].age, '세' => ""))
    dead_[k].age = parse.(Int, replace.(dead_[k].age, '세' => ""))
    mgrn_[k].age = parse.(Int, first.(mgrn_[k].age, 2))
end

try mkdir("G:/figure") catch IOError println("G:/figure already exists") end
try mkdir("G:/figure/subfigure") catch IOError println("G:/figure/subfigure already exists") end
try mkdir("G:/figure/snapshot") catch IOError println("G:/figure/snapshot already exists") end

rslt = rslt_[1]
dead = dead_[1]
mgrn = mgrn_[1]
real = CSV.read("data/KOSIS/population.csv", DataFrame)
vldn = CSV.read("G:/recent/vldn/rslt 2012.csv", DataFrame)
이름_지역 = unique(rslt.location)

df_age = marginal(rslt, :age)

println("Data load done")