using CSV, DataFrames
using LinearAlgebra

using Plots, LaTeXStrings

default(fontfamily = "Computer Modern")
# plot(rand(10), title = "abcdefu")
red = colorant"#C00000"
orange = colorant"#ED7D31"
yellow = colorant"#FFC000"

function ts_sum(df)
    return sum.(eachcol(select(df, Not([:location, :gender, :age]))))
end

function marginal(df, col)
    return combine(groupby(df, col), ["y$t" => sum => "y$t" for t = 2012:2070])[:, Not(col)]
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
yend = 2070

rslt_ = [CSV.read("G:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
dead_ = [CSV.read("G:/recent/dead $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
mgrn_ = [CSV.read("G:/recent/mgrn $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
for k in 1:n_seed
    select!(rslt_[k], 1:(yend - 2017 + 9))
    select!(dead_[k], 1:(yend - 2017 + 9))
    select!(mgrn_[k], 1:(yend - 2017 + 9 + 1))
    rslt_[k].age = parse.(Int, replace.(rslt_[k].age, '세' => ""))
    dead_[k].age = parse.(Int, replace.(dead_[k].age, '세' => ""))
    mgrn_[k].age = parse.(Int, first.(mgrn_[k].age, 2))
end

rslt = rslt_[1]
dead = dead_[1]
mgrn = mgrn_[1]

# # 시작 시뮬레이션 타당성 검토
# k = 1
# t = 4
# MOBILITY = CSV.read("data/KOSIS/mobility.csv", DataFrame)
# rename!(MOBILITY, ["from", "to", "gender", "age", years...])
# MOBILITY = MOBILITY[:, [1,2,3,4,end]]

# migration_matrix = reshape(sum(reshape(MOBILITY.y2012, 17,2,17,17), dims = 1:2), 17, 17)
# for d ∈ 1:17 migration_matrix[d,d] = 0 end
# # 세로 전출, 가로 전출 ex) 서울 → 부산: 13078 = M(2,1)

# heatmap(migration_matrix, size = (800, 800), ticks = [4,9,17], xlabel = "from", ylabel = "to")

# 이동수지 = Int64[]
# for d in 1:17
#     push!(이동수지, sum(migration_matrix[d,:]) - sum(migration_matrix[:,d]))
# end
# 이동수지


# 인구수지 = vec(sum(reshape(rslt_[k][:, t+1], 200, :), dims = 1)) - 
#  vec(sum(reshape(rslt_[k][:, t], 200, :), dims = 1))

# 신생아 = vec(sum(reshape(rslt_[k][:, t+1], 100, 2, :)[1,:,:], dims = 1))
# 사망자 = vec(sum(reshape(dead_[k][:, t], 200, :), dims = 1))
# 이동수지 + 신생아 - 사망자
# plot(인구수지, label = "A")
# plot!(이동수지 + 신생아 - 사망자, label = "B")
# # 끝 시뮬레이션 타당성 검토