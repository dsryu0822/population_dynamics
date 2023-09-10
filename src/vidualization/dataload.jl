using CSV, DataFrames
using LinearAlgebra

n_seed = 10
ybgn = 2012
yend = 2100

interested = lpad.(ybgn:yend, 5, "y")
이름_지역 = ["서울특별시","부산광역시","대구광역시","인천광역시","광주광역시","대전광역시","울산광역시","세종특별자치시","경기도","강원도","충청북도","충청남도","전라북도","전라남도","경상북도","경상남도","제주특별자치도"]
name_location = ["Seoul","Busan","Daegu","Incheon","Gwangju","Daejeon","Ulsan","Sejong","Gyeonggi","Gangwon","Chungbuk","Chungnam","Jeonbuk","Jeonnam","Gyeongbuk","Gyeongnam","Jeju"]

a = 0:99

function as(type, df)
    temp = select(df, interested)
    if type == DataFrame
        return temp
    elseif type == Matrix
        return Matrix(temp)
    end
end

function marginby(df, col)
    gdf = groupby(df, col)
        
    result = deepcopy(gdf[1])
    result[:, interested] .= 0
    for df in gdf
        result[:, interested] .+= as(DataFrame, df)
    end

    if eltype(df[:, col]) <: Number
        result[:, col] .= -1
    elseif eltype(df[:, col]) <: AbstractString
        result[:, col] .= "-"
    end

    return result
end

function ts_sum(df)
    return sum.(eachcol(as(DataFrame, df)))
end

is_capital = (loc -> loc ∈ ["서울특별시", "인천광역시", "경기도"])


rslt_ = [CSV.read("G:/recent/rslt $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
dead_ = [CSV.read("G:/recent/dead $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
mgrn_ = [CSV.read("G:/recent/mgrn $(lpad(seed, 4, '0')).csv", DataFrame) for seed = 1:n_seed]
for k in 1:n_seed
    select!(rslt_[k], ["location", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    select!(dead_[k], ["location", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    select!(mgrn_[k], ["from", "to", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    rslt_[k].age = parse.(Int, replace.(rslt_[k].age, '세' => ""))
    dead_[k].age = parse.(Int, replace.(dead_[k].age, '세' => ""))
    mgrn_[k].age = parse.(Int, first.(mgrn_[k].age, 2))
end

vldn_ = Dict()
for seed = ybgn:2020
    temp = CSV.read("G:/recent/vldn/rslt $(lpad(seed, 4, '0')).csv", DataFrame)
    select!(temp, ["location", "gender", "age", lpad.(ybgn:yend, 5, 'y')...])
    temp.age = parse.(Int, replace.(temp.age, '세' => ""))
    push!(vldn_, seed => temp)
end

try mkdir("G:/figure") catch IOError println("G:/figure already exists") end
try mkdir("G:/figure/subfigure") catch IOError println("G:/figure/subfigure already exists") end
try mkdir("G:/figure/snapshot") catch IOError println("G:/figure/snapshot already exists") end

real = CSV.read("data/KOSIS/population.csv", DataFrame)
rslt = rslt_[1]
dead = dead_[1]
mgrn = mgrn_[1]

# function mg_sum(df, col)
#     return combine(groupby(df, col), ["y$t" => sum => "y$t" for t = ybgn:yend])[:, Not(col)]
# end

# df_age = mg_sum(rslt, :age)

# function ts_age(df)
#     age_l = ts_sum(df[df.age .< 15,:])
#     age_m = ts_sum(df[15 .≤ df.age .< 65,:])
#     age_h = ts_sum(df[65 .≤ df.age,:])
#     return [age_l age_m age_h]
# end

println("Data load done")