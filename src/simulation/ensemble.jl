using Dates
println(Dates.now())
println(Threads.nthreads(), " workers are found")

# using CUDA
using CSV, DataFrames
using Random, Distributions
using Base.Threads

const years = lpad.(2000:2021, 5, 'y')
const female_ratio = 100 / (105 + 100) # 성비는 남:여 = 105:100

const ε = 0 # 이거 근데 어차피 폐기될거고 각 텐서별로 분산 구해서 쓰게될듯 230107 기준으로 0.2 꽤 나쁘지 않았음

include("simulation.jl")

POPULATION = CSV.read("data/KOSIS/population.csv", DataFrame)
POPULATION = ifelse.(ismissing.(POPULATION), 1, POPULATION)
rename!(POPULATION, ["location", "gender", "age", years...])
for year ∈ years
    POPULATION[:, year] = trunc.(Int64, POPULATION[:, year])
    POPULATION[!, year] = convert.(Int64, POPULATION[:, year])
end

MORTALITY = CSV.read("data/KOSIS/mortality.csv", DataFrame)
MORTALITY = ifelse.(ismissing.(MORTALITY), 1, MORTALITY)
rename!(MORTALITY, ["location", "gender", "age", years...])

MOBILITY = CSV.read("data/KOSIS/mobility.csv", DataFrame)
MOBILITY = ifelse.(ismissing.(MOBILITY), 0, MOBILITY)
rename!(MOBILITY, ["from", "to", "gender", "age", years...])

FERTILITY = CSV.read("data/KOSIS/fertility.csv", DataFrame)
FERTILITY = ifelse.(ismissing.(FERTILITY), "0", FERTILITY)
# argyear = findfirst(names(FERTILITY) .== string(tbgn))
# FERTILITY = FERTILITY[Not(1), [1,(argyear:(argyear + 6))...]]
# rename!(FERTILITY, ["location", ("a" .* string.(15:5:45))...])


# println("Simulation start at $(Dates.now()): ")
# @threads for tbgn ∈ 2000:2020
#     simulation(tbgn, tbgn, 2021)
# end

# seed ∈ 2012:2020, simulation(seed, seed, 2100): validation
# seed ∈ [1000]   , simulation(seed, 2021, 2100): 2030년까지 경기도로 유입 인구 2배
@threads for seed ∈ 1001:1017
    simulation(seed, 2021, 2100)
end