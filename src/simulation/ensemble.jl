using CSV, DataFrames
using Dates
using Random, Distributions
using Base.Threads

const years = ["y2012","y2013","y2014","y2015","y2016","y2017","y2018","y2019","y2020","y2021"]
const female_ratio = 20 / 41

include("simulation.jl")

POPULATION = CSV.read("data/KOSIS/population.csv", DataFrame)
rename!(POPULATION, ["location", "gender", "age", years...])
filter!(:age => x -> (x != "100세 이상"), POPULATION)
for year ∈ years
    POPULATION[:, year] = trunc.(Int64, POPULATION[:, year])
    POPULATION[!, year] = convert.(Int64, POPULATION[:, year])
end
POPULATION = POPULATION[:, [1,2,3,end]]
const tensor_population = reshape(POPULATION.y2021, 100, 17, 2)

MORTALITY = CSV.read("data/KOSIS/mortality.csv", DataFrame)
rename!(MORTALITY, ["location", "gender", "age", years...])
MORTALITY = MORTALITY[:, [1,2,3,end]]

const tensor_mortality = reshape(MORTALITY.y2021, 100, 17, 2) ./ tensor_population
# tensor_mortality[a+1,-loc,gen+1]

MOBLITY = CSV.read("data/KOSIS/mobility.csv", DataFrame)
rename!(MOBLITY, ["from", "to", "gender", "age", years...])
MOBLITY = MOBLITY[:, [1,2,3,4,end]]
const tensor_mobility = reshape(MOBLITY.y2021, 17, 2, 17, 17)
# tensor_mobility[a5, gen+1, to, from]

FERTILITY = CSV.read("data/KOSIS/fertility.csv", DataFrame)
FERTILITY = FERTILITY[Not(1), [1,((end-6):end)...]]
rename!(FERTILITY, ["location", ("a" .* string.(15:5:45))...])
const tensor_fertility = parse.(Float64, Matrix(FERTILITY[:, Not(1)]))

Threads.@threads for seed ∈ 1:100
    simulation(POPULATION, seed)
end