using DataDrivenDiffEq
# using ModelingToolkit
using OrdinaryDiffEq
# using DataDrivenSparse
using LinearAlgebra
using DataDrivenDMD
using CSV, DataFrames

cd(@__DIR__); pwd()
data = CSV.read("../data/australia/rainbow/Australia_age_specific_fertility.csv", DataFrame)

# sol = rand(3,10)
ddprob = DataDrivenProblem(Matrix(data)[5:5:35,:])

@parameters t
@variables u(t)[1:7]
Ψ = Basis([u; u[1]^2], u, independent_variable = t)
@time result = solve(ddprob, Ψ, TOTALDMD(), digits = 4)
get_algorithm(result)
println(result)

basis = get_basis(result)
get_parameter_map(basis)