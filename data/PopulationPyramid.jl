using CSV, DataFrames, Plots

cd("data/PopulationPyramid/")
age = 0:5:100

X = []
a1 = @animate for csvfile = readdir()
    year = csvfile[(end-7):(end-4)]
    if csvfile[(end-2):end] != "csv" continue end
    data = CSV.read(csvfile, DataFrame)

    plot(title = "Population Pyramid at " * year)
    plot!(age, data.M, lt = :bar, color = :blue, alpha = .5, label = "male")
    plot!(age, data.F, lt = :bar, color = :red, alpha = .5, label = "female")
    Xₜ = data.M + data.F
    push!(X, Xₜ)
end
gif(a1, "PopulationPyramid.gif", fps = 1)

X = DataFrame(hcat(X...)', :auto)
rename!(X, "age" .* lpad.(string.(age), 2, "0"))
CSV.write("korea_age_structure.csv", X)