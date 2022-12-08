using CSV, DataFrames, Plots

kas = CSV.read("./data/PopulationPyramid/korea_age_structure.csv", DataFrame)

S = []

p1 = plot(xlabel = "age a", ylabel = "survival rate s(a)", legend = :bottomleft)
for t in 1:15
Xt1 = vec(Matrix(kas[[t], Not(end)]))
Xt2 = vec(Matrix(kas[[t+1], Not(1)]))
s_ = Xt2 ./ Xt1
push!(S, s_)
plot!(p1, 5:5:100, s_, alpha = t/15, label = "$(1945 + 5t)")
end
p1

b = rand(21)
xt1 = vec(Matrix(kas[[1], :]))
xt2 = vec(Matrix(kas[[2], :]))

abs2(b'xt1 - xt2[1])