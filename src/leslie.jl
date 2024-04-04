λ = -0.01

y0 = 0.198
o0 = 0.221

balance(λ) = o0*sum(exp.(λ*(1:15))) - y0*sum(exp.(λ*(65:100)))

searchrange = -(0.01:0.000001:0.011)
λ0 = searchrange[argmin(abs.(balance.(searchrange)))]

pop0 = exp.(λ0*(1:100))

using Plots

include("../../setup/Dynamics.jl")
plot(pop0, xticks = [0, 15, 65, 100], ylims = [0, Inf], color = :black, lw = 2)

fertiliry(p) = exp.(-(((1:100) .- p[1])/p[2]) .^ 2)/p[3]
birth(p) = exp(λ0) - sum(pop0 .* fertiliry(p))

birth([15, 1, 1.0])

p0 = findfixed(birth, [30., 1, 10], h = .1)
plot(fertiliry(p0), xticks = [0, 15, 20, 25, 30, 35, 40, 45, 65, 100], color = :green, lw = 2)

Leslie = zeros(100, 100)
Leslie[1,:] = fertiliry(p0)
for i in 2:100
    Leslie[i,i-1] = exp(λ0)
end
using SparseArrays
sparse(Leslie)

pop_ = [pop0]
for _ in 2:1000
    push!(pop_, Leslie * pop_[end])
end
p1 = plot(xticks = [0, 15, 65, 100], ylims = [0, Inf], legend = :best)
p1 = plot(p1, pop_[1], lw = 2, label = "t=1")
p1 = plot(p1, pop_[end], lw = 2, style = :dash, label = "t=1000")