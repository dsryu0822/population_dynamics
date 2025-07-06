using Plots
include("../../setup/Dynamics.jl")

mine(x, α, λ₁, λ₂) = α*exp(-λ₁*x) + (1-α)*λ₂*log(101 - x)

plot(mine.(0:100, 0, 1, 1))
plot(mine.(0:100, 0.5, 0.1, 0.5))
plot(mine.(0:100, 1, 1, 1))

balance(αλ) = o0*sum(mine.(1:15, αλ...)) - y0*sum(mine.(65:100, αλ...))

y0, o0 = 0.198, 0.221
αλ0 = findfixed(balance, [0.0, 0.1, 0.5], h = .0001)
pop0 = mine.(1:100, αλ0...)
plot(pop0, xticks = [0, 15, 65, 100], ylims = [0, Inf], color = :black, lw = 2)

y0, o0 = 0.19185957,    0.005762051
αλ0 = findfixed(balance, [1.0, 0.1, 0.5], h = .0001)
pop0 = mine.(1:100, αλ0...)
plot(pop0, xticks = [0, 15, 65, 100], ylims = [0, Inf], color = :black, lw = 2)

y0, o0 = 0.12348538,	0.261549377
αλ0 = findfixed(balance, [1., 1., 1.], h = .00001)
pop0 = mine.(1:100, αλ0...)
plot(pop0, xticks = [0, 15, 65, 100], ylims = [0, Inf], color = :black, lw = 2)

y0, o0 = 0.30635788,	0.121840009
αλ0 = findfixed(balance, [0.0, 0.1, 0.1], h = .0001)
pop0 = mine.(1:100, αλ0...)
plot(pop0, xticks = [0, 15, 65, 100], ylims = [0, Inf], color = :black, lw = 2)

# ----

mortality = (pop0 ./ circshift(pop0, 1))[2:end]
plot(1 .- mortality, yscale = :log10, color = :red, lw = 2)

fertiliry(p) = exp.(-(((1:100) .- p[1])/p[2]) .^ 2)/p[3]
birth(p) = first(pop0) - sum(pop0 .* fertiliry(p))

p0 = findfixed(birth, [30., 1, 10], h = .01)
plot(fertiliry(p0), xticks = [0, 15, 20, 25, 30, 35, 40, 45, 65, 100], color = :green, lw = 2)

Leslie = zeros(100, 100)
Leslie[1,:] = fertiliry(p0)
for i in 2:100
    Leslie[i,i-1] = mortality[i-1]
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
