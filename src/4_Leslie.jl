circshift(1:4, 1)
pop2022 = findrow(POP, :ISO3 => "KOR", :t => 2022).pop
pop2023 = findrow(POP, :ISO3 => "KOR", :t => 2023).pop
s_ = circshift(pop2023, -1) ./ pop2022
# plot(s_[1:end-2], xticks = 0:20:100)

x_ = Dict([2022 => pop2022[1:end-1]])

gauss(x, μ, k) = k*exp.(-.5((x .- μ)/5).^2)
plot(gauss(1:100, 10, 0.1))
residual(mu, k) = pop2023[1] - x_[2022]'gauss(1:100, mu, k)

μ_ = 20:.1:40
k_ = 0:0.001:1
residuals = reshape([(μ, k, residual(μ, k)) for μ in μ_ for k in k_], length(k_), length(μ_))
argmini, argminj = argmin(abs.(last.(residuals))).I
residual(μ_[argminj], k_[argmini])
f_ = gauss(1:100, μ_[argminj], k_[argmini])


mL = [f_'; [diagm(s_[1:end-2]) zeros(99, 1)]]

[push!(x_, t+1 => mL*x_[t]) for t in 2022:2100]
plot(
    plot(x_[2023]),
    plot(x_[2030]),
    plot(x_[2050]),
    plot(x_[2100]),
    layout = (4, 1), size = [300, 600], ylims = [0, 1000], yticks = [], xticks = [0, 15, 65, 100], xlims = [0, 100]
)