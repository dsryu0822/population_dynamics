using Plots, LaTeXStrings; default(color = :black)

function update(u, birth)
    _u = circshift(u .* death, 1)
    _u[1] = birth'u
    return _u
end

begin
    dt = da = 0.1
    a = t = 0:dt:100; length(a)

    μ(a) = (1 - cos(π*(a^2)/10000))/10
    plot(a, μ.(a), label = "Age-specific mortality")
    png("mortality")

    b(a) = exp((-(a-30)^2)/100)/20
    plot(a, b.(a), label = "Age-specific fertility")
    png("fertility")

    death = exp.(-da*μ.(a))
    birth = da*b.(a)
end

u0 = [(4:0.002:5)...; (5:-0.01:0.01)...] + 0.1rand(length(a))
# u0 = abs.(sin.(π*a/100)) + abs.(cos.(a/5) + 0.05randn(1001))
u_ = [u0]
for tk in 1:1000
    push!(u_, update(u_[end], birth))
end

anim = @animate for (tk, u) in enumerate(u_)
    plot(a, u,
    legend = :none, xlabel = L"a", ylabel = "population density", title = L"t = " * "$(rpad(tk/10, 4, '0'))",
    xlims = (0,100), ylims = (0,5))
end
mp4(anim, "vanilla.mp4")

plot(t, sum.(u_),
xlabel = L"t", ylabel = "total population " * L"\left( \int_\mathbb{R} n(a) da \right) (t)", legend = :none,
leftmargin = 1Plots.cm)
png("vanilla")

# ---

u0 = [(4:0.002:5)...; (5:-0.01:0.01)...] + 0.1rand(length(a))
# u0 = abs.(sin.(π*a/100)) + abs.(cos.(a/5) + 0.05randn(1001))
u_ = [u0]
for tk in 1:1000
    push!(u_, update(u_[end], 
    birth * 2(sum(u_[end][200:650]) / (sum(u_[end][1:199]) + sum(u_[end][651:end])))/3))
end

anim = @animate for (tk, u) in enumerate(u_)
    plot(a, u,
    legend = :none, xlabel = L"a", ylabel = "population density", title = L"t = " * "$(rpad(tk/10, 4, '0'))",
    xlims = (0,100), ylims = (0,5))
end
gif(anim)

mp4(anim, "toy.mp4")