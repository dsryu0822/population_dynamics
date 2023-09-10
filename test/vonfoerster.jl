using Plots, LaTeXStrings; default(color = :black)

begin
    dt = da = 0.1
    a = t = 0:dt:100; length(a)

    μ(a) = (1 - cos(π*(a^2)/10000))/100
    plot(a, μ.(a), label = "Age-specific mortality")

    b(a) = exp((-(a-30)^2)/100)/20
    plot(a, b.(a), label = "Age-specific fertility")

    death = exp.(-μ.(a))
    birth = da*b.(a)
    function update(u)
        _u = circshift(u .* death, 1)
        _u[1] = birth'u
        return _u
    end
end

u0 = [(4:0.002:5)...; (5:-0.01:0.01)...] + 0.1rand(length(a))
# u0 = abs.(sin.(π*a/100)) + abs.(cos.(a/5) + 0.05randn(1001))
u_ = [u0]
for tk in 1:1000
    push!(u_, update(u_[end]))
end

anim = @animate for (tk, u) in enumerate(u_)
    plot(a, u,
    legend = :none, xlabel = L"a", ylabel = "population density", title = L"t = " * "$(rpad(tk/10, 4, '0'))",
    xlims = (0,100), ylims = (0,5))
end
gif(anim)

plot(t, sum.(u_),
xlabel = L"t", ylabel = "total population " * L"\left( \int_\mathbb{R} n(a) da \right) (t)", legend = :none,
leftmargin = 1Plots.cm)

