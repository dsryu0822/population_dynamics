using Distributions, Plots

→ = 0:0.1:5
← = reverse(→)
a1 = @animate for α ∈ [→; ←; -→; -←]
    SN = SkewNormal(0, 1, α)
    x = -3:0.01:3
    plot(x, pdf.(SN, x),
    ylim = (0, 0.8),
    lw = 2, color = :black,
    legend = :none, title = "pdf of Skew Normal Distribution with α = $(lpad(α, 4))")
end
gif(a1, string(@__DIR__) * "/skew_normal.gif")