using CairoMakie, LaTeXStrings

set_theme!(color = :blue)
set_theme!(colors = :darktest)
set_theme!(colormap = :Hiroshige) # https://docs.makie.org/stable/api/index.html#Makie
set_theme!(colormap = :darktest)
set_theme!(colormaps = :darktest)
set_theme!(palette = :darktest)
set_theme!(palette = (patchcolor = [:red, :green, :blue, :yellow, :orange, :pink],)) # https://docs.makie.org/stable/documentation/theming/index.html#palettes
set_theme!(palette = (; patchcolor = cgrad(:Egypt, alpha=0.65))) # https://docs.makie.org/stable/api/index.html#Makie

f = Figure()
ax = Axis(f[1,1])
for _ in 1:10
    lines!(ax, rand(5))
end
f


x = -2pi:0.1:2pi
approx = fill(0.0, length(x))
set_theme!(palette = (; patchcolor = cgrad(:Egypt, alpha=0.65)))
fig, axis, lineplot = lines(x, sin.(x); label = L"sin(x)", linewidth = 3, color = :black,
    axis = (; title = "Polynomial approximation of sin(x)",
        xgridstyle = :dash, ygridstyle = :dash,
        xticksize = 10, yticksize = 10, xtickalign = 1, ytickalign = 1,
        xticks = (-π:π/2:π, ["π", "-π/2", "0", "π/2", "π"])
    ))