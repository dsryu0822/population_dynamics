include("../../DataDrivenModel/core/header.jl")
include("../src/0_datacall.jl")

function pyramid_t(matrix; color1_ = [:gold, :orange, :red], color2_ = color1_)
    X = sum_yng.(eachcol(matrix))
    W = sum_prd.(eachcol(matrix))
    Y = sum_old.(eachcol(matrix))
    plt = plot(xticks = [1950, 1970, 2000, 2020], xlims = [1950, 2023], size = [300, 300], ylims = [0, 4e+8], yticks = [0, 4e+8], grid = false)
    plot!(plt, t_, X,         fillrange = 1,     fc = color1_[1], lw = 0)
    plot!(plt, t_, X + W,     fillrange = X,     fc = color1_[2], lw = 0)
    plot!(plt, t_, X + W + Y, fillrange = X + W, fc = color1_[3], lw = 0)
    return plt
end
function pyramid!(plot, matrix, idx_t; color1_ = [:gold, :orange, :red], color2_ = color1_)
    plot!(  0:15, matrix[  1:16, idx_t], color = color1_[1], lw = 0, fillrange = 1, fillcolor = color2_[1])
    plot!( 15:65, matrix[ 16:66, idx_t], color = color1_[2], lw = 0, fillrange = 1, fillcolor = color2_[2])
    plot!(65:100, matrix[66:101, idx_t], color = color1_[3], lw = 0, fillrange = 1, fillcolor = color2_[3])
    return plot
end
POP = CSV.read("$raw/WPP2024_Population1JanuaryBySingleAgeSex_Medium_1950-2023.csv", DataFrame, select = [4, 13, 16, 20])
rename!(POP, [:ISO3, :t, :age, :pop])
dropmissing!(POP)

t_ = 1950:2023
default(); default(dpi = 300, legend = :none, framestyle = :box, left_margin = 5Plots.mm, color = :black, lw = 2)

for i = ["JPN", "BRA", "USA", "PAK"]
    # color1_ = [1, :black, 2]
    # color2_ = [1, :white, 2]
    # color1_ = color2_ = [:dodgerblue, :royalblue, :indigo]
    # color1_ = color2_ = [:green, :darkgoldenrod2, :brown2]
    # color1_ = color2_ = [:darkblue, :darkgoldenrod2, :brown2]
    # color1_ = color2_ = [:gray64, :gray32, :gray0]
    # color1_ = color2_ = [:limegreen, :goldenrod1, :tomato]
    # color1_ = color2_ = [:blue, :limegreen, :orangered3]
    # color1_ = color2_ = [:red, :purple, :blue]
    color1_ = color2_ = [:dodgerblue, :limegreen, :tomato]
    tgt = filter(row -> row.ISO3 == i, POP)
    matrix_tgt = 1000reshape(tgt.pop, 101, 74)

    plt_pyramidt = pyramid_t(matrix_tgt; color1_)
    # png("pyramid1_$i")
    plt_pyramid1 = plot(xlims = [0, 100], xticks = [0, 15, 65, 100], ylims = [1, 7e+6], yticks = [0, 7e+6])
    pyramid!(plt_pyramid1, matrix_tgt, 1; color1_, color2_)
    plt_pyramid2 = plot(xlims = [0, 100], xticks = [0, 15, 65, 100], ylims = [1, 7e+6], yticks = [0, 7e+6])
    pyramid!(plt_pyramid2, matrix_tgt, 74; color1_, color2_)
    # plot(plt_pyramid1, plt_pyramid2, ticks = [], layout = (2, 1), size = [400, 400])
    # png("pyramid2_$i")
    plot(plt_pyramidt, plt_pyramid1, plt_pyramid2, layout = (:, 1), size = [320, 480], yformatter = :scientific)
    png("pyramid_$i")
end

# layout_prmd = @layout [[a; b] c]
# plot(plt_pyramid1, plt_pyramid2, plt_pyramidt, layout = layout_prmd, yticks = :none, size = [800, 400])
asdf  = factory_lorenz(DataFrame, 28)
plot(asdf.x, asdf.y, asdf.z)

plot(
    scatter(
        exp10.(filter(row -> row.t == 2023, data).x),
        exp10.(filter(row -> row.t == 2023, data).y),
        color = :black
    ),
    scatter(
        filter(row -> row.t == 2023, data).x,
        filter(row -> row.t == 2023, data).y,
        color = :black
    ),
color = :black, legend = :none, size = [800, 400], ticks = [], framestyle = :box)
png("temp")