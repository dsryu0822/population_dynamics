include("dataload.jl")
include("backend_GR.jl")

default()
default(legend = :none, plot_titlelocation = :left, lw = 2, warn_on_unsupported=true)

KLD(p,q) = -sum(p .* log.(q ./ p))
Rsq(y, ŷ) = 1 - sum(abs2, (y .- ŷ)) / sum(abs2, (y .- (sum(y) / length(y))))
# Rsq(y, ŷ, y0) = sum(abs2, (ŷ .- y0)) / sum(abs2, (y .- y0))


vldn = vldn_[2012]


gdf_rslt = groupby(marginby(rslt, :gender), :location)
gdf_vldn = groupby(marginby(vldn, :gender), :location)

plt_compare__ = []
L1__ = []
KL__ = []
R2__ = []
for i in 1:17
    df_rslt = gdf_rslt[i]
    df_vldn = gdf_vldn[i]
    plt_compare_ = []
    L1_ = []
    KL_ = []
    R2_ = []
    for j in 1:9
        # xlabel = L"a", ylabel = L"p(a)", 
        ptemp = plot(xlims = (0, 99), ylims = (0, Inf), grid = false, ticks = false)
        N_a = df_rslt[:, 3 + j]; p_a = N_a ./ sum(N_a)
        N̂_a = df_vldn[:, 3 + j]; p̂_a = N̂_a ./ sum(N̂_a)
        push!(L1_, sum(abs, N_a .- N̂_a))
        push!(KL_, KLD(p_a, p̂_a))
        push!(R2_, Rsq(N_a, N̂_a))

        plot!(ptemp, a, N_a, color = :black)
        plot!(ptemp, a, N̂_a, color = 무지개[i], ls = :dash)
        if j == 9
            xlabel!(ptemp, L"a")
            plot!(ptemp, xticks = [0,50, 99])
        end
        if i == 1 ylabel!(ptemp, "$(j+2011) " *  L"N(a)") end
        if j == 1 title!(ptemp, name_location[i]) end
        push!(plt_compare_, ptemp)
    end
    push!(plt_compare__, plt_compare_)
    push!(L1__, L1_)
    push!(KL__, KL_)
    push!(R2__, R2_)
end
annotate!(plt_compare__[17][1], 99, 13500, text("a", 14, :black))


plt_L1 = plot(titleloc = :right, title = "b", ylabel = L"L_{t}", ylims = (0, 610000), xlims = (2012,2020.5), yticks = [0, 300000, 600000], xformatter = (_...) -> "", yformatter = x -> x / 10^( 5))
plt_KL = plot(titleloc = :right, title = "c", ylabel = L"D_{t}",     ylims = (0, 0.02), xlims = (2012,2020.5), yticks = [0, 0.01, 0.02], yformatter = x -> x / 10^(-2), xlabel = L"t")
plt_R2 = plot(titleloc = :right, title = "d", ylabel = L"R^{2}", ylims = (0, 1),   xlims = (2012,2020.5), yticks = [0, 1])
for i in 1:17
    plot!(plt_L1, 2012:2020, L1__[i], color = 무지개[i])
    plot!(plt_KL, 2012:2020, KL__[i], color = 무지개[i])
    plot!(plt_R2, 2012:2020, R2__[i], color = 무지개[i])
    # scatter!(plt_L1, [2020], [last.(L1__)[i]], mc = 무지개[i], msw = 0, shape = :utriangle, ms = 5)
    scatter!(plt_KL, [2020], [last.(KL__)[i]], mc = 무지개[i], msw = 0, shape = :rect)
    # scatter!(plt_R2, [2020], [last.(R2__)[i]], mc = 무지개[i], msw = 0, shape = :dtriangle, ms = 5)
end
plot!(plt_KL, right_margin = 3mm)

annotate!(plt_L1, 2012, 660000, text(L"\times 10^{5}", :left, 10, :black))
annotate!(plt_KL, 2012, 0.022, text(L"\times 10^{-2}", :left, 10, :black))

annotate!(plt_L1, 2017, 550000, text("Gyounggi", 10, 무지개[9]))
plot!(plt_L1, [2018.5, 2019.5], [550000, 550000], arrow = true, color = 무지개[9])
annotate!(plt_L1, 2014, 350000, text("Seoul", 10, 무지개[1]))
plot!(plt_L1, [2014, 2014], [300000, 150000], arrow = true, color = 무지개[1])

annotate!(plt_KL, 2014, 0.015, text("Sejong", 10, 무지개[8]))
plot!(plt_KL, [2015, 2016.5], [0.015, 0.015], arrow = true, color = 무지개[8])

annotate!(plt_R2, 2014, 0.5, text("Sejong", 10, 무지개[8]))
plot!(plt_R2, [2015, 2016.5], [0.5, 0.5], arrow = true, color = 무지개[8])



grslt = groupby(marginby(rslt, :gender), :location)


whysejong = plot(titleloc = :right, title = "d", legend = :topleft, xlabel = L"t_{0}", xticks = [2012, 2019], xlims = [2011.5, 2019.5], ylims = (0, Inf), yticks = [0, 80000], yformatter = x -> x / 10^(4), ylabel = L"\Delta N" * " of Sejong")
plot!(whysejong, (2012:2020) .- 0.2 , ts_sum(groupby(mgrn, :to)[8]),        lw = 0, bar_width = 0.4, st = :bar, alpha = 0.5, color = :black,  label = "Immigation")
plot!(whysejong, (2012:2020) .+ 0.2 , vec(Matrix(grslt[8][[1], Not(1:3)])), lw = 0, bar_width = 0.4, st = :bar, alpha = 1.0, color = :black, label = "Birth")
annotate!(whysejong, 2011.5, 93000, text(L"\times 10^{4}", :left, 10, :black))

Dvalues = zeros(17, 8)
Rvalues = zeros(17, 8)
for (j, y) in enumerate(2012:2019)
    __vldn_ = deepcopy(marginby(vldn_[y], :gender))
    gvldn = groupby(__vldn_, :location)
    for i in 1:17
        N_a = grslt[i][:, "y2020"]; p_a = N_a ./ sum(N_a)
        N̂_a = gvldn[i][:, "y2020"]; p̂_a = N̂_a ./ sum(N̂_a)
        # Rvalues[i,j] = Rsq(N_a, N̂_a)
        # Rvalues[i,j] = Rsq(p_a, p̂_a)
        Dvalues[i,j] = KLD(p_a, p̂_a)
    end
end

stabilized = plot(titleloc = :right, title = "e", xlabel = L"t_{0}", xticks = [2012, 2019] , xlims = (2011.8, 2019.2), ylabel = L"D_{2020}", yticks = [0, 0.01, 0.02], yformatter = x -> x / 10^(-2), ylims = (0, 0.026))
for (i, lv) in enumerate([eachrow(Dvalues)...])
    plot!(stabilized, 2012:2019, lv, color = 무지개[i], ls = :solid, msw = 0, shape = :rect)
end
annotate!(stabilized, 2011.8, 0.028, text(L"\times 10^{-2}", :left, 10, :black))

annotate!(stabilized, 2015, 0.015, text("Sejong", 10, 무지개[8]))
plot!(stabilized, [2014.4, 2013.8], [0.015, 0.015], arrow = true, color = 무지개[8])
plot!(stabilized, right_margin = 3mm)


plot(
    plot(
        plot(plt_compare__[1][1], left_margin = 8mm), plt_compare__[8][1], plt_compare__[17][1],
        plt_compare__[1][5], plt_compare__[8][5], plt_compare__[17][5],
        plt_compare__[1][9], plt_compare__[8][9], plt_compare__[17][9],
        layout = (3,3)
    ),
    plt_L1,
    plt_KL,
    plot(whysejong, stabilized, layout = @layout [d{0.5w} e]),
    size = (900, 700), layout = @layout [a{0.7w} [b ; c] ; de{0.25h}]
)


png("G:/figure/21_1.png")