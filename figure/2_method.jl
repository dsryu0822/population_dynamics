

data_ = Dict(unique(data.ISO3) .=> collect(groupby(data, :ISO3)) .|> dropmissing)

# tgt = filter(row -> row.ISO3 == "CHN", POP)
# matrix_tgt = reshape(tgt.pop, 101, 74)
# plt_pyramid1 = plot(xlims = [0, 100], xticks = [], ylims = [1, 1000maximum(matrix_tgt)])
# pyramid!(plt_pyramid1, matrix_tgt, 1)
# plt_pyramid2 = plot(xlims = [0, 100], xticks = [], ylims = [1, 1000maximum(matrix_tgt)])
# pyramid!(plt_pyramid2, matrix_tgt, 74)
# plot(plt_pyramid1, plt_pyramid2, ticks = [], layout = (2, 1), size = [400, 400])
# png("pyramid2_CHN")

# plot(data_["CHN"].x, data_["CHN"].y, size = [400, 400])
# plot(data_["CHN"].x, data_["CHN"].y, size = [400, 400], aspect_ratio = 1, ticks = [], xlims = 4.9 .+ [0, 1], ylims = 4.35 .+ [0, 1])
# scatter!(data_["CHN"].x[[end]], data_["CHN"].y[[end]], ms = 10)
# png("data_CHN")

color_ = [:dodgerblue, :limegreen, :tomato]

tgt = filter(row -> row.ISO3 == "IND", POP)
matrix_tgt = reshape(tgt.pop, 101, 74)
plt_pyramid1 = plot(xlims = [0, 100], xticks = [], ylims = [1, maximum(matrix_tgt)], yticks = [])
pyramid!(plt_pyramid1, matrix_tgt, 1, color1_ = color_)
png("temp1")
plt_pyramid2 = plot(xlims = [0, 100], xticks = [], ylims = [1, maximum(matrix_tgt)], yticks = [])
pyramid!(plt_pyramid2, matrix_tgt, 74, color1_ = color_)
png("temp2")

plt_pyramid_ = []
for t = [1950, 1970, 1990, 2010, 2023]
    plt_pyramid = plot(xlims = [0, 100], xticks = [], ylims = [1, maximum(matrix_tgt)], yticks = [], size = [300, 100])
    pyramid!(plt_pyramid, matrix_tgt, t - 1949, color1_ = color_)
    push!(plt_pyramid_, plt_pyramid)
    # png("pyramid_$t")
end
plot(reverse(plt_pyramid_)..., layout = (:, 1), size = [200, 400], margin = 0mm)
png("pyramid_t")

# plot(data_["IND"].x, data_["IND"].y, size = [400, 400])
plt_india = plot(data_["IND"].x, data_["IND"].y, size = [400, 400], aspect_ratio = 1, ticks = [], xlims = 5.5 .+ [-.6, .6], ylims = 4.5 .+ [-.6, .6], framestyle = :none, arrow =  Plots.arrow(:open, :head, 1, 1))
# scatter!(data_["IND"].x[[end]], data_["IND"].y[[end]], ms = 10)
png("data_IND")

trajx = Matrix(CSV.read("data/trajx2.csv", DataFrame))'
trajy = Matrix(CSV.read("data/trajy2.csv", DataFrame))'
plt_vf = plot(lims = [1, 6])
for k in axes(trajx, 1)
    plot!(plt_vf, trajx[k, :], trajy[k, :], color = :black, label = :none, arrow = Plots.arrow(:open, :head, .6, .6))
end
plot(plt_vf, data_["IND"].x, data_["IND"].y, size = [400, 400], aspect_ratio = 1, ticks = [], xlims = 5.5 .+ [-.6, .6], ylims = 4.5 .+ [-.6, .7], framestyle = :none, arrow =  Plots.arrow(:open, :head, 1, 1))
png("temp")