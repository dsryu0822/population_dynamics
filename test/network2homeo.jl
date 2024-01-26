using ProgressMeter
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering, :Distances, :NetworkLayout]
@showprogress for package in packages
    @eval using $(package)
end
mm = Plots.mm
cm = Plots.cm
default(); default(legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))
normalize01(v) = (v .- minimum(v)) / (maximum(v) - minimum(v))
getrratio(v) = sum(v[65:end]) / sum(v[[1:15; 65:end]])

function cosine_matrix(v)
    n = length(v)
    cosM = zeros(n, n)
    for i in 1:n, j in 1:n
        if i < j
            cosM[i,j] = cosine_similarity(v[i], v[j])
        end
    end
    return cosM
end

function padmissing(x)
    if x .|> !ismissing |> all
        return x
    elseif x .|> ismissing |> all
        @error "padmissing: x is all missing"
    end
    x = deepcopy(Vector{Union{Float64, Missing}}(x))
    if first(x) |> ismissing
        k = findfirst(!ismissing, x)
        x[1:k-1] .= x[k]
    end
    if last(x) |> ismissing
        k = findlast(!ismissing, x)
        x[k+1:end] .= x[k]
    end
    y = diff(ismissing.(x))
    heads = findall(y .== 1)
    tails = findall(y .== -1) .+ 1

    argitpl = UnitRange.(heads, tails)
    valitpl = LinRange.(x[heads], x[tails], tails - heads .+ 1)
    x[vcat(argitpl...)] .= vcat(valitpl...)

    return x
end

@time include("datacall.jl")

# t= 1950
# for t in 1950:1:2021
#     colorOp = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), Op_[t])
#     cosM = cosM_[t];
#     # θ = 0.995
#     θ = θ_[t]
#     adjM = cosM .≥ θ;
#     G = Graph(adjM)

#     H_ = []
#     _G = deepcopy(G)
#     deg = degree(_G)
#     while any(deg .> 0)
#         k = argmax(deg)
#         deg[k] = -1

#         push!(H_, k => k)
#         push!(H_, (setdiff(_G.fadjlist[k], findall(deg .< 0)) .=> k)...)
#         deg[_G.fadjlist[k]] .= -2
#     end
#     isolated = findall(iszero.(deg))
#     push!(H_, (isolated .=> isolated)...)
#     H_ = Dict(H_)
#     _vol = [count(values(H_) .== k) for k in 1:n]

#     _adjM = zeros(Int64, n,n)
#     for i in 1:n, j in 1:n
#         _adjM[H_[i],H_[j]] += adjM[i,j]
#     end
#     _adjM[diagind(_adjM)] .= 0
#     _adjM = Symmetric(_adjM)
#     findall(deg .== -1)

#     represent = sort(unique(values(H_)))
#     y_ = 10Op_[t]
#     x_ = 10normalize01(log10.(Gc_[t]))

#     i_ = floor.(Int64, y_)
#     j_ = floor.(Int64, x_)
    
#     ijk_ = Dict(zip(i_, j_) .=> 1:first(size(G)))
#     rx_ = deepcopy(x_)
#     rx_[deg .== -2] .= -10
#     # x_[isolated] .= 19 .+ 5rand(length(isolated))

#     if isfile("g32/$(t).png")
#         @warn "g32/$(t).png already exists"
#         break
#     end
#     ISO3vol = ISO3 .* string.(_vol)
#     g1 = graphplot(_adjM, x = rx_, y = y_, ma = 0, la = 0.1, ew = sqrt.(_adjM)) # ew = log10.(_adjM .+ 1))
#     scatter!(rx_, y_, text = ISO3vol, fontsize = 10, color = colorOp, ms = 10sqrt.(_vol),
#     msw = 0, shape = :rect)

#     g2 = graphplot(adjM, x = x_, y = y_, ma = 0, la = 0.1)
#     # g3 = deepcopy(g2)
#     g3 = plot()
#     scatter!(g2, x_, y_, fontsize = 10, text = ISO3, color = colorOp, ms = 15,
#     msw = 0, shape = :rect)

#     for ijk ∈ ijk_
#         k = ijk.second
#         j = ijk.first[2]
#         i = ijk.first[1]
#         _j = j .+ 0.9(0:0.01:1.0)
#         _i = i .+ 0.9normalize01(vec_Pop_[t][k])
#         # plot!(g3, _j, _i, color = colorOp[k], lw = 1.5, fill = i, fa = 0.5)
#         plot!(g3, _j, _i, color = :black, lw = 1.5, fill = i, fa = 0.5)
#     end
    
#     gargs = (; xlims = [-1, 11], ylims = [-1, 11], showaxis = true)
#     plot!(g1; gargs...)
#     plot!(g2; gargs...)
#     plot!(g3; gargs...)

#     plot(g1, g2,  layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
#     # png("g12/$(t).png")
#     plot(g3, g2, layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
#     png("g32/$(t).png")

# end
# p5 = plot(legend = :none, size = [800, 400], margin = 4mm)
# plot!(p5, 1950:2021, [θ_[t] for t in 1950:2021], lw = 2, color = :black,
# xlabel = L"t", ylabel = L"θ^{\ast} (t)", xticks = [1950, 1973, 1997, 2008, 2020])
# plot!(twinx(), 1961:2022, DATA_GDPgrw.GDPgrw, lw = 2, color = :red, ylabel = "GDP growth")
# png("p5-2.png")

# for t in 1950:2021
#     if isfile("yo/$(t).png")
#         @warn "yo/$(t).png already exists"
#         break
#     end
#     scatter(eachcol(yo[t])..., txt = ISO3,
#     title = "t=$t", xlabel = "young", ylabel = "old",
#     # xlims = [0, 0.5], ylims = [0, 0.35],
#     ms = 0, size = [600, 600])
#     png("yo/$(t).png")
# end

# for k in 1:129
#     if isfile("vf/k = $k, $(ISO3[k]).png")
#         @warn "vf/k = $k, $(ISO3[k]).png already exists"
#         break
#     end
#     plot(eachcol(tnsr_yo[:, k, :])...,
#     title = "k = $k, $(ISO3[k])", xlabel = "young", ylabel = "old",
#     color = :black, alpha = (1:72) ./ 72, arrow = true,
#     xlims = [0, 0.5], ylims = [0, 0.35],
#     size = [600, 600])
#     png("vf/k = $k, $(ISO3[k]).png")
# end

include("../../setup/Portables.jl")

cxcy = cat(eachslice(head(tnsr_yo), dims = 3)..., dims = 2)
dxdy = cat(eachslice(diff(tnsr_yo, dims = 3), dims = 3)..., dims = 2)

anchor = eachcol(round.(Int64, 100cxcy))
anchors = unique(anchor)
idx_ancr_ = findall.([[ac == ancr for ac in anchor] for ancr in anchors])

pos = stack(anchors) ./ 100
dir = 3hcat([mean(dxdy[:, idx_ancr], dims = 2) for idx_ancr in idx_ancr_]...)

# quiver(cxcy[1,:], cxcy[2,:], quiver = (dxdy[1,:], dxdy[2,:]),
# color = :black, alpha = 0.1;
# qargs...); png("vf-1")
# quiver(pos[1,:], pos[2,:], quiver = (dir[1,:], dir[2,:]),
# color = :black, arrow = arrow(:closed, :head, 0.5);
# qargs...); png("vf-2")

include("../../DataDrivenModel/src/DDM.jl")

f = SINDy(cxcy', dxdy', λ = 0.01, K = 3)
print(f, ["young", "old"])

test_ = [SINDy(cxcy', dxdy', λ = 0.01, K = 1)]
for k in 2:20
    push!(test_, SINDy(cxcy', dxdy', λ = 0.01, K = k))
end
test_[end]
test_[1] |> propertynames
plot(getproperty.(test_, :MSE), xticks = 1:20, xlabel = "polynomial order (k)", ylabel = "MSE")
print(test_[4], ["y", "o"])
print(test_[6], ["y", "o"])


pos = Base.product(0:.01:.55, 0:.01:0.35) |> collect .|> collect |> vec |> stack
dir = f.(eachcol(pos)) |> stack
qargs = (; xlims = [0, 0.55], ylims = [0, 0.33], xlabel = "young", ylabel = "old", size = [600, 600])
quiver(pos[1,:], pos[2,:], quiver = (dir[1,:], dir[2,:]),
color = :black, arrow=arrow(0.01);
qargs...)
png("vf-3")


_pos = Base.product(0:.001:.55, 0:.001:.35) |> collect .|> collect |> vec |> stack
_dir = f.(eachcol(_pos)) |> stack
# _argfix = findall(norm.(eachcol(_dir)) .< 0.001)

# scatter(_pos[1,_argfix], _pos[2,_argfix], color = :black, ms = 1; qargs...)
# png("fixedpoints")

# heatmap(reshape(log10.(norm.(eachcol(_dir))), 551, 351)', size = [700, 600], legend = :right)
# png("log10dir")

_argx0 = findall(abs.(_dir[1,:]) .< 0.00002)
_argy0 = findall(abs.(_dir[2,:]) .< 0.00002)
ncln = plot(; qargs...)
scatter!(ncln, _pos[1,_argx0], _pos[2,_argx0], color = RGB(0.7, 0, 0), msw = 0, alpha =0.5, ms = 1; qargs...)
scatter!(ncln, _pos[1,_argy0], _pos[2,_argy0], color = RGB(0, 0, 0.7), msw = 0, alpha =0.5, ms = 1; qargs...)
png("nullcleins")

_argx0 = findall(abs.(_dir[1,:]) .< 0.00002);
_argy0 = findall(abs.(_dir[2,:]) .< 0.00002);
_pos[:, _argx0 ∩ _argy0] # 0.248, 0.176
x = 0.266
y = 0.156
xy = x*y

J = ([ 0 1 0 2x y 0 3x^2 2xy y^2 0
       0 0 1 0 x 2y 0 x^2 2xy 3y^2 ] * f.matrix)
eigen(J).values
real.(eigen(J).values)


ic_ = [
    [.5, .01], [.5, .2], [.53, .03],
    [.4, .03],
    [.3, .125], [.33, .1], [.35, .15],
    [.2, .01], [.2, .1], [.2, .2], [.22, .22], [.24, .26], [.23, .21], [.2, .15], [.255, .155], [.275, .155],
    [.1, .05], [.17, .05], [.15, .015], [.19, .02],
    [.05, .1],[.09, .005], [.09, .01],
]
traj__ = []
for ic in ic_
    traj_ = [ic]
    for t in 1:5000
        push!(traj_, traj_[end] + 0.01f(traj_[end]))
    end
    push!(traj__, traj_)
end

d3 = plot3d(lims = [0,1], camera = [90+45, 45], xlabel = "young", ylabel = "old", zlabel = "middle", formatter = (_...) -> "", size = [600, 600])
_ncln = deepcopy(ncln)
for k in eachindex(traj__)
    _x = first.(traj__[k])
    _y = last.(traj__[k])
    plot!(_ncln, _x, _y, arrow = true, color = :black)
    _z = 1 .- (_x .+ _y)
    plot!(d3, _x, _y, _z, arrow = true, color = :black)
end
_ncln
plot(d3)
plot(d3, [1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], color = :black, lw = 2)
k = 1
png("temp")

M_ = Symmetric.(cosine_matrix.([[vec_Pop_[t][k] for t in 1950:2021] for k in 1:129]))
_argsdn = argmin(minimum.(diag.(M_, 1)))
ISO3[_argsdn]
θ = minimum(diag(M_[_argsdn], 1))
diameter_ = []
for k = 1:129
    G = Graph(M_[k] .≥ θ)
    push!(diameter_, diameter(G))
    graphplot(G, names = 1950:2021, size = [600, 600], ms = 0, la = .1,
    title = "k = $k, $(ISO3[k]) | diameter = $(diameter_[k])")
    png("dm/k = $k, $(ISO3[k]) diameter = $(diameter_[k])")
end
scatter(diameter_, Gc_[2021], txt = ISO3, scale = :log10, ma = 0,
ylabel = "GDP per capita in 2021", xlabel = "diameter", size = [600, 600]); png("3")
histogram(diameter_, xlabel = "diameter", color = :black); png("4")
histogram(Gc_[2021], xlabel = "GDP per capita in 2021", color = :black); png("5")

plot.([vec_Pop_[t][argmin(Op_[t])] for t in 1950:2021])[70]
scatter(1950:2021, [ISO3[argmin(Op_[t])] for t in 1950:2021], color = :black, xlabel = "t", title = "youngest")
scatter(1950:2021, [ISO3[argmax(Op_[t])] for t in 1950:2021], color = :black, xlabel = "t", title = "oldest")