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
YO(pop) = sum.([pop[1:15], pop[65:end]]) ./ sum(pop)

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

vec_Pop_ = Dict([t => getproperty.(collect(POPy[t]), :PopTotal) for t in 1950:2021])
cosM_ = Dict([t => vec_Pop_[t] |> cosine_matrix |> Symmetric for t in 1950:2021])
Op_ = Dict([t => getrratio.(vec_Pop_[t]) |> normalize01 for t in 1950:2021])

θ_ = Dict([t => (cosM_[t] |> eachcol .|> maximum |> minimum) for t in 1950:2021])

# t= 1950
for t in 1950:1:2021
    colorOp = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), Op_[t])
    cosM = cosM_[t];
    # θ = 0.995
    θ = θ_[t]
    adjM = cosM .≥ θ;
    G = Graph(adjM)

    H_ = []
    _G = deepcopy(G)
    deg = degree(_G)
    while any(deg .> 0)
        k = argmax(deg)
        deg[k] = -1

        push!(H_, k => k)
        push!(H_, (setdiff(_G.fadjlist[k], findall(deg .< 0)) .=> k)...)
        deg[_G.fadjlist[k]] .= -2
    end
    isolated = findall(iszero.(deg))
    push!(H_, (isolated .=> isolated)...)
    H_ = Dict(H_)
    _vol = [count(values(H_) .== k) for k in 1:n]

    _adjM = zeros(Int64, n,n)
    for i in 1:n, j in 1:n
        _adjM[H_[i],H_[j]] += adjM[i,j]
    end
    _adjM[diagind(_adjM)] .= 0
    _adjM = Symmetric(_adjM)
    findall(deg .== -1)

    represent = sort(unique(values(H_)))
    y_ = 10Op_[t]
    x_ = 10normalize01(log10.(Gc_[t]))

    i_ = floor.(Int64, y_)
    j_ = floor.(Int64, x_)
    
    ijk_ = Dict(zip(i_, j_) .=> 1:first(size(G)))
    rx_ = deepcopy(x_)
    rx_[deg .== -2] .= -10
    # x_[isolated] .= 19 .+ 5rand(length(isolated))

    if isfile("g32/$(t).png")
        @warn "g32/$(t).png already exists"
        break
    end
    ISO3vol = ISO3 .* string.(_vol)
    g1 = graphplot(_adjM, x = rx_, y = y_, ma = 0, la = 0.1, ew = sqrt.(_adjM)) # ew = log10.(_adjM .+ 1))
    scatter!(rx_, y_, text = ISO3vol, fontsize = 10, color = colorOp, ms = 10sqrt.(_vol),
    msw = 0, shape = :rect)

    g2 = graphplot(adjM, x = x_, y = y_, ma = 0, la = 0.1)
    # g3 = deepcopy(g2)
    g3 = plot()
    scatter!(g2, x_, y_, fontsize = 10, text = ISO3, color = colorOp, ms = 15,
    msw = 0, shape = :rect)

    for ijk ∈ ijk_
        k = ijk.second
        j = ijk.first[2]
        i = ijk.first[1]
        _j = j .+ 0.9(0:0.01:1.0)
        _i = i .+ 0.9normalize01(vec_Pop_[t][k])
        # plot!(g3, _j, _i, color = colorOp[k], lw = 1.5, fill = i, fa = 0.5)
        plot!(g3, _j, _i, color = :black, lw = 1.5, fill = i, fa = 0.5)
    end
    
    gargs = (; xlims = [-1, 11], ylims = [-1, 11], showaxis = true)
    plot!(g1; gargs...)
    plot!(g2; gargs...)
    plot!(g3; gargs...)

    plot(g1, g2,  layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
    # png("g12/$(t).png")
    plot(g3, g2, layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
    png("g32/$(t).png")

end
p5 = plot(legend = :none, size = [800, 400], margin = 4mm)
plot!(p5, 1950:2021, [θ_[t] for t in 1950:2021], lw = 2, color = :black,
xlabel = L"t", ylabel = L"θ^{\ast} (t)", xticks = [1950, 1973, 1997, 2008, 2020])
plot!(twinx(), 1961:2022, DATA_GDPgrw.GDPgrw, lw = 2, color = :red, ylabel = "GDP growth")
png("p5-2.png")


yo = Dict([t => stack(YO.(vec_Pop_[t])) for t in 1950:2021])
findfirst(ISO3 .== "KOR") # 67 = KOR, 63 = JPN
tnsr_yo = cat([yo[t] for t in 1950:2021]..., dims = 3)
for t in 1950:2021
    if isfile("yo/$(t).png")
        @warn "yo/$(t).png already exists"
        break
    end
    scatter(eachrow(yo[t])..., txt = ISO3,
    title = "t=$t", xlabel = "young", ylabel = "old",
    # xlims = [0, 0.5], ylims = [0, 0.35],
    ms = 0, size = [600, 600])
    png("yo/$(t).png")
end

for k in 1:129
    if isfile("vf/k = $k, $(ISO3[k]).png")
        @warn "vf/k = $k, $(ISO3[k]).png already exists"
        break
    end
    plot(eachrow(tnsr_yo[:, k, :])...,
    title = "k = $k, $(ISO3[k])", xlabel = "young", ylabel = "old",
    color = :black, alpha = (1:72) ./ 72, arrow = true,
    xlims = [0, 0.5], ylims = [0, 0.35],
    size = [600, 600])
    png("vf/k = $k, $(ISO3[k]).png")
end

head(x) = eval(Meta.parse("$(x)[$(":, " ^ (length(size(x))-1)) 1:(end-1)]"))
tail(x) = eval(Meta.parse("$(x)[$(":, " ^ (length(size(x))-1)) 2: end   ]"))

cxcy = cat(eachslice(head(tnsr_yo), dims = 3)..., dims = 2)
dxdy = cat(eachslice(diff(tnsr_yo, dims = 3), dims = 3)..., dims = 2)

anchor = eachcol(round.(Int64, 100cxcy))
anchors = unique(anchor)
idx_ancr_ = findall.([[ac == ancr for ac in anchor] for ancr in anchors])

pos = stack(anchors) ./ 100
dir = 3hcat([mean(dxdy[:, idx_ancr], dims = 2) for idx_ancr in idx_ancr_]...)

qargs = (; xlims = [0, 0.55], ylims = [0, 0.33], xlabel = "young", ylabel = "old", size = [600, 600])
quiver(cxcy[1,:], cxcy[2,:], quiver = (dxdy[1,:], dxdy[2,:]),
color = :black, alpha = 0.1;
qargs...); png("vf-1")
quiver(pos[1,:], pos[2,:], quiver = (dir[1,:], dir[2,:]),
color = :black, arrow = arrow(:closed, :head, 0.5);
qargs...); png("vf-2")


tempΘ(x) = Θ(x, K = 3)
_Ξ = STLSQ(tempΘ(cxcy'), dxdy')
Ξ = _Ξ.matrix

# nullspace(Matrix(Ξ'))
pos_ = Base.product(0:.001:.55, 0:.001:0.35) |> collect .|> collect |> vec |> stack |> transpose
dir_ = (tempΘ(pos_) * Ξ)

_argfix = findall(norm.(eachrow(dir_)) .< 0.001)
scatter(pos_[_argfix,1], pos_[_argfix,2], color = :black, ms = 1; qargs...)
png("fixedpoints")

heatmap(reshape(log10.(norm.(eachrow(dir_))), 551, 351)', size = [700, 600], legend = :right)
png("log10dir")

quiver(pos_[:,1], pos_[:,2], quiver = (dir_[:,1], dir_[:,2]),
color = :black, arrow=arrow(0.01);
qargs...)
png("vf-3")

using Combinatorics, LinearAlgebra, SparseArrays, DataFrames

struct STLSQresult
    matrix::AbstractMatrix
    MSE::Float64
end
function Base.show(io::IO, ed::STLSQresult)
    show(io, "text/plain", sparse(ed.matrix))
    # println()
    print(io, "\nMSE = $(ed.MSE)")    
end

function STLSQ(ΘX, Ẋ; λ = 10^(-6), verbose = false)::STLSQresult
    if ΘX isa AbstractSparseMatrix
        ΘX = Matrix(ΘX)
    end
    Ξ = ΘX \ Ẋ
    dim = size(Ξ, 2)
    __🚫 = 0
    
    while true
        verbose && print(".")
        🚫 = abs.(Ξ) .< λ
        Ξ[🚫] .= 0
        for j in 1:dim
            i_ = .!🚫[:, j]
            Ξ[i_, j] = ΘX[:,i_] \ Ẋ[:,j]
        end
        if __🚫 == 🚫 verbose && println("Stopped!"); break end # Earl_X stopping
        __🚫 = deepcopy(🚫)
    end
    Ξ =  sparse(Ξ)
    MSE = sum(abs2, Ẋ - ΘX * Ξ) / length(Ẋ)
    verbose && println("MSE = $MSE")

    return STLSQresult(Ξ, MSE)
end
function STLSQ(df::AbstractDataFrame, Ysyms::AbstractVector{T}, Xsyms::AbstractVector{T};
    K = 1, M = 0, f_ = Function[],
    λ = 10^(-6), verbose = false) where T <: Union{Integer, Symbol}
    X = Θ(df[:, Xsyms], K = K, M = M, f_ = f_)
    Y = Matrix(df[:, Ysyms])
    return STLSQ(X, Y, λ = λ, verbose = verbose)
end
function Θ(X::AbstractMatrix; K = 1, M = 0, f_ = Function[])
    dim = size(X, 2)
    ansatz = []

    for k in 0:K
        for case = collect(multiexponents(dim, k))
            push!(ansatz, prod(X .^ case', dims = 2))
        end
    end
    ΘX = hcat(ansatz...)
    for f in f_
        ΘX = [ΘX f.(X)]
    end
    for m in 1:M
        ΘX = [ΘX cospi.(m*X) sinpi.(m*X)]
    end

    return ΘX
end
Θ(X::AbstractVector; K = 1, M = 0, f_ = Function[]) = 
 Θ(reshape(X, 1, :), K = K, M = M, f_ = f_)
Θ(X::AbstractDataFrame; K = 1, M = 0, f_ = Function[]) = 
        Θ(Matrix(X), K = K, M = M, f_ = f_)
  Θ(X::DataFrameRow; K = 1, M = 0, f_ = Function[]) = 
       Θ(collect(X), K = K, M = M, f_ = f_)

println("")

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