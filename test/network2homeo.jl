using ProgressBars
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering, :Distances, :NetworkLayout]
for package in ProgressBar(packages)
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

vec_Pop_ = Dict([t => getproperty.(collect(POPy[t]), :PopTotal) for t in 1950:2021])
cosM_ = Dict([t => vec_Pop_[t] |> cosine_matrix |> Symmetric for t in 1950:2021])
Op_ = Dict([t => getrratio.(vec_Pop_[t]) |> normalize01 for t in 1950:2021])

θ_ = Dict([t => (cosM_[t] |> eachcol .|> maximum |> minimum) for t in 1950:2021])

plot(vec_Pop_[1950][1] ./ maximum(vec_Pop_[1950][1]), fill = 0, fa = 0.5)

t= 1950
# for t in 1950:10:2021
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
    y_ = 30Op_[t]
    x_ = 10log10.(Gc_[t]) .- 25

    rx_ = deepcopy(x_)
    rx_[deg .== -2] .= -10
    # x_[isolated] .= 19 .+ 5rand(length(isolated))

    if !isfile("$(t).png")
        ISO3vol = ISO3 .* string.(_vol)
        g1 = graphplot(_adjM, x = rx_, y = y_, curvature = 0.5, ma = 0, la = 0.1, ew = sqrt.(_adjM)) # ew = log10.(_adjM .+ 1))
        scatter!(rx_, y_, fontsize = 10, text = ISO3vol, color = colorOp, ms = 10sqrt.(_vol),
        xlims = [0, 30], ylims = [-1, 31], showaxis = true,
        msw = 0, shape = :rect, size = [800, 900])

        g2 = graphplot(adjM, x = x_, y = y_, curvature = 0.5, ma = 0, la = 0.1)
        scatter!(x_, y_, fontsize = 10, text = ISO3, color = colorOp, ms = 15,
        xlims = [0, 30], ylims = [-1, 31], showaxis = true,
        msw = 0, shape = :rect, size = [800, 900])

        plot(g1, g2, layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
        png("$(t).png")
    else
        @warn "$(t).png already exists"
    end
# end
p5 = plot(legend = :none, size = [800, 400], margin = 4mm)
plot!(p5, 1950:2021, [θ_[t] for t in 1950:2021], lw = 2, color = :black,
xlabel = L"t", ylabel = L"θ^{\ast} (t)", xticks = [1950, 1973, 1997, 2008, 2020])
plot!(twinx(), 1961:2022, DATA_GDPgrw.GDPgrw, lw = 2, color = :red, ylabel = "GDP growth")
png("p5-2.png")

