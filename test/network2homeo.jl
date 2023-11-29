using ProgressBars
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering, :Distances, :NetworkLayout]
for package in ProgressBar(packages)
    @eval using $(package)
end
cm = Plots.cm
default(); default(legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))
normalize(v) = (v .- minimum(v)) / (maximum(v) - minimum(v))
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


cd(@__DIR__)
cd("//155.230.155.221/ty/Population")

DATA_summary = CSV.read("G:/summary.csv", DataFrame)

@time DATA_eco = CSV.read("G:/world_GDPpc.csv", DataFrame);
dropmissing!(DATA_eco)
rename!(DATA_eco, "Code" => :ISO3)
rename!(DATA_eco, 4 => :GDPpc)
filter!(:ISO3 => x -> x ∈ DATA_summary.ISO3, DATA_eco)

fISO3 = combine(groupby(DATA_eco, :ISO3), nrow)
filter!(:nrow => x -> x == 32, fISO3)
rECO = filter(:ISO3 => x -> (x ∈ fISO3.ISO3), DATA_eco)

ECOy = collect(groupby(rECO, :Year));
ECOy = Dict(1990:2021 .=> ECOy)
Gc_ = Dict([t => getproperty(ECOy[t], :GDPpc) for t in 1990:2021])

# ---

@time DATA_pop = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale, :PopTotal]);
dropmissing!(DATA_pop)
rename!(DATA_pop, :ISO3_code => :ISO3)
filter!(:LocTypeName => x -> x == "Country/Area", DATA_pop)
@time sort!(DATA_pop, :ISO3)

rPOP = filter(:ISO3 => x -> (x ∈ fISO3.ISO3), DATA_pop)
areaname = unique(rPOP.Location); ISO3 = unique(rPOP.ISO3); n = length(ISO3)

POPy = collect(groupby(rPOP, :Time));
POPy = Dict(1950:2021 .=> groupby.(POPy, :ISO3, sort = false))
# sort = false is important to keep the order of ISO3
# sort = ture is no nessary

vec_Pop_ = Dict([t => getproperty.(collect(POPy[t]), :PopTotal) for t in 1950:2021])
cosM_ = Dict([t => vec_Pop_[t] |> cosine_matrix |> Symmetric for t in 1950:2021])
Op_ = Dict([t => getrratio.(vec_Pop_[t]) |> normalize for t in 1950:2021])


t = 1990
t = 2021
for t in 1990:2021
θ = 0.995
colorOp = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), Op_[t])
cosM = cosM_[t];
adjM = cosM .> θ;
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

ISO3vol = ISO3 .* string.(_vol)
g1 = graphplot(_adjM, x = rx_, y = y_, curvature = 0.5, ma = 0, ew = log10.(_adjM .+ 1))
scatter!(rx_, y_, fontsize = 10, text = ISO3vol, color = colorOp, ms = 15,
xlims = [0, 30], ylims = [-1, 31], showaxis = true,
msw = 0, shape = :rect, size = [800, 900])

g2 = graphplot(adjM, x = x_, y = y_, curvature = 0.5, ma = 0, ew = log10.(adjM .+ 1))
scatter!(x_, y_, fontsize = 10, text = ISO3, color = colorOp, ms = 15,
xlims = [0, 30], ylims = [-1, 31], showaxis = true,
msw = 0, shape = :rect, size = [800, 900])

plot(g1, g2, layout = (1,2), size = [1600, 900], framestyle = :box, plot_title = "t = $t")
png("$(t).png")
end