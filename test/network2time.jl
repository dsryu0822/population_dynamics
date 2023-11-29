using ProgressBars
packages = [:CSV, :DataFrames, :LinearAlgebra,
            :Plots, :StatsBase, :LaTeXStrings,
            :Graphs, :GraphRecipes, :SparseArrays, :Colors, :ColorSchemes,
            :Random, :Clustering, :Distances]
for package in ProgressBar(packages)
    @eval using $(package)
end
cm = Plots.cm
default(); default(legend = :none)

cosine_similarity(x, y) = (x'y) / (norm(x)*norm(y))
normalize(v) = (v .- minimum(v)) / (maximum(v) - minimum(v))
getrratio(v) = sum(v[65:end]) / sum(v[[1:15; 65:end]])
getaratio(v) = sum(v[65:end]) / sum(v[1:15])

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

cd(@__DIR__)
DATA_summary = CSV.read("summary.csv", DataFrame)
ISO3 = DATA_summary.ISO3; areaname = DATA_summary.areaname

@time DATA = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale, :PopTotal]);
dropmissing!(DATA)
filter!(:LocTypeName => x -> x == "Country/Area", DATA)

rDATA = filter(:ISO3_code => x -> (x ∈ ISO3), DATA)
ISO3 = unique(rDATA.ISO3_code); areaname = unique(rDATA.Location)

DATAy = collect(groupby(rDATA, :Time));
DATAy = Dict(1950:2021 .=> groupby.(DATAy, :ISO3_code))

vec_Pop_ = getproperty.(collect(DATAy[1950]), :PopTotal);
Op = getrratio.(vec_Pop_) |> normalize

cosM = vec_Pop_ |> cosine_matrix |> Symmetric
Random.seed!(0); stdx = rand(161); stdy = 3Op

color1 = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), Op)
colorV = color1

graphplot(cosM .> 0.995, names = ISO3, x = stdx, y = stdy,
msa = 0, la = 0.1, nodeshape = :rect, nodecolor = colorV,
xlims = [0,1], ylims = [0, 3], fontsize = 10, size = [600, 1200])

cd("//155.230.155.221/ty/Population")

vec_Pop_ = getproperty.(collect(DATAy[2021]), :PopTotal);
cosM = vec_Pop_ |> cosine_matrix |> Symmetric;
findall((cosM .> 0.995)[:, ISO3 .== "ITA"])

ranking = []
Op_ = []
for t ∈ ProgressBar(1950:2021)
    vec_Pop_ = getproperty.(collect(DATAy[t]), :PopTotal);
    cosM = vec_Pop_ |> cosine_matrix |> Symmetric
    Op = getrratio.(vec_Pop_) |> normalize
    push!(Op_, Op)
    push!(ranking, sortperm(Op))

    stdy = 3Op
    color1 = get.(Ref(ColorSchemes.diverging_bwr_55_98_c37_n256), Op)
    colorV = color1
    graphplot(cosM .> 0.995, names = ISO3, x = stdx, y = stdy,
    msa = 0, la = 0.1, nodeshape = :rect, nodecolor = colorV,
    xlims = [0,1], ylims = [0, 3], fontsize = 10, size = [600, 1200])
    png("$t.png")
end
_Op = stack(Op_)
_ranking = stack(ranking)
_ranking[62, :]
p4 = plot(xticks = 1950:20:2030, xlims = [1950, 2025], xlabel = L"t", ylabel = L"O_{p} (t)")
for k in findall(ISO3 .∈ Ref(["KOR", "JPN", "FRA", "SWE", "IND", "CHN", "USA"]))
    plot!(p4, 1950:2021, _Op[k, :], label = ISO3[k], lw = 2)
    annotate!(p4, 2022, _Op[k, end], text(ISO3[k], 7, :left))
end
p4
title!(p4, "Time evolution of the relative old rate of countries")
plot!(p4, dpi = 200); png("Op.png")


vec_Pop_ = Dict([t => getproperty.(collect(DATAy[t]), :PopTotal) for t in 1950:2021])
cosM_ = Dict(1950:2021 .=> vec_Pop_ |> values .|> cosine_matrix .|> Symmetric)
θ_ = 10. .^ (-0.5:-0.01:-4)
Λ_ = []
for θ ∈ ProgressBar(θ_)
    λ_ = []
    for t ∈ 1950:2021
        push!(λ_, eigen(cosM_[t] .> (1 - θ)).values)
    end
    push!(Λ_, pairwise(Euclidean(), stack(λ_)))
end
logΛ = log.(maximum.(Λ_))
plot(1 .- θ_, logΛ,
xlims = [minimum(1 .- θ_), 1], ylims = [3, 4.3], xformatter = x -> "$x",
xticks = [[0.7, 0.8, 0.99]; round.((1 .- θ_)[[argmax(logΛ)]], digits = 3)],
xlabel = L"\theta", ylabel = L"\log \Lambda_{\theta}",
xscale = :log10, color = :black, lw = 2, dpi = 300)
vline!((1 .- θ_)[[argmax(logΛ)]], color = :red)
png("p5")