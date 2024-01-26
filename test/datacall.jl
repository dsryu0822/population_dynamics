YO(pop) = sum.([pop[1:15], pop[65:end]]) ./ sum(pop)

# cd(@__DIR__)
# cd("//155.230.155.221/ty/Population")

DATA_summary = CSV.read("G:/summary.csv", DataFrame)

DATA_eco = CSV.read("G:/world_GDPpc.csv", DataFrame);
@info "DATA_eco loaded"
dropmissing!(DATA_eco)
rename!(DATA_eco, "Code" => :ISO3)
rename!(DATA_eco, 4 => :GDPpc)
filter!(:ISO3 => x -> x ∈ DATA_summary.ISO3, DATA_eco)
sort!(DATA_eco, :ISO3)

fISO3 = combine(groupby(DATA_eco, :ISO3), nrow)
filter!(:nrow => x -> x == 32, fISO3)
rECO = filter(:ISO3 => x -> (x ∈ fISO3.ISO3), DATA_eco)

ECOy = collect(groupby(rECO, :Year));
ECOy = Dict(1990:2021 .=> ECOy)
Gc_ = Dict([t => getproperty(ECOy[t], :GDPpc) for t in 1990:2021])
push!(Gc_, [t => getproperty(ECOy[1990], :GDPpc) for t in 1950:1989]...)

# ---

DATA_pop = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :Location, :AgeGrp, :PopMale, :PopFemale, :PopTotal]);
@info "DATA_pop loaded"
dropmissing!(DATA_pop)
rename!(DATA_pop, :ISO3_code => :ISO3)
filter!(:LocTypeName => x -> x == "Country/Area", DATA_pop)
sort!(DATA_pop, :ISO3)

rPOP = filter(:ISO3 => x -> (x ∈ fISO3.ISO3), DATA_pop)
areaname = unique(rPOP.Location); ISO3 = unique(rPOP.ISO3); n = length(ISO3)

POPy = collect(groupby(rPOP, :Time));
POPy = Dict(1950:2021 .=> groupby.(POPy, :ISO3, sort = false))
# sort = false is important to keep the order of ISO3
# sort = ture is no nessary

# ---

# DATA_GDPgrw = DataFrame(year = 1961:2022, GDPgrw = CSV.read("C:/world_GDPgrw.csv", DataFrame)[:, 2])
# DATA_GDPgrw

vec_Pop_ = Dict([t => getproperty.(collect(POPy[t]), :PopTotal) for t in 1950:2021])
# cosM_ = Dict([t => vec_Pop_[t] |> cosine_matrix |> Symmetric for t in 1950:2021])
# Op_ = Dict([t => getrratio.(vec_Pop_[t]) |> normalize01 for t in 1950:2021])
# θ_ = Dict([t => (cosM_[t] |> eachcol .|> maximum |> minimum) for t in 1950:2021])

yo = Dict([t => stack(YO.(vec_Pop_[t])) for t in 1950:2021])
findfirst(ISO3 .== "KOR") # 67 = KOR, 63 = JPN
tnsr_yo = cat([yo[t] for t in 1950:2021]..., dims = 3)