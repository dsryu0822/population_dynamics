using DataFrames, CSV

YO(pop) = sum.([pop[1:15], pop[65:end]]) ./ sum(pop)
normalize01(x) = (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
# cd(@__DIR__)
# cd("//155.230.155.221/ty/Population")

DATA_summary = CSV.read("G:/summary.csv", DataFrame)
# ---
DATA_gdp = CSV.read("G:/GDPpc1960.csv", DataFrame, drop = ["Country Name"])
@info "DATA_gdp loaded"
rename!(DATA_gdp, "Country Code" => :ISO3)
filter!(:ISO3 => x -> (x ∈ DATA_summary.ISO3), DATA_gdp)
DATA_gdp = dropmissing(stack(DATA_gdp, string.(1960:2022)))
rename!(DATA_gdp, :variable => :Time, :value => :ecnm)
DATA_gdp.Time = parse.(Int64, DATA_gdp.Time)
# ---
DATA_pop = CSV.read("G:/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :AgeGrp, :PopTotal]);
dropmissing!(DATA_pop)
filter!(:LocTypeName => x -> x == "Country/Area", DATA_pop); select!(DATA_pop, Not(:LocTypeName))
@info "DATA_pop loaded"
rename!(DATA_pop, :ISO3_code => :ISO3)
filter!(:ISO3 => x -> (x ∈ DATA_summary.ISO3), DATA_pop)
# ---
GDP = deepcopy(DATA_gdp)
POP = deepcopy(DATA_pop)

transform!(groupby(GDP, :Time), :ecnm => normalize01 => :ecnm)

cd = DataFrame()
for k in eachindex(DATA_summary.ISO3)
    df_k = collect(groupby.(collect(groupby(POP, :ISO3)), :Time)[k])
    iso3 = df_k[1].ISO3[1]
    time = df_k[1].Time[1]:df_k[end-1].Time[1]
    cc = stack(YO.(getproperty.(df_k, :PopTotal)))
    dd = diff(cc, dims = 2)
    append!(cd, DataFrame(; ISO3 = iso3, Time = time, yng = cc[1, 1:(end-1)], old = cc[2, 1:(end-1)], dyng = dd[1, :], dold = dd[2, :]))
end
data = dropmissing(outerjoin(cd, GDP, on = [:ISO3, :Time]))