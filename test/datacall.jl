using DataFrames, CSV

if "data/cached.csv" |> isfile
    data = CSV.read("data/cached.csv", DataFrame)
    @info "cached.csv loaded"
else
    YO(pop) = sum.([pop[1:15], pop[66:end]]) ./ sum(pop)
    normalize01(x) = (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
    # cd(@__DIR__)
    # cd("//155.230.155.221/ty/Population")

    DATA_summary = CSV.read("G:/population/summary.csv", DataFrame)
    # ---
    DATA_gdp = CSV.read("G:/population/GDPpc1960.csv", DataFrame, drop = ["Country Name"])
    @info "DATA_gdp loaded"
    rename!(DATA_gdp, "Country Code" => :ISO3)
    filter!(:ISO3 => x -> (x ∈ DATA_summary.ISO3), DATA_gdp)
    DATA_gdp = stack(DATA_gdp, string.(1960:2022))
    rename!(DATA_gdp, :variable => :Time, :value => :ecnm)
    DATA_gdp.Time = parse.(Int64, DATA_gdp.Time)
    # ---
    DATA_pop = CSV.read("G:/population/world_population.csv", DataFrame, select = [:Time, :LocTypeName, :ISO3_code, :AgeGrp, :PopTotal]);
    dropmissing!(DATA_pop)
    filter!(:LocTypeName => x -> x == "Country/Area", DATA_pop); select!(DATA_pop, Not(:LocTypeName))
    @info "DATA_pop loaded"
    rename!(DATA_pop, :ISO3_code => :ISO3)
    filter!(:ISO3 => x -> (x ∈ DATA_summary.ISO3), DATA_pop)
    # ---
    GDP = deepcopy(DATA_gdp)
    POP = deepcopy(DATA_pop)

    dropmissing!(GDP)
    # transform!(groupby(GDP, :Time), :ecnm => (x -> normalize01(log10.(x))) => :ecnm)
    GDP = [GDP zeros(nrow(GDP))]; rename!(GDP, :x1 => :decnm)
    for i1 = 1:nrow(GDP)
        iso3 = GDP.ISO3[i1]
        time = GDP.Time[i1]
        i2 = findfirst((GDP.ISO3 .== iso3) .&& (GDP.Time .== (time+1)))
        if i2 |> isnothing
            GDP.decnm[i1] = NaN
        else
            GDP.decnm[i1] = GDP.ecnm[i2] - GDP.ecnm[i1]
        end
    end
    GDP = GDP[.!isnan.(GDP.decnm), :]
    
    cd = DataFrame()
    for k = eachindex(DATA_summary.ISO3)
        df_k = collect(groupby.(collect(groupby(POP, :ISO3)), :Time)[k])
        iso3 = df_k[1].ISO3[1]
        time = df_k[1].Time[1]:df_k[end].Time[1]
        cc = stack(YO.(getproperty.(df_k, :PopTotal)))
        dd = diff(cc, dims = 2)
        append!(cd, DataFrame(; ISO3 = iso3, Time = time[1:(end-1)], yng = cc[1, 1:(end-1)], old = cc[2, 1:(end-1)], dyng = dd[1, :], dold = dd[2, :]))
        # include("../../DataDrivenModel/core/FDM.jl")
        # dd = stack(fdiff.(eachrow(cc); stencil = 3))
        # append!(cd, DataFrame(; ISO3 = iso3, Time = time[2:(end-1)], yng = cc[1, 2:(end-1)], old = cc[2, 2:(end-1)], dyng = dd[:, 1], dold = dd[:, 2]))
    end
    data = outerjoin(GDP, cd, on = [:ISO3, :Time])
    sort!(data, [:ISO3, :Time])
    CSV.write("data/cached.csv", data)
end

using Plots; default(color = :black, legend = :none)

for df = groupby(GDP, :Time)
    t = df.Time[1]
    histogram(df.ecnm, title = t, xlims = [0, 1.2e+5])
    png("G:/population/eda/GDPpc_hist(t)/t=$(t).png")
end

gdf = groupby(GDP, :ISO3)
p1_ = []
p1 = plot()
lp1 = plot()
for df = gdf
    plot(df.Time, df.ecnm, title = last(df.ISO3))
    png("G:/population/eda/GDPpc_timeseries/$(last(df.ISO3)).png")
    plot!(p1, df.Time, df.ecnm, alpha = .1)
    plot!(lp1, df.Time, df.ecnm, alpha = .1, yscale = :log10)
end
p1_[11]
p1
png(lp1, "logtimeseries.png")
lp1