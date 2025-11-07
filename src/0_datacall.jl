using CSV, DataFrames
cd(@__DIR__); cd("..")

function findrow(df, conditions...)
    cols = first.(conditions)
    queries = last.(conditions)
    bit_ = []
    for (col, query) in zip(cols, queries)
        push!(bit_, (df[:, col] .== query))
    end
    return df[reduce(.*, bit_), :]
end

function diff2(v)
    dv = []
    for k in eachindex(v)
        try push!(dv, v[k+1] - v[k])
        catch e push!(dv, missing) end
    end
    return dv
end
sum_yng(popbyage) = sum(popbyage[1:15])
sum_prd(popbyage) = sum(popbyage[16:65])
sum_old(popbyage) = sum(popbyage[66:end])
# ratio_yng(popbyage) = sum_yng(popbyage) / sum(popbyage)
# ratio_old(popbyage) = sum_old(popbyage) / sum(popbyage)

raw = "G:/population/data/raw"
_ISO3 = CSV.read("$raw/countries_codes_and_coordinates.csv", DataFrame)
rename!(_ISO3, [:country, :alpha2, :ISO3, :country_code, :lat, :lon])
ISO3 = _ISO3[:, 3]
country = Dict(_ISO3[:, 3] .=> _ISO3[:, 1])

if isfile("data/data.csv")
    data = CSV.read("data/data.csv", DataFrame)
else
    POP = CSV.read("$raw/WPP2024_Population1JanuaryBySingleAgeSex_Medium_1950-2023.csv", DataFrame, select = [4, 13, 16, 20])
    rename!(POP, [:ISO3, :t, :age, :pop])
    dropmissing!(POP)
    sort!(POP, [:ISO3, :t, :age])
    # filter!(:ISO3 => x -> x in ISO3, POP)
    @assert all(unique(combine(groupby(POP, [:ISO3, :t]), nrow => :n).n) .== 101)
    # POP = combine(groupby(POP, [:ISO3, :t]), :pop => ratio_yng => :x, :pop => ratio_old => :y)
    # POP = combine(groupby(POP, [:ISO3]), :t, :x, :y,  :x => diff2 => :dx, :y => diff2 => :dy)
    POP = combine(groupby(POP, [:ISO3, :t]), :pop => sum_yng => :x, :pop => sum_old => :y, :pop => sum_prd => :w)
    POP.x .= log10.(POP.x); POP.y .= log10.(POP.y); POP.w .= log10.(POP.w)
    POP = combine(groupby(POP, [:ISO3]), :t, :x, :w, :y, :x => diff2 => :dx, :y => diff2 => :dy)

    GDP = CSV.read("$raw/easy_API_NY.GDP.PCAP.CD_DS2_en_csv_v2_38293.csv", DataFrame)
    GDP = sort(stack(GDP, 2:66))
    rename!(GDP, [:ISO3, :t, :z])
    # dropmissing!(GDP)
    filter!(:ISO3 => x -> x in ISO3, GDP)
    GDP.z .= log10.(GDP.z)
    GDP.t .= parse.(Int, GDP.t)
    GDP = combine(groupby(GDP, :ISO3), :t, :z, :z => diff2 => :dz)

    data = outerjoin(POP, GDP, on = [:ISO3, :t])[:, [:ISO3, :t, :x, :w, :y, :z, :dx, :dy, :dz]]
    sort!(data, [:ISO3, :t])
    CSV.write("data/data.csv", data)
end
