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

raw = "G:/population/data/raw"
_ISO3 = CSV.read("$raw/countries_codes_and_coordinates.csv", DataFrame)
rename!(_ISO3, [:country, :alpha2, :ISO3, :country_code, :lat, :lon])
ISO3 = _ISO3[:, 3]
country = Dict(_ISO3[:, 3] .=> _ISO3[:, 1])

csvfile = filter(isfile, ["data/data.csv", "C:/Users/rmsms/OneDrive/바탕 화면/world_population.csv"])[1]
csvfile0 = filter(isfile, ["data/data0.csv", "C:/Users/rmsms/OneDrive/바탕 화면/world_population.csv"])[1]
if isfile(csvfile)
    data = CSV.read(csvfile, DataFrame)
    data0 = CSV.read(csvfile0, DataFrame)
else
    POP = CSV.read("$raw/WPP2024_Population1JanuaryBySingleAgeSex_Medium_1950-2023.csv", DataFrame, select = [4, 13, 16, 20])
    rename!(POP, [:ISO3, :t, :age, :pop])
    dropmissing!(POP)
    sort!(POP, [:ISO3, :t, :age])
    @assert all(unique(combine(groupby(POP, [:ISO3, :t]), nrow => :n).n) .== 101)
    POP.pop .*= 1000
    POP = combine(groupby(POP, [:ISO3, :t]), :pop => sum_yng => :x, :pop => sum_old => :y, :pop => sum_prd => :w)


    """'''''''''''''''''''''''''''''''''''''''''''''''''''''

                        log scaled numbers

    '''''''''''''''''''''''''''''''''''''''''''''''''''''"""
    POP.x .= log10.(POP.x); POP.y .= log10.(POP.y); POP.w .= log10.(POP.w)
    POP = combine(groupby(POP, [:ISO3]), :t, :x, :w, :y, :x => diff2 => :dx, :w => diff2 => :dw, :y => diff2 => :dy)

    GDP = CSV.read("$raw/easy_API_NY.GDP.PCAP.CD_DS2_en_csv_v2_38293.csv", DataFrame)
    GDP = sort(stack(GDP, 2:66))
    rename!(GDP, [:ISO3, :t, :z])
    # dropmissing!(GDP)
    filter!(:ISO3 => x -> x in ISO3, GDP)
    GDP.z .= log10.(GDP.z)
    GDP.t .= parse.(Int, GDP.t)
    GDP = combine(groupby(GDP, :ISO3), :t, :z, :z => diff2 => :dz)

    data = outerjoin(POP, GDP, on = [:ISO3, :t])[:, [:ISO3, :t, :x, :w, :y, :z, :dx, :dw, :dy, :dz]]
    sort!(data, [:ISO3, :t])
    CSV.write("data/data.csv", data)

    """'''''''''''''''''''''''''''''''''''''''''''''''''''''

                            ratio

    '''''''''''''''''''''''''''''''''''''''''''''''''''''"""    
    matrix_M = Matrix(POP[:, [:x, :w, :y]])
    M =  DataFrame(matrix_M ./ sum(matrix_M, dims = 2), [:X, :W, :Y])
    dM = combine(groupby([POP M], :ISO3), :ISO3, :t, :X => diff2 => :dX, :W => diff2 => :dW, :Y => diff2 => :dY)
    data0 = outerjoin([POP M], dM, on = [:ISO3, :t])

    sort!(data0, [:ISO3, :t])
    CSV.write("data/data0.csv", data0)
end
