using CSV, DataFrames, JLD2

"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Population

''''''''''''''''''''''''''''''''''''''''''''''''''"""

rawdata = CSV.read("data/unpopulation_dataportal_20261021214127.csv", DataFrame)

temp = reshape(rawdata.Value, 21, :)
txyz = rename([rawdata[21:21:end, [:Iso3, :Time]] DataFrame(log10.([sum(temp[1:3, :], dims=1); sum(temp[4:13, :], dims=1); sum(temp[14:21, :], dims=1)])', ["x", "y", "z"])], :Iso3 => :ISO3, :Time => :t)
iso3 = Dict(unique(txyz.ISO3) .=> DataFrame.(collect(groupby(txyz, :ISO3))))

for (k, v) in iso3
    if sum(exp10.([v.x[75], v.y[75], v.z[75]])) < 1e+6
        delete!(iso3, k)
    end
end
ISO3 = sort(collect(keys(iso3)))


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                GDP per capita

''''''''''''''''''''''''''''''''''''''''''''''''''"""

gdppc = rename(CSV.read("data/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_1529.csv", DataFrame, header = 4)[:, [2, end-2]], ["ISO3", "g2024"])
gdppc.ISO3 .= String.(gdppc.ISO3)
gdppc = gdppc[gdppc.ISO3 .∈ Ref(ISO3), :]

gdppc.g2024 .= tryparse.(Float64, gdppc.g2024)
gdppc = gdppc[.!isnothing.(gdppc.g2024), :]


"""''''''''''''''''''''''''''''''''''''''''''''''''''

                    Land area

''''''''''''''''''''''''''''''''''''''''''''''''''"""

landarea = rename(CSV.read("data/API_AG.LND.TOTL.K2_DS2_en_csv_v2_2282.csv", DataFrame, header = 4)[:, [2, end-2]], ["ISO3", "a2024"])
landarea.ISO3 .= String.(landarea.ISO3)
landarea = landarea[landarea.ISO3 .∈ Ref(ISO3), :]

landarea.a2024 .= tryparse.(Float64, landarea.a2024)
landarea = landarea[.!isnothing.(landarea.a2024), :]


JLD2.@save "data.jld2" ISO3 iso3 gdppc landarea