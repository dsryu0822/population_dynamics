using CSV, DataFrames, JLD2

rawdata = CSV.read("data/unpopulation_dataportal_20261021214127.csv", DataFrame)

ISO3 = unique(rawdata.Iso3)
location = unique(rawdata.Location)

temp = reshape(rawdata.Value, 21, :)
txyz = rename([rawdata[21:21:end, [:Iso3, :Time]] DataFrame(log10.([sum(temp[1:3, :], dims=1); sum(temp[4:13, :], dims=1); sum(temp[14:21, :], dims=1)])', ["x", "y", "z"])], :Iso3 => :ISO3, :Time => :t)
iso3 = Dict(ISO3 .=> DataFrame.(collect(groupby(txyz, :ISO3))))

JLD2.@save "data.jld2" ISO3 iso3