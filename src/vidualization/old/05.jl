using Shapefile, GeoInterface, Proj
using Plots, ColorSchemes

function shp_cache(shptable)
    plgn_ = []
    for (k, plgn) in enumerate(shptable.geometry)
        sub_plgn_ = []
        for sub_plgn ∈ GeoInterface.coordinates(plgn)
            push!(sub_plgn_, hcat(collect.(trans.(sub_plgn[1]))...))
        end
        push!(plgn_, sub_plgn_)
    end
    return plgn_;
end

function normalize01(v)
    v = (v .- minimum(v))
    return v / maximum(v)
end

trans = Proj.Transformation(
    "+proj=tmerc +lat_0=38 +lon_0=127.5 +k=0.9996 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +units=m +no_defs",
    "+proj=lonlat +ellps=WGS84 +datum=WGS84 +nodefs")
    
shptable = Shapefile.Table("src/simulation/EDA/map/ctp_rvn.shp")
plgn_ = shp_cache(shptable);

for t in [2012, 2031, 2041, 2050]
pop = combine(groupby(rslt_[1], :location), "y$t" => sum)[:, 2] |> normalize01
colormap = get(ColorSchemes.algae, pop)
pp1 = plot(legend = :bottomright, title = "year: $t", size = (500, 600),
    xticks = 125:130, yticks = 33:39, aspect_ratio = 1.2)
for k = (17:-1:1)
    sub_plgn_ = plgn_[k]
    for sub_plgn ∈ sub_plgn_
        plot!(pp1, eachrow(sub_plgn)..., fillrange = minimum(sub_plgn),
            color = :black, fc = colormap[k], label = :none)
    end
end
# pp1
png(pp1, "G:/figure/05 y$t.png")
println("G:/figure/05 y$t.png")
end
