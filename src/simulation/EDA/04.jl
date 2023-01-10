include("basic.jl")

diridx = [1,4,9,2,3,5,6,7,8,(10:17)...]
invidx = [1,4,5,2,6,7,8,9,3,(10:17)...]

rslt = rslt_[1]
dead = dead_[1]
mgrn = mgrn_[1]

rslt_gdf = groupby(rslt, :location)
dead_gdf = groupby(dead, :location)

for k ∈ 1:17
    df = rslt_gdf[k]
    # age_l = filter(:age => age ->      age < 15, df) |> ts_sum
    # age_m = filter(:age => age -> 15 ≤ age < 65, df) |> ts_sum
    # age_h = filter(:age => age -> 65 ≤ age     , df) |> ts_sum
    # age_c = cumsum([age_l age_m age_h], dims = 2)
    # age_r = cumsum([(age_l ./ age_m) (age_h ./ age_m)], dims = 2)
    # t_super_aging = 2020 + findfirst((age_h ./ age_c[:,3]) .> .21)
    # push!(t_super_aging_, t_super_aging)

    아기 = filter(:age => age -> age == 0, df) |> ts_sum    
    ddf = dead_gdf[k]
    시체 = ddf |> ts_sum
    
    mtrx_mgrn_ = Dict()
    유입 = Int64[]
    유출 = Int64[]
    for t ∈ 2021:yend
        vctr_mgrn = mgrn[:, "y$t"]
        mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
        mtrx_mgrn .*= (1 .- I(17))
        temp = mtrx_mgrn[diridx, diridx]
        temp[1:3,1:3] .= 0
        temp[4:17,4:17] .= 0
        mtrx_mgrn = temp[invidx, invidx]
        push!(mtrx_mgrn_, t => mtrx_mgrn)
        push!(유입, sum(mtrx_mgrn[k, :]))
        push!(유출, sum(mtrx_mgrn[:, k]))
    end
    
    pp1 = plot(xticks = 0:30:99, title = "$k $(name_location[k])", xlabel = "age", ylabel = "Population", legend = :outertopright)
    plot!(pp1, 0:99, combine(groupby(df, :age), :y2021 => sum => :pyramid).pyramid, label = "2021", lw = 2, color = :black)
    plot!(pp1, 0:99, combine(groupby(df, :age), :y2050 => sum => :pyramid).pyramid, label = "2050", lw = 2, color =   red)
    push!(pyramid_, pp1)

    pp3 = plot(legend = :outertopright)
    plot!(pp3, 2021:yend, 아기 - 시체, lw = 2, fa = .5, fillrange = 0, label = "net growth")
    plot!(pp3, 2021:yend, 유입 - 유출, lw = 2, fa = .5, fillrange = 0, label = "net migration")
    
    plot(pp1, pp3, layout = (2,1), size = (600, 900))
    png("D:/figure/04 $k")
end

# coordinate = DataFrame([
# "서울" 126.9783882  37.5666103
# "부산" 129.0750223  35.1798160
# "대구" 128.6017630  35.8713900
# "인천" 126.7051505  37.4559418
# "광주" 126.8513380  35.1600320
# "대전" 127.3849508  36.3504396
# "울산" 129.3112994  35.5394773
# "세종" 127.2894325  36.4803512
# "경기" 127.5508020  37.4363177
# "강원" 128.3115261  37.8603672
# "충북" 127.6551404  36.7853718
# "충남" 126.8453965  36.6173379
# "전북" 127.2368291  35.6910153
# "전남" 126.9571667  34.9007274
# "경북" 128.9625780  36.6308397
# "경남" 128.2417453  35.4414209
# "제주" 126.5758344  33.4273366
# ], ["locaiton", "lon", "lat"])

# scatter(coordinate.lon, coordinate.lat, text = 1:17,
#     aspect_ratio = 1, size = (400, 600),
#     color = :white, ms = 20)