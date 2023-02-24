
# d2 = ""
# mtrx_mgrn .*= (1 .- I(17))
# temp = mtrx_mgrn .÷ 20000
# for i = 1:16
#     for j = (i+1):17
#         □ = "<->"
#         # for _ in 1:1
#         for k in 1:min(temp[j,i], temp[i,j])
#             d2 *= "$(name_location[i]) $□ $(name_location[j])\n"
#         end
#         if temp[i,j] - temp[j, i] > 0
#             □ = "<-"
#         else
#             □ = "->"
#         end
#         # for _ in 1:1
#             for k in 1:abs(temp[i,j] - temp[j, i])
#             d2 *= "$(name_location[i]) $□ $(name_location[j])\n"
#         end
#     end
#     println()
# end
# d2 *= "\n"
# for i = 1:17
#     d2 *= """
# $(name_location[i]): {
#     shape: circle
#     height: $(10trunc(Int64, vrtx_wght[i]))
# }
# """
# end
# println(d2)
# write("2012.d2", d2)

clusters = [
    ["서울특별시", "인천광역시", "경기도"]
    , ["대구광역시", "울산광역시", "대전광역시", "세종특별자치시", "광주광역시", "부산광역시"]
    , ["전라남도", "전라북도", "경상북도", "강원도", "충청북도", "충청남도", "경상남도", "제주특별자치도"]
]
@assert sum(length.(clusters)) == 17
ncluster = length(clusters)
namecluster = ["CAP", "MET", "ETC"]

# d2 = "2012 -> 2030 -> 2050 -> 2070\n"
d2 = ""
for y = [5,23,43,63]
    rslt_y = rslt[:, [1,2,3,  4]]
    mgrn_y = mgrn[:, [1,2,3,4,y]]
    # CM = zeros(Int64, ncluster,ncluster)
    POP = []
    CMY = zeros(Int64, ncluster,ncluster)
    CMO = zeros(Int64, ncluster,ncluster)

    for i in 1:ncluster
        push!(POP, trunc(Int64, sqrt(sum(filter(:location => location -> location ∈ clusters[i], rslt_y)[:, end]) / 1000)) )
        for j in 1:ncluster
            flow = filter(:to => to -> to ∈ clusters[j],
                filter(:from => from -> from ∈ clusters[i], mgrn_y)
            )        
            # CM[i,j] = sum(flow[:, end])
            CMY[i,j] = sum(flow[flow.age .< 35, end])
            CMO[i,j] = sum(flow[.!(flow.age .< 35), end])
        end
    end
    # CM - CM'

    d2 *= "$(y+2007) : {\n"
    for i in 1:ncluster
        for j in (i+1):ncluster
            rep = abs((CMO - CMO')[i,j]) ÷ 10000
            if (CMO - CMO')[i,j] > 0
                d2 *= "$(namecluster[i]) -> $(namecluster[j]): {stroke:\"#C00000\"}\n"^rep
            else
                d2 *= "$(namecluster[i]) <- $(namecluster[j]): {stroke:\"#C00000\"}\n"^rep
            end
            rep = abs((CMY - CMY')[i,j]) ÷ 5000
            if (CMY - CMY')[i,j] > 0
                d2 *= "$(namecluster[i]) -> $(namecluster[j]): {stroke:\"#0E64B2\"}\n"^rep
            else
                d2 *= "$(namecluster[i]) <- $(namecluster[j]): {stroke:\"#0E64B2\"}\n"^rep
            end
        end
    end
    d2 *= "\n\n"
    for k in 1:ncluster
        d2 *= "$(namecluster[k]).shape: circle \n"
        d2 *= "$(namecluster[k]).width: $(POP[k]) \n"
    end
    d2 *= "}\n\n"
end
write("temp.d2", d2)
