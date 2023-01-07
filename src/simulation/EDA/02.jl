
MOBILITY = CSV.read("data/KOSIS/mobility.csv", DataFrame)
rename!(MOBILITY, ["from", "to", "gender", "age", years...])
MOBILITY = MOBILITY[:, [1,2,3,4,end]]

migration_matrix = reshape(sum(reshape(MOBILITY.y2021, 17,2,17,17), dims = 1:2), 17, 17)
for d ∈ 1:17 migration_matrix[d,d] = 0 end
# 세로 전출, 가로 전출 ex) 서울 → 전출: 13078 in 2,1 of matrix

heatmap(migration_matrix, size = (800, 800), ticks = [4,9,17], xlabel = "from", ylabel = "to")

이동수지 = Int64[]
for d in 1:17
    push!(이동수지, sum(migration_matrix[d,:]) - sum(migration_matrix[:,d]))
end
이동수지

# 시작 시뮬레이션 타당성 검토
k = 1
t = 4
인구수지 = vec(sum(reshape(rslt_[k][:, t+1], 200, :), dims = 1)) - 
 vec(sum(reshape(rslt_[k][:, t], 200, :), dims = 1))

신생아 = vec(sum(reshape(rslt_[k][:, t+1], 100, 2, :)[1,:,:], dims = 1))
사망자 = vec(sum(reshape(dead_[k][:, t], 200, :), dims = 1))
이동수지 + 신생아 - 사망자
plot(인구수지, label = "A")
plot!(이동수지 + 신생아 - 사망자, label = "B")
# 끝 시뮬레이션 타당성 검토