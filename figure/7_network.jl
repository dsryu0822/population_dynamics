include("../../DataDrivenModel/core/header.jl")
include("../src/0_datacall.jl")

cos_similarity(a, b) = a'b / (norm(a) * norm(b))
POP = CSV.read("$raw/WPP2024_Population1JanuaryBySingleAgeSex_Medium_1950-2023.csv", DataFrame, select = [4, 13, 16, 20])
rename!(POP, [:ISO3, :t, :age, :pop])
dropmissing!(POP)
sort!(POP, [:ISO3, :t, :age])

data_2022 = findrow(data, :t => 2022)

pop_ = []
for is3 in unique(data.ISO3)
    push!(pop_, findrow(POP, :ISO3 => is3, :t => 2022).pop)
end

mW = zeros(length(pop_), length(pop_))
for i in eachindex(pop_)
    for j in eachindex(pop_)
        if i ≥ j continue end
        mW[i, j] = cos_similarity(pop_[i], pop_[j])
    end
end
mW = Symmetric(mW)

θ = .99

A = mW .> θ
xx_ = data_2022.x; yy_ = data_2022.y;
# xx_ = vec(sum(A, dims = 2)); yy_ = data_2022.y;
plt_ntwk = plot(size = [2000, 2000], legend = :none)
for i in axes(A, 1)
    for j in axes(A, 2)
        if i ≥ j continue end
        if A[i, j]
            plot!(plt_ntwk, [xx_[i], xx_[j]], [yy_[i], yy_[j]], lw = .5, c = :black, alpha = .3)
        end
    end
end
scatter(plt_ntwk, xx_, yy_, text = data_2022.ISO3, ms = 0, shape = :rect, msw = 0)
