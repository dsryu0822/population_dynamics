using CSV, DataFrames, GLM, Plots, StatsPlots

cd("./data/australia")

# ---

AFR = CSV.read("rainbow/Australia_age_specific_fertility.csv", DataFrame)
heatmap(Matrix(AFR), size = (400,400))

try
    select!(AFR, Symbol.(1979:2000))
catch
end
Y = DataFrame(Matrix(Matrix(AFR)'), :auto)
rename!(Y, "a" .* string.(15:49))
y = vec(sum(Matrix(Matrix(AFR)'), dims = 2))

# ---

PDN = CSV.read("macrotrends/australia-population-2022-12-06.csv", DataFrame, header = 15)
pop = PDN[30:51, 2]

# ---

FLF = CSV.read("Australian Bureau of Statistics/ABS_LF_1.0.0_M9.3+2..10.AUS.M.csv", DataFrame)
# show(stdout, "text/plain", FLF[12:12:270,[8, 9]])
temp1 = FLF[(1+11):12:(537-262), [8, 3, 9]]
temp2 = FLF[(538+11):12:(end-262), [8, 3, 9]]
# lbf = temp1[:, 3] ./ temp2[:, 3]
lbf = temp1[:, 3] / 10000
cd(@__DIR__); cd(".."); pwd()

# ---

FHE = CSV.read("Australian Government Department of Education/Australian_female_education.csv", DataFrame)
edu = FHE[31:end, 3] / 100

# ------
# ------

year_ = 1979:2000
α_ = (year_ .- 1978) ./ 22
data = DataFrame(; year, y, pop, lbf, edu)

@df data corrplot(cols(2:5), grid = false, size = (800, 800))

# ------

include("lemma.jl")

Pearsonian1.(x, 1, 0.5, 0.5)

predicted_asf = []
begin
p1 = plot(legend = :outerright, xlabel = "age", ylabel = "Age specific fertility", title = "Real Data(Age specific fertility only)")
for (year, ec, α) = zip(year_, eachcol(AFR), α_)
    plot!(p1, 15:49, collect(ec)
    , color = :black, xticks = 15:5:50, label = "$(year)", alpha = α)
end
p2 = plot(legend = :outerright, xlabel = "age", ylabel = "Age specific fertility", title = "Theoretical(Not include age specific fertility)")
# for (year, y0, m1, m2, α) = zip(year_
#     , 10(15 .- data.pop)
#     , (data.lbf * 10)
#     , (data.edu * 10)
#     , α_)
for (year, y0, m1, m2, α) = zip(year_
    , 10(15 .- data.pop)
    , (5 .- (data.edu) * 5)
    , (9 .- (data.lbf) * 15)
    , α_)
# for (year, y0, m1, m2, α) = zip(year_
#     ,  9.70069009401122(14.785514621365198 .- data.pop)
#     , (data.lbf * 9.375111584950243)
#     , (data.edu * 10.713707418033632)
#     , α_)
    push!(predicted_asf, Pearsonian1.(x, y0, m1, m2))
    plot!(p2, x30, Pearsonian1.(x, y0, m1, m2)
    , color = :black, xticks = 15:5:50, label = "$(year)", alpha = α)
end
plot(p1, p2, layout = (2,1), size = (700, 700))
end

TSM_real = Matrix(AFR)
TSM_pred = hcat(predicted_asf...)[1:10:(end-1),:]
p3 = plot(legend = :outerright, size = 75 .* (16, 9))
plot!(p3, [2000],[0], label = "real", color = :black)
plot!(p3, [2000],[0], label = "pred", color = :black, style = :dash)
for i in 5:5:35
plot!(p3, 1979:2000, TSM_real[i,:], label = "age $(i+15)", color = i ÷ 5, lw = 2)
plot!(p3, 1979:2000, TSM_pred[i,:], label = :none, color = i ÷ 5, style = :dash)
end
p3

function L(y0a,y0b,m1a,m1b,m2a,m2b)
    temp = [Pearsonian1.(x, y0a*(y0b - data[i, 3]), m1a - data[i, 4] * m1b, m2a - data[i, 5] * m2b) for i in 1:nrow(data)]
    return sum(abs2, Matrix(AFR) - hcat(temp...)[1:10:(end-1),:])
end
function ∂(f, x...)
    return [
        f(x[1] + h, x[2], x[3], x[4], x[5], x[6]) - f(x...),
        f(x[1], x[2] + h, x[3], x[4], x[5], x[6]) - f(x...),
        f(x[1], x[2], x[3] + h, x[4], x[5], x[6]) - f(x...),
        f(x[1], x[2], x[3], x[4] + h, x[5], x[6]) - f(x...),
        f(x[1], x[2], x[3], x[4], x[5] + h, x[6]) - f(x...),
        f(x[1], x[2], x[3], x[4], x[5], x[6] + h) - f(x...)
    ] / h
end
h = 0.0001
γ = 0.00000001

optimizer = [10, 15, 5, 5, 9, 15]
L(optimizer...)
for i in 1:100
    optimizer -= γ*∂(L, optimizer...)
    println(L(optimizer...))
end
optimizer