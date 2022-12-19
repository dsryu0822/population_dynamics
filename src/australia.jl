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
TFR = vec(sum(Matrix(Matrix(AFR)'), dims = 2))

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

# ---

FHE = CSV.read("Australian Government Department of Education/Australian_female_education.csv", DataFrame)
edu = FHE[31:end, 3] / 100
cd(@__DIR__); cd(".."); pwd()

# ------
# ------

year_ = 1979:2000
α_ = (year_ .- 1978) ./ 22
data = DataFrame(; year_, TFR, pop, lbf, edu)

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

age = LinRange(-14.9, 19.9, 35)
age = -10:10


p1_1 = plot(xlabel = "age", ylabel = "Age specific fertility", title = "Real Data", legend = :outerright)
p2_1 = plot(xlabel = "age", ylabel = "Age specific fertility", title = "Theoretical, L2 optimized", legend = :outerright)
w_ = []
for i = 1:22
    X = hcat(
        repeat([log(data.pop[i])], length(age))
        , data.lbf[i] .* log.(1 .+ age / a₁)
        , data.edu[i] .* log.(1 .- age / a₂)
    )
    y = log.(AFR[5:(end-10),i])
    w = inv(X'X)*(X')*y
    # plot(y - X*w)
    push!(w_, w)

    plot!(p1_1, 15:49, AFR[:,i]
        , alpha = α_[i], label = "$(year_[i])", color = :black)
    plot!(p2_1, x30, Pearsonian1.(x, data.pop[i] ^ w[1], data.lbf[i] * w[2], data.edu[i] * w[3])
        , alpha = α_[i], label = "$(year_[i])", color = :black)
end
plot(p1_1, p2_1, layout = (2,1), size = (700, 700))

W = hcat(w_...)

using LaTeXStrings
plot(
    plot(year_, W[1,:], color = :black, title = L"\beta_0 (t)"),
    plot(year_, W[2,:], color = :black, title = L"\beta_1 (t)"),
    plot(year_, W[3,:], color = :black, title = L"\beta_2 (t)"),
    layout = (3,1), xlabel = L"t", legend = :none
)
plot(
    plot(year_, data.TFR, color = :black, title = "TFR"),
    plot(year_, data.pop, color = :black, title = "pop"),
    plot(year_, data.lbf, color = :black, title = "lbf"),
    plot(year_, data.edu, color = :black, title = "edu"),
    size = 50 .* (20,9), legend = :none
)
plot(
    plot(year_, data.pop .^ W[1,:], color = :black, title = L"b_0 = pop^{\beta_0(t)}"),
    plot(year_, W[2,:] .* data.lbf, color = :black, title = L"m_1 = \beta_1(t)  \cdot lbf"),
    plot(year_, W[3,:] .* data.edu, color = :black, title = L"m_2 = \beta_2(t)  \cdot edu"),
    layout = (3,1), xlabel = L"t", legend = :none, size = (700, 700)
)

println("end")