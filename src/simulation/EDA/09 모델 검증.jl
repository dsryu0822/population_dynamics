include("basic.jl")
# pp09_01 = plot(title = "Total population of Korea", ylabel = "Population", xticks = [2012, 2017, 2021])
# plot!(pp09_01, 2012:2021, ts_sum(rslt)[1:10], label = "Real Data", color = :black, style = :solid, lw = 2)
# plot!(pp09_01, 2012:2021, ts_sum(vldn_rslt), label = "Simulation", color = red   , style = :dash , lw = 2)

# pp09_02 = plot(title = L"L_1" * " error", xlabel = "year", xticks = [2012, 2017, 2021])
# plot!(pp09_02, 2012:2021, abs.(ts_sum(rslt)[1:10] - ts_sum(vldn_rslt)), legend = :none,
# ylabel = "Number (#)", color = red, lw = 2, fillrange = 0, fa = 0.25)
# plot!(twinx(), 2012:2021, 100abs.(ts_sum(rslt)[1:10] - ts_sum(vldn_rslt)) ./ ts_sum(rslt)[1:10], legend = :none,
# ylabel = "Percentage (%)", color = red, lw = 2, fillrange = 0, fa = 0.25)

# plot(pp09_01, pp09_02, layout = (2, 1), size = (800, 400),
# margin = 5mm, right_margin = 5mm, bottom_margin = 5mm, dpi = 200)
# png("G:/figure/09.png")

real = CSV.read("data/KOSIS/population.csv", DataFrame)
real = ifelse.(ismissing.(real), 0, real)
real = trunc.(Int64, vec(sum(Matrix(real[:, Not(1:3)]), dims = 1)))

vldn_ = Dict()
v_ = Dict()
L_ = []
for y ∈ 2000:2020
    push!(vldn_, y => CSV.read("G:/recent/vldn/rslt $y.csv", DataFrame))
    v_[y] = (vldn_[y] |> ts_sum)
    push!(L_, 100maximum(abs.((real - v_[y]) ./ real)))
end

logocolors = Colors.JULIA_LOGO_COLORS
br = range(logocolors.red, logocolors.blue, length=21)
pp09_01 = plot(ylabel = "Total Population")
plot!(pp09_01, 2000:2021, real, shape = :+, style = :solid, label = "Real", color = :black, la = 0)
for (k, y) in enumerate(2000:2020)
    if y == 2000
        _label = "Start in 2000"
    elseif y == 2020
        _label = "Start in 2020"
    else
        _label = :none
    end
    plot!(pp09_01, y:2021, v_[y][k:end], alpha = 0.5, lw = 2, label = _label, color = br[k])
end

pp09_02 = plot(2000:2020, L_, st = :bar, color = br,
    label = :none, xlabel = "Year", ylabel = L"L_{\infty}" * " error (%)", ylims = (0, Inf))

plot(pp09_01, pp09_02, xlims = (1999,2021), layout = (2, 1), xticks = [2000, 2010, 2021])
png("G:/figure/09.png")
