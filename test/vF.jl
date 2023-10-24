using Plots, LaTeXStrings; default(color = :black, legend = :none)
using CSV, DataFrames

data = CSV.read.("./data/KOSIS/population.csv", DataFrame)

function interpolate(arr)
    _arr = collect.([LinRange(arr[k], arr[k+1], 11) for k in 1:(length(arr)-1)])
    pop!.(_arr)
    return [vcat(_arr...); arr[end]]
end

begin
    function update(u, birth)
        _u = circshift(u .* death, 1)
        _u[1] = birth'u
        return _u
    end

    dt = da = 0.1
    a_ = 0:da:100;
    t_ = 0:dt:100;

    μ(a) = (1 - cos(π*(a^2)/10000))/50
    # plot(a, μ.(a), label = "Age-specific mortality")
    # png("mortality")

    b(a, r=30) = exp((-(a-r)^2)/100)
    # plot(a, b.(a), label = "Age-specific fertility")
    # png("fertility")

    death = exp.(-da*μ.(a_))
    birth = da*b.(a_)
end


# u0 = [vec(stack([combine(groupby(data, :연령별), "2012 년" => sum => "u0").u0 for _ in 1:10])'); 0] ./ 10
# u0 = [(10:-0.01:0.00)...]

u1 = collect(100000*(1:-0.001:0))

_data = [combine(groupby(data, :연령별), "2020 년" => sum => "u0").u0; 0]
u2 = interpolate(_data)*da; sum(u2)

a2 = plot()
r_ = 0:0.01:1
yaxis_ = []
xaxis_ = r_
for r = r_
    # r = 1
    # u0 = r*u1 + (1-r)*u2

    u_ = [u2]
    ODR_ = [sum(u_[end][651:end]) / sum(u_[end][151:650])]
    # 노년부양비(Old-age dependency ratio)
    _birth = birth * u_[end][1]/(((ODR_[end]^(-r))*birth)'u_[end])

    for t = t_[1:(end-1)]
        ODR = sum(u_[end][651:end]) / sum(u_[end][151:650])

        push!(ODR_, ODR)
        push!(u_, update(u_[end], (ODR_[end]^(-r))*_birth))
        # push!(u_, update(u_[end], birth*(ODR^r)))
    end
    a1 = plot(xticks = [0, 20, 65, 100], ylims = (0, maximum(stack(u_))), legend = :best)
    for k in 1:5
        plot!(a1, a_, u_[1 + 200k], alpha = k/5, color = k, label = 1+25*(k-1))
    end
    a1
    plot(t_, ODR_)
    plot(t_, sum.(u_))

    # push!(yaxis_, sum(u_[end]))
    push!(yaxis_, ODR_[end])
    # push!(xaxis_, fill(r, length(yaxis_[end])))
end
scatter!(a2, xaxis_, yaxis_, ylims = (0, Inf), color = 2,
xlabel = "r", ylabel = "ODR after 100 years")