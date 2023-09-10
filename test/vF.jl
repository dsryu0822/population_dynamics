using Plots, LaTeXStrings; default(color = :black, legend = :none)
using CSV, DataFrames

data = CSV.read.("./data/KOSIS/population.csv", DataFrame)

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
_u0 = [combine(groupby(data, :연령별), "2012 년" => sum => "u0").u0; 0]
__u0 = collect.([LinRange(_u0[k], _u0[k+1], 11) for k in 1:100])
pop!.(__u0); u0 = [vcat(__u0...); 0]*da

u1 = collect(100000*(1:-0.001:0))

r_ = 0:0.01:1
yaxis_ = []
xaxis_ = []
for r = r_
    # r = 0

    u_ = [u0]
    foo_ = [sum(u_[end][201:650]) / sum(u_[end][651:end])]
    _birth = birth * u0[1]/(((foo_[end]^r)*birth)'u0)

    for t = t_[1:(end-1)]
        # birth = da*b.(a_)
        # birth = da*b.(a_, r)
        foo = sum(u_[end][201:650]) / sum(u_[end][651:end])

        push!(foo_, foo)
        push!(u_, update(u_[end], (foo_[end]^r)*_birth))
        # push!(u_, update(u_[end], birth*(foo^r)))
    end
    a1 = plot(xticks = [0, 20, 65, 100], ylims = (0, maximum(stack(u_))), legend = :best)
    for k in 1:5
        plot!(a1, a_, u_[1 + 200k], alpha = k/5, color = k, label = 1+25*(k-1))
    end
    a1
    plot(t_, foo_)
    plot(t_, sum.(u_))

    push!(yaxis_, [sum(u_[end][1:200])])
    # push!(yaxis_, foo_[500:100:end])
    push!(xaxis_, fill(r, length(yaxis_[end])))
end
scatter(xaxis_, yaxis_, ylims = (0, Inf))