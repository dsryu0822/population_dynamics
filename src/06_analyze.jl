using DataFrames, LinearAlgebra, Statistics, Plots

# Shared by all models. Change data selection here only.
data_config = (nobs = 21, minpop = 1e6, vars = [:x, :y, :z])

# Recent annual observations, selected before constructing delayed samples.
function population_data(iso3; nobs = data_config.nobs,
    minpop = data_config.minpop, vars = data_config.vars)
    @assert nobs >= 3 "At least three annual observations are required."
    data = Dict()
    for (k, raw) in iso3
        v = sort(raw, :t)
        if nrow(v) >= nobs && sum(exp10.(collect(v[end, vars]))) > minpop
            v = v[end-nobs+1:end, [:t; vars]]
            @assert all(diff(v.t) .== 1) "Expected annual observations: $k"
            data[k] = v
        end
    end
    isempty(data) && error("No countries satisfy the data selection.")
    return data
end

# NGRC features [1, current, previous] predict next; never join countries first.
function delayed(data; vars = [:x, :y, :z])
    pairs = map(sort(collect(keys(data)))) do k
        U = Matrix(data[k][:, vars])
        R = selfprod([ones(size(U, 1)-2) U[2:end-1, :] U[1:end-2, :]])
        (R, U[3:end, :])
    end
    return vcat(first.(pairs)...), vcat(last.(pairs)...)
end

# U contains observed states in rows. Return n future annual states, without seeds.
# Continuous and lags=1 maps start at U[end, :]; lags=2 maps use two rows.
function forecast(f, U, n; discrete = false, lags = 2, dt = 0.01)
    @assert n >= 0 && dt > 0
    @assert lags in (1, 2)
    @assert size(U, 1) >= (discrete ? lags : 1)
    pred = fill(NaN, n, size(U, 2))
    u = Float64.(U[end, :])
    previous = Float64.(U[max(1, size(U, 1)-1), :])
    nsteps = ceil(Int, 1/dt)
    h = 1/nsteps
    for j in 1:n
        if discrete
            value = lags == 1 ? vec(f(u)) : vec(f(vcat(previous', u')))
            previous, u = u, value
        else
            for _ in 1:nsteps
                k1 = f(u)
                k2 = f(u + h/2 * k1)
                k3 = f(u + h/2 * k2)
                k4 = f(u + h * k3)
                u += h/6 * (k1 + 2k2 + 2k3 + k4)
                all(isfinite, u) || break
            end
        end
        all(isfinite, u) || break
        pred[j, :] = u
    end
    return pred
end

function numerical_jacobian(f, u)
    J = zeros(length(f(u)), length(u))
    for j in eachindex(u)
        h = cbrt(eps(Float64)) * max(1, abs(u[j]))
        offset = zeros(length(u)); offset[j] = h
        J[:, j] = (f(u + offset) - f(u - offset)) / (2h)
    end
    return J
end

function fixedpoint_newton(g, start; tol = 1e-12, maxiters = 100, jacobian = nothing)
    u = Float64.(collect(start))
    for _ in 1:maxiters
        value = g(u)
        all(isfinite, value) || break
        norm(value, Inf) <= tol && return u
        J = isnothing(jacobian) ? numerical_jacobian(g, u) : jacobian(u)
        all(isfinite, J) || break
        step = pinv(J) * value
        accepted = false
        for j in 0:20
            candidate = u - 0.5^j * step
            trial = g(candidate)
            if all(isfinite, trial) && norm(trial, Inf) < norm(value, Inf)
                u = candidate
                accepted = true
                break
            end
        end
        accepted || break
    end
    return u
end

function stability(f, u; discrete = false, lags = 2, tol = 1e-6)
    @assert lags >= 1
    if discrete && lags > 1
        n = length(u)
        step(v) = [v[n+1:end]; vec(f(permutedims(reshape(v, n, lags))))]
        values = eigvals(numerical_jacobian(step, repeat(u, lags)))
        growth = log(maximum(abs.(values)))
    elseif discrete
        values = eigvals(numerical_jacobian(f, u))
        growth = log(maximum(abs.(values)))
    else
        values = eigvals(numerical_jacobian(f, u))
        growth = maximum(real.(values))
    end
    label = growth < -tol ? "stable" : growth > tol ? "unstable" : "inconclusive"
    return values, growth, label
end

# Multistart search does not guarantee all roots. Stability is local only.
# growth_rate = max Re(eigenvalue) for ODEs, log(spectral radius) for annual maps.
function find_fixedpoints(f, data; discrete = false, lags = 2, vars = [:x, :y, :z],
    ngrid = 7, residual_tol = 1e-8, merge_tol = 1e-5, stability_tol = 1e-6)
    @assert lags >= 1
    g(u) = !discrete ? f(u) : (lags == 1 ? vec(f(u)) : vec(f(repeat(u', lags, 1)))) - u
    bounds = [extrema(data[!, name]) for name in vars]
    starts = Iterators.product([range(lo, hi, length = ngrid) for (lo, hi) in bounds]...)
    roots = Vector{Vector{Float64}}()
    for start in starts
        u = fixedpoint_newton(g, start)
        if all(isfinite, u) && norm(g(u), Inf) <= residual_tol
            all(norm(u - root, Inf) > merge_tol for root in roots) && push!(roots, u)
        end
    end
    result = DataFrame([name => Float64[] for name in vars])
    result.residual = Float64[]
    result.in_data_bounds = Bool[]
    result.eigenvalues = Vector{ComplexF64}[]
    result.growth_rate = Float64[]
    result.stability = String[]
    for u in roots
        values, growth, label = stability(f, u; discrete, lags, tol = stability_tol)
        in_bounds = all(lo <= x <= hi for (x, (lo, hi)) in zip(u, bounds))
        push!(result, (u..., norm(g(u), Inf), in_bounds, ComplexF64.(values), growth, label))
    end
    return sort!(result, vars)
end

# State-to-next-state samples, constructed within each country.
function map_data(data; vars = data_config.vars)
    states = [Matrix(data[k][:, vars]) for k in sort(collect(keys(data)))]
    X = vcat([U[1:end-1, :] for U in states]...)'
    Y = vcat([U[2:end, :] for U in states]...)'
    return Matrix(X), Matrix(Y)
end

mlp_loss(model, x, y) = mean(abs2, model(x) - y)

# Paired rows keep country/year identity and share targets across model families.
function transition_data(countries; vars = data_config.vars)
    tables = []
    for k in sort(collect(keys(countries)))
        v = countries[k]
        X, Y = Matrix(v[1:end-1, vars]), Matrix(v[2:end, vars])
        push!(tables, hcat(DataFrame(ISO3 = fill(string(k), size(X, 1)), t = v.t[1:end-1]),
            DataFrame(X, vars), DataFrame(Y, "next_" .* string.(vars)),
            DataFrame((Y-X) ./ diff(v.t), "d" .* string.(vars))))
    end
    return vcat(tables...)
end

# Fit scaling on training rows only. Save weights/statistics rather than closures.
function fit_mlp(X, Y; epochs = 5000, seed = 0, width = 32, rate = 1e-3)
    μX, σX = vec(mean(X, dims = 2)), max.(vec(std(X, dims = 2)), 1e-8)
    μY, σY = vec(mean(Y, dims = 2)), max.(vec(std(Y, dims = 2)), 1e-8)
    Xn, Yn = Float32.((X .- μX) ./ σX), Float32.((Y .- μY) ./ σY)
    Random.seed!(seed)
    model = Flux.Chain(Flux.Dense(size(X, 1) => width, tanh),
        Flux.Dense(width => width, tanh), Flux.Dense(width => size(Y, 1)))
    optimizer = Flux.setup(Flux.Adam(rate), model)
    best = deepcopy(model)
    losses = [mlp_loss(model, Xn, Yn)]
    best_loss = first(losses)
    for epoch in 1:epochs
        Flux.train!(mlp_loss, model, [(Xn, Yn)], optimizer)
        loss = mlp_loss(model, Xn, Yn)
        if loss < best_loss
            best, best_loss = deepcopy(model), loss
        end
        push!(losses, loss)
        epoch % 1000 == 0 && @info "MLP training" epoch loss best_loss
    end
    return (model = Flux.f64(best), input_mean = μX, input_scale = σX,
        output_mean = μY, output_scale = σY, losses = losses)
end

function mlp_map(run)
    return Flux.Chain(x -> (x .- run.input_mean) ./ run.input_scale,
        run.model, y -> y .* run.output_scale .+ run.output_mean)
end

# Greedy one-to-one matching in the original log10 coordinates, within each fold.
function match_fixedpoints(points; tolerance = 0.15)
    matches = DataFrame(fold = Int[], mlp_id = Int[], sindy_id = Int[], distance = Float64[])
    for fold in sort(unique(points.fold))
        a = points[(points.fold .== fold) .& (points.model .== "MLP"), :]
        b = points[(points.fold .== fold) .& (points.model .== "SINDy"), :]
        candidates = [(norm(collect(a[i, [:x,:y,:z]])-collect(b[j, [:x,:y,:z]])), i, j)
            for i in 1:nrow(a) for j in 1:nrow(b)]
        used_a, used_b = Set{Int}(), Set{Int}()
        for (distance, i, j) in sort(candidates)
            if distance <= tolerance && !(i in used_a) && !(j in used_b)
                push!(matches, (fold, a.id[i], b.id[j], distance))
                push!(used_a, i); push!(used_b, j)
            end
        end
    end
    return matches
end

# Main figures show stable roots inside the common observed-variable ranges.
# All roots, including roots outside these ranges, remain in the saved tables.
function fold_fixedpoint_plots(states, points, matches; vars = [:x,:y,:z], kfold = 10)
    bounds = [extrema(states[!, v]) for v in vars]
    stable = points[(points.stability .== "stable") .& points.in_data_bounds, :]
    background = states[1:5:end, :]
    panels = []
    for fold in 1:kfold
        a = stable[(stable.fold .== fold) .& (stable.model .== "MLP"), :]
        b = stable[(stable.fold .== fold) .& (stable.model .== "SINDy"), :]
        outside = count((points.fold .== fold) .& (points.stability .== "stable") .& .!points.in_data_bounds)
        plt = scatter(background[!,vars[1]], background[!,vars[2]], background[!,vars[3]],
            color = :gray, alpha = 0.08, ms = 1, msw = 0, label = :none,
            xlims = bounds[1], ylims = bounds[2], zlims = bounds[3],
            xlabel = string(vars[1]), ylabel = string(vars[2]), zlabel = string(vars[3]),
            camera = (35,25), title = "fold $fold | M $(nrow(a)), S $(nrow(b)) | outside $outside",
            titlefontsize = 9, guidefontsize = 9, tickfontsize = 7,
            legend = fold == 1 ? :topleft : false)
        scatter!(plt, a.x, a.y, a.z, color = :dodgerblue, shape = :circle, ms = 4, label = "MLP")
        scatter!(plt, b.x, b.y, b.z, color = :orangered, shape = :utriangle, ms = 4, label = "SINDy")
        for row in eachrow(matches[matches.fold .== fold, :])
            p = points[findfirst(==(row.mlp_id), points.id), :]
            q = points[findfirst(==(row.sindy_id), points.id), :]
            plot!(plt, [p.x,q.x], [p.y,q.y], [p.z,q.z], color = :black,
                linestyle = :dash, lw = 1, label = :none)
        end
        push!(panels, plt)
    end
    by_fold = plot(panels..., layout = (2,5), size = [1800,850], margin = 4Plots.mm)
    panels = []
    colors = palette(:tab10, kfold)
    for model in ("MLP", "SINDy"), (u,v) in ((:x,:y),(:x,:z),(:y,:z))
        plt = scatter(background[!,u], background[!,v], color = :gray, alpha = 0.08,
            ms = 1, msw = 0, label = :none, title = "$model: $u / $v",
            xlabel = string(u), ylabel = string(v), xlims = extrema(states[!,u]),
            ylims = extrema(states[!,v]),
            legend = model == "MLP" && u == :x && v == :y ? :outerright : false,
            legendfontsize = 7)
        for fold in 1:kfold
            rows = stable[(stable.fold .== fold) .& (stable.model .== model), :]
            scatter!(plt, rows[!,u], rows[!,v], color = colors[fold],
                shape = model == "MLP" ? :circle : :utriangle, ms = 5, label = "fold $fold")
        end
        push!(panels, plt)
    end
    projections = plot(panels..., layout = (2,3), size = [1400,850],
        left_margin = 5Plots.mm, bottom_margin = 4Plots.mm)
    return by_fold, projections
end

# Cut an existing BallTree into exactly nleaves terminal regions.
# Split the largest available ball first; sample counts may differ by region.
function ball_leaves(tree, nleaves::Integer)
    @assert 1 <= nleaves <= length(NearestNeighbors.leaves(tree))
    nodes = [NearestNeighbors.treeroot(tree)]
    while length(nodes) < nleaves
        candidates = [i for i in eachindex(nodes) if !isempty(NearestNeighbors.children(nodes[i]))]
        i = candidates[argmax([NearestNeighbors.treeregion(nodes[j]).r for j in candidates])]
        children = NearestNeighbors.children(nodes[i])
        deleteat!(nodes, i)
        append!(nodes, children)
    end
    return nodes
end

function ball_members(node)
    children = NearestNeighbors.children(node)
    isempty(children) && return collect(NearestNeighbors.leaf_point_indices(node))
    return vcat(ball_members.(children)...)
end

# Named callable keeps recipe functions serializable when f4 is saved.
if !isdefined(@__MODULE__, :GaussianRBF)
    struct GaussianRBF
        center::Vector{Float64}
        width::Float64
    end
end
(g::GaussianRBF)(u) = exp(-sum(abs2, u .- g.center) / (2g.width^2))

# Extend the existing cook/Θ recipe with multivariate Gaussian candidates.
function cook_rbf(vrbl, centers, widths; poly = 0:1)
    @assert size(centers, 2) == length(widths) && all(widths .> 0)
    cnfg = cook(vrbl; poly)
    xnames = string.(last(vrbl))
    @assert size(centers, 1) == length(xnames)
    for (j, (c, w)) in enumerate(zip(eachcol(centers), widths))
        term = "RBF$j"
        distance = join(["($(xnames[i])-($(c[i])))^2" for i in eachindex(c)], "+")
        tex = "exp(-($distance)/(2*($w)^2))"
        push!(cnfg, (nrow(cnfg)+1, term, GaussianRBF(collect(c), w), eachrow,
            collect(eachindex(xnames)), tex))
    end
    return cnfg
end

# v and all prediction matrices share the same rows (forecast years only).
function prediction_plot(v, predictions; vars = [:x, :y, :z], title = "")
    panels = []
    for (j, name) in enumerate(vars)
        plt = plot(v.t, v[!, name], color = :black, label = "actual",
            title = "$title, $name", xlabel = "year", ylabel = "log10 population",
            left_margin = 5Plots.mm, bottom_margin = 3Plots.mm)
        for (label, pred) in predictions
            plot!(plt, v.t, pred[:, j], label = label)
        end
        push!(panels, plt)
    end
    return plot(panels..., layout = (1, length(vars)), size = [1100, 350])
end
