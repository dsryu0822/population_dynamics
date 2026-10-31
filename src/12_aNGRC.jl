using DataFrames, LinearAlgebra, Statistics

struct aNGRCresult
    recipe::AbstractDataFrame
    matrix::AbstractMatrix
    mse::Float64
    aic::Float64
    r2::Float64
    vrbl::Vector{String}
    delay::Vector{Int}
    λ::Float64
end

"""
    bake(vrbl; poly = 0:2, delay = 0:1)

Build polynomial terms from delayed variables: x0, y0, x1, y1, ... .
Only the requested degrees and delays are included. Delay is measured in rows.
"""
function bake(vrbl; poly = 0:2, delay = 0:1)
    vrbl = string.(collect(vrbl))
    poly, delay = collect(poly), collect(delay)
    isempty(vrbl) && throw(ArgumentError("vrbl must not be empty."))
    length(unique(vrbl)) == length(vrbl) || throw(ArgumentError("Variable names must be unique."))
    for (name, values) in (("poly", poly), ("delay", delay))
        !isempty(values) && all(v -> v isa Integer && v >= 0, values) ||
            throw(ArgumentError("$name must contain nonnegative integers."))
    end
    poly, delay = sort!(unique(Int.(poly))), sort!(unique(Int.(delay)))
    xnames = [v * string(d) for d in delay for v in vrbl]
    length(unique(xnames)) == length(xnames) || throw(ArgumentError("Delayed variable names overlap."))
    recipe = DataFrame(deg = Int[], term = String[], vecv = Vector{Int}[], tex = String[])
    # Generate intermediate degrees before selecting poly, as in cook.
    terms = [Int[]]
    for d in 0:maximum(poly)
        if d > 0
            terms = [[t; i] for t in terms for i in (isempty(t) ? 1 : last(t)):length(xnames)]
        end
        d in poly || continue
        for t in terms
            push!(recipe, (d, isempty(t) ? "1" : join(xnames[t], "⋅"), t,
                isempty(t) ? "" : join(xnames[t], "*")))
        end
    end
    insertcols!(recipe, 1, :index => 1:nrow(recipe))
    select!(recipe, :index, :term, :vecv, :tex)
    metadata!(recipe, "vrbl", copy(vrbl))
    metadata!(recipe, "delay", copy(delay))
    return recipe
end

"""
    Π(X, recipe)

Expand a chronological M×n time series using a bake recipe. Returns
(M-maximum(delay))×p features; its first row predicts time maximum(delay)+2.
The last row predicts the state immediately after X.
"""
function Π(X::AbstractMatrix, recipe::AbstractDataFrame)
    vrbl, delay = metadata(recipe, "vrbl"), metadata(recipe, "delay")
    size(X, 2) == length(vrbl) || throw(DimensionMismatch("X must have $(length(vrbl)) state columns."))
    d = maximum(delay)
    size(X, 1) > d || throw(ArgumentError("At least $(d+1) chronological rows are required."))
    U = Float64.(X)
    all(isfinite, U) || throw(ArgumentError("State data must be finite."))
    H = hcat([U[d+1-j:end-j, :] for j in delay]...)
    P = ones(size(H, 1), nrow(recipe))
    for (j, r) in enumerate(eachrow(recipe)), i in r.vecv
        P[:, j] .*= H[:, i]
    end
    return P
end
Π(X::AbstractDataFrame, recipe::AbstractDataFrame) = Π(Matrix(X[:, metadata(recipe, "vrbl")]), recipe)

# Group indices preserve the supplied row order within each series.
function angrc_data(data, vrbl, tag)
    if data isa AbstractDataFrame
        column = tag === nothing ? nothing : tag isa Integer ? names(data)[tag] : string(tag)
        column === nothing || column in names(data) || throw(ArgumentError("Unknown tag column: $column"))
        column in vrbl && throw(ArgumentError("The tag column must not be a state variable."))
        U = Matrix(data[:, vrbl])
        tags = column === nothing ? nothing : data[!, column]
    elseif data isa AbstractMatrix
        tag === nothing || tag isa Integer || throw(ArgumentError("A matrix tag must be a column index."))
        tag === nothing || tag in axes(data, 2) || throw(BoundsError(data, (:, tag)))
        U = tag === nothing ? data : data[:, setdiff(collect(axes(data, 2)), [tag])]
        tags = tag === nothing ? nothing : data[:, tag]
    else
        throw(ArgumentError("data must be a matrix or DataFrame."))
    end
    size(U, 2) == length(vrbl) || throw(DimensionMismatch("State columns must match vrbl."))
    groups = tags === nothing ? [collect(axes(U, 1))] : [findall(t -> isequal(t, k), tags) for k in unique(tags)]
    return U, groups
end

"""
    aNGRC(data, vrbl, recipe; tag = nothing, λ = 1e-5)

Fit next-state values by ridge regression, without crossing tag boundaries.
DataFrame tag accepts a column name or position; matrix tag accepts a position.
Matrix state columns (excluding tag) follow vrbl. Rows run oldest to newest.
Groups with fewer than maximum(delay)+2 rows cannot train and are skipped.
"""
function aNGRC(data, vrbl, recipe::AbstractDataFrame; tag = nothing, λ = 1e-5)
    vrbl = string.(collect(vrbl))
    vrbl == metadata(recipe, "vrbl") || throw(ArgumentError("vrbl must match bake's variable names and order."))
    isfinite(λ) && λ >= 0 || throw(ArgumentError("λ must be finite and nonnegative."))
    U, groups = angrc_data(data, vrbl, tag)
    delay = metadata(recipe, "delay")
    d = maximum(delay)
    features, targets = Matrix{Float64}[], Matrix{Float64}[]
    for idx in groups
        length(idx) >= d+2 || continue
        V = Float64.(U[idx, :])
        all(isfinite, V) || throw(ArgumentError("State data must be finite."))
        push!(features, Π(V[1:end-1, :], recipe))
        push!(targets, V[d+2:end, :])
    end
    isempty(features) && throw(ArgumentError("No series has the $(d+2) rows required for training."))
    X, Y = vcat(features...), vcat(targets...)
    all(isfinite, X) || throw(ArgumentError("Polynomial features overflowed."))
    matrix = Symmetric(X'X + λ*I) \ (X'Y)
    error = Y - X*matrix
    mse = sum(abs2, error) / length(Y)
    aic = length(Y)*log(mse) + 2length(matrix)
    total = sum(abs2, Y .- mean(Y, dims = 1))
    r2 = iszero(total) ? NaN : 1 - sum(abs2, error)/total
    return aNGRCresult(deepcopy(recipe), matrix, mse, aic, r2, copy(vrbl), copy(delay), Float64(λ))
end

(s::aNGRCresult)(X::AbstractMatrix) = Π(X, s.recipe) * s.matrix
(s::aNGRCresult)(X::AbstractDataFrame) = Π(X, s.recipe) * s.matrix
Base.show(io::IO, s::aNGRCresult) = print(io, "aNGRC_with_$(nrow(s.recipe))_terms")
function Base.print(io::IO, s::aNGRCresult)
    table = hcat(select(s.recipe, :index, :term), DataFrame(s.matrix, "next_" .* s.vrbl))
    show(io, MIME"text/plain"(), table; allrows = true, allcols = true)
end
Base.print(s::aNGRCresult) = print(stdout, s)

# Example (definitions only; including this file does not fit a model):
# vrbl = [:x, :y, :z]
# cnfg = bake(vrbl; poly = 0:3, delay = [0, 2, 4])
# f = aNGRC(data, vrbl, cnfg; tag = :ISO3, λ = 1e-2)
# print(f)
# pred = f(country[:, vrbl]) # (number of rows - 4) × 3, including one future state

if abspath(PROGRAM_FILE) == @__FILE__
    vrbl = [:x, :y, :z]
    cnfg = bake(vrbl; poly = 0:3, delay = [0, 2, 4])
    f = aNGRC(data, vrbl, cnfg; tag = :ISO3, λ = 1e-2)
    print(f)
    pred = f(Matrix(country[:, vrbl]))
end
