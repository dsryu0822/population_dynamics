# Evaluate saved folds only; no fitting and no fixed-point calculations.
include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core", join=true)[[1,2]])
if !isdefined(@__MODULE__, :STLSQresult)
    include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_SINDy.jl")
end
include(joinpath(@__DIR__, "06_analyze.jl"))
using Flux

saved_dir = joinpath(@__DIR__, "../archive/SINDy_MLP")
output_dir = joinpath(saved_dir, "prediction")
mkpath(output_dir)
saved = JLD2.load(joinpath(saved_dir, "results.jld2"))
data, vars, kfold = saved["data"], saved["data_config"].vars, saved["KFOLD"]
countries = sort(unique(data.ISO3))
scores = DataFrame(ISO3=String[], fold=Int[], model=String[], horizon=Int[],
    rmse=Float64[], failed=Bool[], persistence_rmse=Float64[])
one_step = DataFrame(ISO3=String[], fold=Int[], model=String[], nvalid=Int[],
    rmse=Float64[], failed=Bool[], persistence_rmse=Float64[])
trajectories = DataFrame(ISO3=String[], fold=Int[], model=String[], t=Float64[],
    x=Float64[], y=Float64[], z=Float64[])
truth = Dict()

for fold in 1:kfold
    checkpoint = JLD2.load(joinpath(saved_dir, "fold_$(lpad(fold,2,'0')).jld2"))
    # Compile the saved polynomial; avoid rebuilding its recipe at every RK4 stage.
    polynomial = define(Function, checkpoint["sindy"]; fname="saved_sindy_$fold")
    sindy = let polynomial=polynomial
        u -> begin
            du = zeros(length(u))
            polynomial(du,u,nothing,0.0)
            du
        end
    end
    for row in (1,div(nrow(data),2),nrow(data))
        u = collect(data[row,vars])
        @assert isapprox(sindy(u),checkpoint["sindy"](u);atol=1e-8,rtol=1e-7)
    end
    models = ("MLP" => mlp_map(checkpoint["mlp"]), "SINDy" => sindy)
    @info "Loaded saved models" fold
    @assert checkpoint["valid_indices"] == findall(data.fold .== fold)
    for country in countries
        rows = sort(data[data.ISO3 .== country, :], :t)
        U = vcat(Matrix(rows[:,vars]), Matrix(rows[end:end,"next_" .* string.(vars)]))
        years = [rows.t; rows.t[end]+1]
        actual = U[2:end,:]
        truth[country] = DataFrame(hcat(years,U), [:t; vars])
        persistence = repeat(U[1:1,:], size(actual,1), 1)
        valid = rows.fold .== fold
        for (label, model) in models
            discrete = label == "MLP"
            pred = forecast(model, U[1:1,:], size(actual,1); discrete, lags=1)
            for h in (1,5,10,20)
                @assert h <= size(actual,1)
                failed = !all(isfinite, pred[1:h,:])
                # Nonfinite rollouts count as infinite error, never disappear from averages.
                error = failed ? Inf : sqrt(mean(abs2, pred[1:h,:]-actual[1:h,:]))
                push!(scores, (country,fold,label,h,error,failed,
                    sqrt(mean(abs2,persistence[1:h,:]-actual[1:h,:]))))
            end
            for j in axes(pred,1)
                push!(trajectories, (country,fold,label,Float64(years[j+1]),pred[j,:]...))
            end
            if any(valid)
                indices = findall(valid)
                next = vcat([forecast(model,U[j:j,:],1;discrete,lags=1) for j in indices]...)
                target = actual[indices,:]
                failed = !all(isfinite,next)
                push!(one_step,(country,fold,label,length(indices),
                    failed ? Inf : sqrt(mean(abs2,next-target)), failed,
                    sqrt(mean(abs2,U[indices,:]-target))))
            end
        end
    end
    @info "Evaluated saved fold" fold
end

CSV.write(joinpath(output_dir,"rollout_scores.csv"),scores)
CSV.write(joinpath(output_dir,"heldout_one_step.csv"),one_step)
CSV.write(joinpath(output_dir,"trajectories.csv"),trajectories)
summary = combine(groupby(scores,[:model,:horizon]),
    :rmse => median => :median_rmse, :rmse => mean => :mean_rmse,
    :failed => sum => :failures, nrow => :n,
    [:rmse,:persistence_rmse] => ((a,b)->mean(a .< b)) => :beats_persistence_fraction)
heldout_summary = combine(groupby(one_step,:model),
    [:rmse,:nvalid] => ((e,n)->sqrt(sum(n .* e.^2)/sum(n))) => :pooled_rmse,
    :failed => sum => :failed_country_folds, :nvalid => sum => :nvalid,
    [:persistence_rmse,:nvalid] => ((e,n)->sqrt(sum(n .* e.^2)/sum(n))) => :persistence_rmse)
country_scores = combine(groupby(scores[scores.horizon .== 20,:],[:ISO3,:model]),
    :rmse => median => :median_rmse, :rmse => minimum => :best_rmse,
    :rmse => maximum => :worst_rmse, :failed => sum => :failed_folds,
    :persistence_rmse => first => :persistence_rmse)
country_comparison = unstack(select(country_scores,:ISO3,:model,:median_rmse),:model,:median_rmse)
country_comparison.winner = [m < s ? "MLP" : s < m ? "SINDy" : "tie"
    for (m,s) in zip(country_comparison.MLP,country_comparison.SINDy)]
CSV.write(joinpath(output_dir,"rollout_summary.csv"),summary)
CSV.write(joinpath(output_dir,"heldout_summary.csv"),heldout_summary)
CSV.write(joinpath(output_dir,"country_scores.csv"),country_scores)
CSV.write(joinpath(output_dir,"country_comparison.csv"),country_comparison)
display(summary); display(heldout_summary)
display(combine(groupby(country_comparison,:winner),nrow=>:countries))

plt = scatter(country_comparison.MLP,country_comparison.SINDy,
    xlabel="MLP median 20-year RMSE",ylabel="SINDy median 20-year RMSE",
    xscale=:log10,yscale=:log10,
    label="country / region",ms=3,alpha=0.6,title="Saved folds: recursive prediction (log10 population)")
finite_errors = filter(isfinite,[country_comparison.MLP;country_comparison.SINDy])
hi = maximum(finite_errors)
lo = minimum(filter(>(0),finite_errors))
plot!(plt,[lo,hi],[lo,hi],color=:black,linestyle=:dash,label="equal error")
savefig(plt,joinpath(output_dir,"country_comparison.png"))

# Save every country's ten trajectories, with shared vertical limits per coordinate.
mkpath(joinpath(output_dir,"countries"))
for country in countries
    panels = []
    v = truth[country]
    for label in ("MLP","SINDy"), name in vars
        values = [v[!,name]; trajectories[trajectories.ISO3 .== country,name]]
        low, high = extrema(filter(isfinite,values))
        padding = max(0.03,0.05*(high-low))
        p = plot(v.t,v[!,name],color=:black,lw=2,label="observed",
            title="$country $label $name",xlabel="year",ylabel="log10 population",
            ylims=(low-padding,high+padding),legend=false)
        for fold in 1:kfold
            rows = trajectories[(trajectories.ISO3 .== country) .&
                (trajectories.model .== label) .& (trajectories.fold .== fold),:]
            plot!(p,rows.t,rows[!,name],color=palette(:tab10)[fold],alpha=0.6,lw=1)
        end
        push!(panels,p)
    end
    p = plot(panels...,layout=(2,3),size=(1200,650),margin=4Plots.mm)
    savefig(p,joinpath(output_dir,"countries","$country.png"))
end
JLD2.@save joinpath(output_dir,"evaluation.jld2") scores one_step summary heldout_summary country_scores country_comparison trajectories truth
