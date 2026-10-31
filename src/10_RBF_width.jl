# Width-by-leaf map experiment. Leave the active 04_RBF.jl model unchanged.
source_path=joinpath(@__DIR__,"09_RBF_complexity.jl")
prefix=first(split(read(source_path,String),"@load joinpath(@__DIR__,\"../data.jld2\") iso3";limit=2))
include_string(@__MODULE__,prefix,source_path)
using Logging

function run_rbf_width_sweep()
    saved=JLD2.load(joinpath(@__DIR__,"../archive/RBF_complexity/configuration.jld2"))
    original=saved["configuration"]
    data=original.data
    countries=saved["countries"]
    bounds=saved["bounds"]
    vars=original.data_config.vars
    sysms=("next_" .* string.(vars),string.(vars))
    X=Matrix(data[:,vars])'
    states=Matrix(vcat([countries[k] for k in sort(collect(keys(countries)))]...)[:,vars])
    tree=BallTree(X;leafsize=1,reorder=false)
    nleaves_values=collect(10:10:100)
    width_scales=[0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,4.0,6.0]
    folder=joinpath(@__DIR__,"../archive/RBF_width")
    mkpath(folder)
    configuration=(data=data,data_config=original.data_config,nleaves_values=nleaves_values,
        width_scales=width_scales,lambda=original.lambda,discrete=true,
        standard_starts=original.standard_starts,width_floor=1e-6)
    config_path=joinpath(folder,"configuration.jld2")
    if isfile(config_path)
        @assert isequal(JLD2.load(config_path,"configuration"),configuration)
    else
        JLD2.@save config_path configuration countries bounds
    end
    summaries,root_tables,country_tables=DataFrame[],DataFrame[],DataFrame[]
    previous_leaf=Dict{Float64,Vector{Vector{Float64}}}()
    original_threads=BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        for nleaves in nleaves_values
            balls=treeregion.(ball_leaves(tree,nleaves))
            centers=stack(ball.center for ball in balls)
            radii=[ball.r for ball in balls]
            previous_width=Vector{Vector{Float64}}()
            for width_scale in width_scales
                identifier="L$(lpad(nleaves,3,'0'))_W$(replace(string(width_scale),'.'=>'p'))"
                path=joinpath(folder,"$identifier.jld2")
                if isfile(path)
                    checkpoint=JLD2.load(path)
                    summary,roots,prediction_scores=checkpoint["summary"],checkpoint["roots"],checkpoint["prediction_scores"]
                else
                    started=time()
                    widths=width_scale .* max.(radii,1e-6)
                    recipe=cook_rbf(sysms,centers,widths;poly=0:1)
                    f4,design=with_logger(NullLogger()) do
                        sweep_stlsq(data,sysms,recipe;lambda=original.lambda)
                    end
                    f=SweepRBF(centers,widths,Matrix(f4.matrix))
                    discrepancies=Float64[]
                    for j in (1,div(size(X,2),2),size(X,2))
                        features=[only(dr.func.(dr.funh(reshape(X[:,j],1,:)[:,dr.vecv]))) for dr in eachrow(recipe)]
                        push!(discrepancies,norm(f(X[:,j])-vec(features'*Matrix(f4.matrix)),Inf))
                    end
                    starts=[original.standard_starts;collect(eachcol(centers));previous_width;get(previous_leaf,width_scale,Vector{Float64}[])]
                    roots=sweep_roots(f,starts,bounds;standard_count=length(original.standard_starts),discrete=true)
                    insertcols!(roots,1,:nleaves=>nleaves,:width_scale=>width_scale)
                    roots.nearest_observed_distance=[minimum(norm(collect(r[vars])-collect(s)) for s in eachrow(states)) for r in eachrow(roots)]
                    roots.relaxation_years=[r.stability=="stable" ? -1/r.growth_rate : NaN for r in eachrow(roots)]
                    prediction_scores=DataFrame(ISO3=String[],nleaves=Int[],width_scale=Float64[],rmse=Float64[],failed=Bool[])
                    for k in sort(collect(keys(countries)))
                        v=countries[k]
                        actual=Matrix(v[2:end,vars])
                        pred=forecast(f,Matrix(v[1:1,vars]),nrow(v)-1;discrete=true,lags=1)
                        failed=!all(isfinite,pred) || !isfinite(sum(abs2,pred-actual))
                        push!(prediction_scores,(string(k),nleaves,width_scale,failed ? Inf : sqrt(mean(abs2,pred-actual)),failed))
                    end
                    singular=svdvals(design)
                    rank=count(>(maximum(size(design))*eps(Float64)*first(singular)),singular)
                    condition=rank<nrow(recipe) ? Inf : first(singular)/last(singular)
                    stable=roots.stability .== "stable"
                    summary=DataFrame(nleaves=nleaves,width_scale=width_scale,
                        roots=nrow(roots),stable=count(stable),stable_in_bounds=count(stable .& roots.in_data_bounds),
                        stable_outside=count(stable .& .!roots.in_data_bounds),
                        standard_stable_in_bounds=count(stable .& roots.in_data_bounds .& roots.standard_search),
                        training_rmse=sqrt(f4.mse),rollout_median_rmse=median(prediction_scores.rmse),
                        rollout_mean_rmse=mean(prediction_scores.rmse),rollout_failures=count(prediction_scores.failed),
                        candidates=nrow(recipe),active_terms=nrow(f4.recipeF),numerical_rank=rank,condition_number=condition,
                        coefficient_norm=norm(f4.matrix),evaluation_disagreement=maximum(discrepancies),
                        search_starts=length(starts),elapsed_seconds=time()-started)
                    JLD2.@save path f4 centers widths radii roots summary prediction_scores singular
                end
                previous_width=[collect(r[vars]) for r in eachrow(roots)]
                previous_leaf[width_scale]=previous_width
                push!(summaries,summary);push!(root_tables,roots);push!(country_tables,prediction_scores)
                CSV.write(joinpath(folder,"summary.csv"),vcat(summaries...))
                CSV.write(joinpath(folder,"fixedpoints.csv"),vcat(root_tables...))
                CSV.write(joinpath(folder,"country_scores.csv"),vcat(country_tables...))
                @info "Width experiment completed" nleaves width_scale stable_inside=summary.stable_in_bounds[1] stable_outside=summary.stable_outside[1] rmse=summary.rollout_median_rmse[1]
            end
        end
    finally
        BLAS.set_num_threads(original_threads)
    end
    roots=vcat(root_tables...)
    CSV.write(joinpath(folder,"stable_fixedpoints_all.csv"),roots[roots.stability .== "stable",:])
    CSV.write(joinpath(folder,"stable_fixedpoints.csv"),roots[(roots.stability .== "stable") .& roots.in_data_bounds,:])
    return vcat(summaries...)
end

width_summary=run_rbf_width_sweep()
