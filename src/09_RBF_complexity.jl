# Same RBF/STLSQ model as 04; vary only the exact BallTree frontier size.
ENV["GKSwstype"] = "100"
include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core",join=true)[[1,2]])
if !isdefined(@__MODULE__, :STLSQresult)
    include("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core/3_SINDy.jl")
end
include(joinpath(@__DIR__,"06_analyze.jl"))
using NearestNeighbors, Random

# Direct evaluation and exact derivatives of the saved constant + linear + Gaussian recipe.
# Verified against the original STLSQresult and central differences for every model.
struct SweepRBF
    centers::Matrix{Float64}
    widths::Vector{Float64}
    coefficients::Matrix{Float64}
end
function rbf_value_jacobian(f::SweepRBF,u)
    d=length(u)
    value=copy(f.coefficients[1,:])
    J=copy(f.coefficients[2:d+1,:]')
    for i in 1:d, j in 1:d
        value[i] += f.coefficients[j+1,i]*u[j]
    end
    for k in axes(f.centers,2)
        distance=0.0
        for j in 1:d
            distance += (u[j]-f.centers[j,k])^2
        end
        phi=exp(-distance/(2*f.widths[k]^2))
        iszero(phi) && continue
        for i in 1:d
            weight=f.coefficients[k+d+1,i]*phi
            value[i] += weight
            for j in 1:d
                J[i,j] -= weight*(u[j]-f.centers[j,k])/f.widths[k]^2
            end
        end
    end
    return value,J
end
function (f::SweepRBF)(u)
    d=length(u)
    value=copy(f.coefficients[1,:])
    for i in 1:d,j in 1:d
        value[i]+=f.coefficients[j+1,i]*u[j]
    end
    for k in axes(f.centers,2)
        distance=0.0
        for j in 1:d
            distance+=(u[j]-f.centers[j,k])^2
        end
        phi=exp(-distance/(2*f.widths[k]^2))
        iszero(phi) && continue
        for i in 1:d
            value[i]+=f.coefficients[k+d+1,i]*phi
        end
    end
    return value
end

function sweep_ode_forecast(f,u,n;reltol=1e-8,abstol=1e-10,dtmax=0.1)
    problem=ODEProblem((u,p,t)->f(u),Float64.(u),(0.0,Float64(n)))
    solution=solve(problem,Tsit5();saveat=1.0:1.0:n,save_start=false,
        reltol,abstol,dtmax,maxiters=1_000_000)
    success=string(solution.retcode)=="Success" && length(solution.u)==n
    prediction=success ? permutedims(stack(solution.u)) : fill(NaN,n,length(u))
    return prediction,string(solution.retcode),solution.stats.naccept,solution.stats.nreject
end

function rbf_finite_jacobian(f::SweepRBF,u)
    J=zeros(length(u),length(u))
    for j in eachindex(u)
        h=min(cbrt(eps(Float64))*max(1,abs(u[j])),minimum(f.widths)/100)
        offset=zeros(length(u));offset[j]=h
        J[:,j]=(f(u+offset)-f(u-offset))/(2h)
    end
    return J
end

# Same STLSQ threshold/refit equations as the local library. Share initial QR across
# outputs and refit only outputs whose active set changed; do not refactor unchanged
# systems. This matters near the sample-count limit. Check against saved library fits.
function sweep_stlsq(data,sysms,recipe;lambda=1e-5)
    started=time()
    input=Matrix(data[:,last(sysms)])
    # stack avoids specializing hcat on hundreds/thousands of positional columns.
    design=stack(dr.func.(dr.funh(input[:,dr.vecv])) for dr in eachrow(recipe))
    @info "Built RBF design" terms=nrow(recipe) seconds=time()-started
    target=Matrix(data[:,first(sysms)])
    scales=norm.(eachcol(design))
    A=design ./ scales'
    coefficients=A\target
    @info "Initial STLSQ solve" terms=nrow(recipe) seconds=time()-started
    previous=falses(size(coefficients))
    converged=false
    for iteration in 1:length(coefficients)+1
        mask=abs.(coefficients) .< lambda .* scales
        if mask==previous
            converged=true
            break
        end
        coefficients[mask].=0
        groups=Dict{BitVector,Vector{Int}}()
        for j in axes(coefficients,2)
            mask[:,j]==previous[:,j] && continue
            active=BitVector(.!mask[:,j])
            push!(get!(groups,active,Int[]),j)
        end
        for (active,columns) in groups
            any(active) && (coefficients[active,columns]=A[:,active]\target[:,columns])
        end
        @info "STLSQ refit" iteration terms=nrow(recipe) active_coefficients=count(!iszero,coefficients) seconds=time()-started
        previous=mask
    end
    @assert converged
    matrix=sparse(coefficients ./ scales)
    retained=[any(!iszero,row) for row in eachrow(matrix)]
    fitted=design*matrix
    mse=mean(abs2,target-fitted)
    result=STLSQresult("SINDy",recipe,recipe[retained,:],matrix,matrix[retained,:],
        mse,length(target)*log(mse)+2*nrow(recipe),rsq(vec(fitted),vec(target)),first(sysms),last(sysms))
    return result,design
end

function sweep_roots(f, starts, bounds; standard_count, discrete=true)
    roots=DataFrame(x=Float64[],y=Float64[],z=Float64[],residual=Float64[],
        in_data_bounds=Bool[],spectral_radius=Float64[],growth_rate=Float64[],
        stability=String[],standard_search=Bool[],jacobian_difference=Float64[],max_real_eigenvalue=Float64[])
    g(u)=discrete ? f(u)-u : f(u)
    Jg(u)=discrete ? last(rbf_value_jacobian(f,u))-I : last(rbf_value_jacobian(f,u))
    for (i,start) in enumerate(starts)
        u=fixedpoint_newton(g,start;jacobian=Jg)
        all(isfinite,u) && norm(g(u),Inf)<=1e-8 || continue
        duplicate=findfirst(j->norm(u-collect(roots[j,[:x,:y,:z]]),Inf)<=1e-5,1:nrow(roots))
        isnothing(duplicate) || continue
        J=last(rbf_value_jacobian(f,u))
        discrepancy=norm(J-rbf_finite_jacobian(f,u),Inf)
        values=eigvals(J)
        rho=discrete ? maximum(abs.(values)) : NaN
        growth=discrete ? log(rho) : maximum(real.(values))
        label=growth < -1e-6 ? "stable" : growth > 1e-6 ? "unstable" : "inconclusive"
        inside=all(lo<=v<=hi for (v,(lo,hi)) in zip(u,bounds))
        push!(roots,(u...,norm(g(u),Inf),inside,rho,growth,label,i<=standard_count,discrepancy,maximum(real.(values))))
    end
    return sort!(roots,[:x,:y,:z])
end

@load joinpath(@__DIR__,"../data.jld2") iso3
vars=data_config.vars
countries=population_data(iso3;data_config...)
X,Y=map_data(countries;vars)
states=vcat([countries[k] for k in sort(collect(keys(countries)))]...)
bounds=[extrema(states[!,v]) for v in vars]
tree=BallTree(X;leafsize=1,reorder=false)
ticks=unique([collect(10:10:200);collect(250:50:500);750;1000;1500;2000;2500;3000;3200;size(X,2)])
ticks=sort(filter(n->n<=length(NearestNeighbors.leaves(tree)),ticks))
width_scale,lambda=1.0,1e-5
discrete=false # ODE comparison experiment; the active model in 04_RBF.jl is an annual map.
sysms=((discrete ? "next_" : "d") .* string.(vars),string.(vars))
data=hcat(DataFrame(X',last(sysms)),DataFrame((discrete ? Y : Y-X)',first(sysms)))
grid=[collect(p) for p in Iterators.product([range(lo,hi,length=7) for (lo,hi) in bounds]...)] |> vec
observed=[collect(states[i,vars]) for i in unique(round.(Int,range(1,nrow(states),length=200)))]
Random.seed!(0)
random_starts=[[lo+rand()*(hi-lo) for (lo,hi) in bounds] for _ in 1:128]
standard_starts=[grid;observed;random_starts]
output_dir=joinpath(@__DIR__,discrete ? "../archive/RBF_complexity" : "../archive/RBF_ODE_complexity")
mkpath(output_dir)
configuration=(data_config=data_config,ticks=ticks,width_scale=width_scale,lambda=lambda,
    standard_starts=standard_starts,data=data,discrete=discrete,
    integration=(method="Tsit5",reltol=1e-8,abstol=1e-10,dtmax=0.1,maxiters=1_000_000))
config_path=joinpath(output_dir,"configuration.jld2")
if isfile(config_path)
    @assert isequal(JLD2.load(config_path,"configuration"),configuration) "Saved sweep configuration differs. Use a new output directory."
else
    JLD2.@save config_path configuration countries bounds
end

previous_roots=Vector{Vector{Float64}}()
summaries=DataFrame[]
all_roots=DataFrame[]
all_predictions=DataFrame[]
blas_threads=BLAS.get_num_threads()
BLAS.set_num_threads(4)
try
    for nleaves in ticks
        checkpoint=joinpath(output_dir,"leaves_$(lpad(nleaves,4,'0')).jld2")
        if isfile(checkpoint)
            saved=JLD2.load(checkpoint)
            summary,roots,prediction_scores=saved["summary"],saved["roots"],saved["prediction_scores"]
            @info "Resuming completed RBF complexity" nleaves
        else
            @info "Starting RBF complexity" nleaves
            started=time()
            leaves=ball_leaves(tree,nleaves)
            balls=NearestNeighbors.treeregion.(leaves)
            centers=stack(ball.center for ball in balls)
            radii=[ball.r for ball in balls]
            widths=width_scale .* max.(radii,1e-6)
            counts=length.(ball_members.(leaves))
            recipe=cook_rbf(sysms,centers,widths;poly=0:1)
            @info "Prepared BallTree basis" nleaves seconds=time()-started
            f4,design=sweep_stlsq(data,sysms,recipe;lambda)
            coefficients=Matrix(f4.matrix)
            f=SweepRBF(centers,widths,coefficients)
            for j in (1,div(size(X,2),2),size(X,2))
                u=X[:,j]
                features=[only(dr.func.(dr.funh(reshape(u,1,:)[:,dr.vecv]))) for dr in eachrow(recipe)]
                @assert isapprox(f(u),vec(features'*coefficients);atol=1e-8,rtol=1e-7)
                @assert isapprox(last(rbf_value_jacobian(f,X[:,j])),rbf_finite_jacobian(f,X[:,j]);atol=1e-5,rtol=1e-4)
            end
            @info "Fitted RBF; searching fixed points" nleaves elapsed=time()-started
            starts=[standard_starts;collect(eachcol(centers));previous_roots]
            roots=sweep_roots(f,starts,bounds;standard_count=length(standard_starts),discrete)
            insertcols!(roots,1,:nleaves=>nleaves)
            roots.nearest_observed_distance=[minimum(norm(collect(r[vars])-collect(s)) for s in eachrow(Matrix(states[:,vars]))) for r in eachrow(roots)]
            prediction_scores=DataFrame(ISO3=String[],nleaves=Int[],rmse=Float64[],failed=Bool[],
                solver_status=String[],accepted_steps=Int[],rejected_steps=Int[])
            for k in sort(collect(keys(countries)))
                v=countries[k]
                actual=Matrix(v[2:end,vars])
                if discrete
                    pred=forecast(f,Matrix(v[1:1,vars]),size(actual,1);discrete=true,lags=1)
                    status,accepted,rejected="map",size(actual,1),0
                else
                    pred,status,accepted,rejected=sweep_ode_forecast(f,collect(v[1,vars]),size(actual,1))
                end
                failed=!all(isfinite,pred)
                push!(prediction_scores,(string(k),nleaves,failed ? Inf : sqrt(mean(abs2,pred-actual)),failed,status,accepted,rejected))
            end
            # Singular values describe numerical identifiability, not statistical generalization.
            singular=svdvals(design)
            tolerance=maximum(size(design))*eps(Float64)*first(singular)
            numerical_rank=count(>(tolerance),singular)
            condition=numerical_rank<size(design,2) ? Inf : first(singular)/last(singular)
            stable=roots.stability .== "stable"
            summary=DataFrame(nleaves=nleaves,candidates=nrow(recipe),active_terms=nrow(f4.recipeF),
                nonzero_coefficients=count(!iszero,coefficients),training_rmse=sqrt(f4.mse),
                rollout_median_rmse=median(prediction_scores.rmse),rollout_mean_rmse=mean(prediction_scores.rmse),
                rollout_failures=count(prediction_scores.failed),roots=nrow(roots),
                stable=count(stable),stable_in_bounds=count(stable .& roots.in_data_bounds),
                stable_outside=count(stable .& .!roots.in_data_bounds),
                standard_stable_in_bounds=count(stable .& roots.in_data_bounds .& roots.standard_search),
                unstable=count(roots.stability .== "unstable"),inconclusive=count(roots.stability .== "inconclusive"),
                min_growth=isempty(roots.growth_rate) ? NaN : minimum(roots.growth_rate),
                min_samples_per_leaf=minimum(counts),max_samples_per_leaf=maximum(counts),
                min_radius=minimum(radii),median_radius=median(radii),
                numerical_rank=numerical_rank,condition_number=condition,
                coefficient_norm=norm(coefficients),search_starts=length(starts),elapsed_seconds=time()-started)
            JLD2.@save checkpoint f4 centers widths counts radii roots summary prediction_scores singular
            CSV.write(joinpath(output_dir,"roots_$(lpad(nleaves,4,'0')).csv"),roots)
            @info "Completed RBF complexity" nleaves stable=summary.stable[1] in_bounds=summary.stable_in_bounds[1] roots=nrow(roots) seconds=time()-started
        end
        global previous_roots=[collect(r[vars]) for r in eachrow(roots)]
        push!(summaries,summary);push!(all_roots,roots);push!(all_predictions,prediction_scores)
        CSV.write(joinpath(output_dir,"summary.csv"),vcat(summaries...))
        CSV.write(joinpath(output_dir,"fixedpoints.csv"),vcat(all_roots...))
        CSV.write(joinpath(output_dir,"country_scores.csv"),vcat(all_predictions...))
    end
finally
    BLAS.set_num_threads(blas_threads)
end
