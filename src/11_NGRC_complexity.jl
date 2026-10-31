# Separate NGRC experiment. Do not include or modify the user's 02_NGRC.jl.
ENV["GKSwstype"]="100"
include.(readdir("C:/Users/rmsms/OneDrive/lab/DataDrivenModel/core",join=true)[[1,2,5]])
include(joinpath(@__DIR__,"06_analyze.jl"))
using SHA

struct PolynomialBasis
    parents::Vector{Int}
    variables::Vector{Int}
    terms::Vector{Tuple}
end
function polynomial_basis(nvariables,degree)
    terms=Tuple[()]; parents=[0]; variables=[0]
    for d in 1:degree
        previous=findall(t->length(t)==d-1,terms)
        for i in previous, v in (isempty(terms[i]) ? 1 : last(terms[i])):nvariables
            push!(terms,(terms[i]...,v));push!(parents,i);push!(variables,v)
        end
    end
    return PolynomialBasis(parents,variables,terms)
end
function polynomial_features(b,h)
    phi=ones(eltype(h),length(b.parents))
    for i in 2:length(phi)
        phi[i]=phi[b.parents[i]]*h[b.variables[i]]
    end
    return phi
end
function polynomial_value_jacobian(b,W,h)
    phi=polynomial_features(b,h)
    value=vec(W'*phi)
    J=zeros(size(W,2),length(h))
    adjoint=copy(W)
    for i in length(phi):-1:2
        p,v=b.parents[i],b.variables[i]
        for j in axes(W,2)
            J[j,v]+=adjoint[i,j]*phi[p]
            adjoint[p,j]+=adjoint[i,j]*h[v]
        end
    end
    return value,J
end
struct PolynomialNGRC
    lags::Int
    degree::Int
    basis::PolynomialBasis
    coefficients::Matrix{Float64}
    fixed_basis::PolynomialBasis
    fixed_coefficients::Matrix{Float64}
    lambda::Float64
end
function polynomial_ngrc(lags,degree,basis,W,lambda)
    fixed=polynomial_basis(3,degree)
    lookup=Dict(t=>i for (i,t) in enumerate(fixed.terms))
    coefficients=zeros(length(fixed.terms),3)
    for (i,t) in enumerate(basis.terms)
        term=Tuple(sort([mod1(v,3) for v in t]))
        coefficients[lookup[term],:]+=W[i,:]
    end
    return PolynomialNGRC(lags,degree,basis,W,fixed,coefficients,lambda)
end
(f::PolynomialNGRC)(h)=vec(f.coefficients'*polynomial_features(f.basis,h))
fixed_map(f,u)=vec(f.fixed_coefficients'*polynomial_features(f.fixed_basis,u))
function delayed_jacobian(f,u)
    _,first_block=polynomial_value_jacobian(f.basis,f.coefficients,repeat(u,f.lags))
    J=zeros(3*f.lags,3*f.lags)
    J[1:3,:]=first_block
    if f.lags>1
        J[4:end,1:end-3]=Matrix{Float64}(I,3*f.lags-3,3*f.lags-3)
    end
    return J
end
function ngrc_samples(countries,lags;start_index=4,vars=[:x,:y,:z])
    histories,targets=[],[]
    for k in sort(collect(keys(countries)))
        U=Matrix(countries[k][:,vars])
        push!(histories,stack([vcat([U[t-j,:] for j in 0:lags-1]...) for t in start_index:size(U,1)-1];dims=1))
        push!(targets,U[start_index+1:end,:])
    end
    return vcat(histories...),vcat(targets...)
end
function ngrc_design(basis,H)
    R=ones(size(H,1),length(basis.terms))
    for i in 2:size(R,2)
        R[:,i]=R[:,basis.parents[i]].*H[:,basis.variables[i]]
    end
    return R
end
function polynomial_rollout(f,U,n)
    h=vcat([collect(U[end-j,:]) for j in 0:f.lags-1]...)
    result=fill(NaN,n,3)
    for i in 1:n
        next=f(h)
        all(isfinite,next) || break
        result[i,:]=next
        h=[next;h[1:end-3]]
    end
    return result
end
function ngrc_fixedpoints(f,starts,bounds)
    roots=DataFrame(x=Float64[],y=Float64[],z=Float64[],residual=Float64[],
        in_data_bounds=Bool[],spectral_radius=Float64[],growth_rate=Float64[],
        stability=String[],collapsed_spectral_radius=Float64[],collapse_disagreement=Float64[])
    g(u)=fixed_map(f,u)-u
    Jg(u)=last(polynomial_value_jacobian(f.fixed_basis,f.fixed_coefficients,u))-I
    for start in starts
        u=fixedpoint_newton(g,start;jacobian=Jg)
        all(isfinite,u) && norm(g(u),Inf)<=1e-8 || continue
        residual=norm(f(repeat(u,f.lags))-u,Inf)
        isfinite(residual) && residual<=1e-8 || continue
        all(norm(u-collect(r[[:x,:y,:z]]),Inf)>1e-5 for r in eachrow(roots)) || continue
        J=delayed_jacobian(f,u)
        all(isfinite,J) || continue
        rho=maximum(abs.(eigvals(J)))
        growth=log(rho)
        label=growth < -1e-6 ? "stable" : growth > 1e-6 ? "unstable" : "inconclusive"
        collapsed=last(polynomial_value_jacobian(f.fixed_basis,f.fixed_coefficients,u))
        inside=all(lo<=v<=hi for (v,(lo,hi)) in zip(u,bounds))
        push!(roots,(u...,residual,inside,rho,growth,label,maximum(abs.(eigvals(collapsed))),norm(f(repeat(u,f.lags))-fixed_map(f,u),Inf)))
    end
    return sort!(roots,[:x,:y,:z])
end

function run_ngrc_complexity()
    original_path=joinpath(@__DIR__,"02_NGRC.jl")
    original_hash=bytes2hex(sha256(read(original_path)))
    @load joinpath(@__DIR__,"../data.jld2") iso3
    countries=population_data(iso3;data_config...)
    vars=data_config.vars
    states=vcat([countries[k] for k in sort(collect(keys(countries)))]...)
    bounds=[extrema(states[!,v]) for v in vars]
    Random.seed!(0)
    starts=[vec([collect(p) for p in Iterators.product([range(lo,hi,length=7) for (lo,hi) in bounds]...)]);
        [collect(states[i,vars]) for i in unique(round.(Int,range(1,nrow(states),length=200)))];
        [[lo+rand()*(hi-lo) for (lo,hi) in bounds] for _ in 1:128]]
    structures=[(lags=k,degree=p,window="common",start_index=4) for k in 1:4 for p in (k==1 ? [1,2,3,4,5,6,8] : k==2 ? collect(1:6) : k==3 ? collect(1:5) : collect(1:4))]
    pushfirst!(structures,(lags=2,degree=2,window="original",start_index=2))
    lambdas=10.0 .^ [-12,-10,-8,-6,-4,-2,0,2,4]
    folder=joinpath(@__DIR__,"../archive/NGRC_complexity")
    mkpath(folder)
    configuration=(data_config=data_config,structures=structures,lambdas=lambdas,starts=starts,
        input_scaling="none",basis="all monomials degree 0:p",original_hash=original_hash)
    path=joinpath(folder,"configuration.jld2")
    if isfile(path)
        @assert isequal(JLD2.load(path,"configuration"),configuration)
    else
        JLD2.@save path configuration countries bounds
    end
    summaries,root_tables,country_tables=DataFrame[],DataFrame[],DataFrame[]
    prior_roots=Dict{Tuple{String,Int,Int},Vector{Vector{Float64}}}()
    original_threads=BLAS.get_num_threads()
    BLAS.set_num_threads(4)
    try
        for structure in structures
            k,p,window,start_index=structure.lags,structure.degree,structure.window,structure.start_index
            @info "Preparing NGRC structure" window lags=k degree=p features=binomial(3*k+p,p)
            H,Y=ngrc_samples(countries,k;start_index,vars)
            basis=polynomial_basis(3*k,p)
            R=ngrc_design(basis,H)
            @assert size(R,2)==binomial(3*k+p,p)
            F=svd(R)
            projected=F.U'*Y
            rank=count(>(maximum(size(R))*eps(Float64)*first(F.S)),F.S)
            condition=rank<size(R,2) ? Inf : first(F.S)/last(F.S)
            for lambda in lambdas
                identifier="$(window)_K$(k)_P$(p)_E$(round(Int,log10(lambda)))"
                path=joinpath(folder,"$identifier.jld2")
                if isfile(path)
                    saved=JLD2.load(path)
                    summary,roots,prediction_scores=saved["summary"],saved["roots"],saved["prediction_scores"]
                else
                    started=time()
                    W=F.V*((F.S ./ (F.S.^2 .+ lambda)).*projected)
                    f=polynomial_ngrc(k,p,basis,W,lambda)
                    evaluation_disagreement=0.0
                    for j in (1,div(size(H,1),2),size(H,1))
                        evaluation_disagreement=max(evaluation_disagreement,norm(f(vec(H[j,:]))-vec((R[j:j,:]*W)'),Inf))
                    end
                    extra=get(prior_roots,(window,k,p),Vector{Float64}[])
                    roots=ngrc_fixedpoints(f,[starts;extra],bounds)
                    insertcols!(roots,1,:window=>window,:lags=>k,:degree=>p,:lambda=>lambda)
                    roots.nearest_observed_distance=[minimum(norm(collect(r[vars])-collect(q[vars])) for q in eachrow(states)) for r in eachrow(roots)]
                    roots.relaxation_years=[r.stability=="stable" ? -1/r.growth_rate : NaN for r in eachrow(roots)]
                    prediction_scores=DataFrame(ISO3=String[],window=String[],lags=Int[],degree=Int[],lambda=Float64[],rmse=Float64[],failed=Bool[],persistence_rmse=Float64[])
                    for country in sort(collect(keys(countries)))
                        v=countries[country]
                        U=Matrix(v[:,vars])
                        actual=U[start_index+1:end,:]
                        pred=polynomial_rollout(f,U[1:start_index,:],size(actual,1))
                        failed=!all(isfinite,pred) || !isfinite(sum(abs2,pred-actual))
                        persistence=repeat(U[start_index:start_index,:],size(actual,1),1)
                        push!(prediction_scores,(string(country),window,k,p,lambda,failed ? Inf : sqrt(mean(abs2,pred-actual)),failed,sqrt(mean(abs2,persistence-actual))))
                    end
                    stable=roots.stability .== "stable"
                    summary=DataFrame(window=window,lags=k,degree=p,lambda=lambda,features=size(R,2),samples=size(R,1),horizon=21-start_index,
                        training_rmse=sqrt(mean(abs2,R*W-Y)),rollout_median_rmse=median(prediction_scores.rmse),rollout_mean_rmse=mean(prediction_scores.rmse),
                        rollout_failures=count(prediction_scores.failed),beats_persistence_fraction=mean(prediction_scores.rmse .< prediction_scores.persistence_rmse),
                        roots=nrow(roots),stable=count(stable),stable_in_bounds=count(stable .& roots.in_data_bounds),stable_outside=count(stable .& .!roots.in_data_bounds),
                        numerical_rank=rank,condition_number=condition,effective_df=sum(F.S.^2 ./ (F.S.^2 .+ lambda)),coefficient_norm=norm(W),
                        evaluation_disagreement=evaluation_disagreement,elapsed_seconds=time()-started)
                    JLD2.@save path f roots summary prediction_scores
                end
                prior_roots[(window,k,p)]=[collect(r[vars]) for r in eachrow(roots)]
                push!(summaries,summary);push!(root_tables,roots);push!(country_tables,prediction_scores)
                CSV.write(joinpath(folder,"summary.csv"),vcat(summaries...;cols=:union))
                CSV.write(joinpath(folder,"fixedpoints.csv"),vcat(root_tables...))
                CSV.write(joinpath(folder,"country_scores.csv"),vcat(country_tables...))
                @info "NGRC experiment completed" window lags=k degree=p lambda stable_inside=summary.stable_in_bounds[1] rmse=summary.rollout_median_rmse[1] failures=summary.rollout_failures[1]
            end
        end
    finally
        BLAS.set_num_threads(original_threads)
    end
    @assert bytes2hex(sha256(read(original_path)))==original_hash
    roots=vcat(root_tables...)
    CSV.write(joinpath(folder,"stable_fixedpoints.csv"),roots[(roots.stability .== "stable") .& roots.in_data_bounds,:])
    CSV.write(joinpath(folder,"stable_fixedpoints_all.csv"),roots[roots.stability .== "stable",:])
end

run_ngrc_complexity()
