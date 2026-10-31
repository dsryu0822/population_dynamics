using Flux, NearestNeighbors, SHA, SparseArrays

function model_periods(iso3, ids; train_years = 2000:2019, test_years = 2020:2024, vars = [:x,:y,:z])
    @assert length(ids)==161 && length(unique(ids))==161
    @assert length(train_years)>=3 && !isempty(test_years) && last(train_years)+1==first(test_years)
    train, test = Dict{String,DataFrame}(), Dict{String,DataFrame}()
    for id in ids
        raw = sort(iso3[id], :t)
        train[id] = raw[in.(raw.t,Ref(train_years)), [:t;vars]]
        test[id] = raw[in.(raw.t,Ref(test_years)), [:t;vars]]
        @assert train[id].t==collect(train_years) && test[id].t==collect(test_years) "Missing/duplicate year: $id"
        @assert all(isfinite,Matrix(train[id][:,vars])) && all(isfinite,Matrix(test[id][:,vars]))
    end
    return train, test
end

function model14_pairs(train; vars = [:x,:y,:z])
    return map_data(train;vars)
end

if !isdefined(@__MODULE__,:RCmap14)
    struct RCmap14{T}
        payload::T
        r::Vector{Float64}
    end
end
function rc_state14(payload,r,u)
    m=payload.model
    return (1-m.α).*r.+m.α.*tanh.(m.A*r+m.Win*((u.-payload.input_mean)./payload.input_scale).+m.ξ)
end
function (f::RCmap14)(u::AbstractVector)
    f.r .= rc_state14(f.payload,f.r,u)
    return vec(f.payload.model.Wout*f.r+f.payload.model.c)
end
reset_rc14!(f::RCmap14)=(fill!(f.r,0.0);f)
Base.show(io::IO,f::RCmap14)=print(io,"RC($(f.payload.model.N) nodes, α=$(f.payload.model.α))")

if !isdefined(@__MODULE__,:EnsembleMap14)
    struct EnsembleMap14{T}
        models::T
    end
end

# Copy the reservoir state so ensemble forecasts do not alter the original RC.
function ensemble(models...)
    isempty(models)&&throw(ArgumentError("At least one model is required."))
    any(f->f isa EnsembleMap14 && rc_memory14(f)!==nothing,models)&&
        throw(ArgumentError("Pass RC directly to ensemble rather than nesting a stateful ensemble."))
    count(f->f isa RCmap14,models)<=1||throw(ArgumentError("One RC reservoir per ensemble is supported."))
    members=map(f->f isa RCmap14 ? RCmap14(f.payload,copy(f.r)) : f,models)
    return EnsembleMap14(members)
end
function mapvalue14(f,u)
    value=vec(Float64.(collect(applicable(f,u) ? f(u) : f(u...))))
    length(value)==3||throw(DimensionMismatch("The model must return three coordinates."))
    return value
end
(f::EnsembleMap14)(u::AbstractVector)=sum(mapvalue14(model,u) for model in f.models)/length(f.models)
(f::EnsembleMap14)(x::Real,y::Real,z::Real)=f([x,y,z])
function rc_memory14(f)
    f isa RCmap14 && return f
    if f isa EnsembleMap14
        i=findfirst(model->model isa RCmap14,f.models)
        return i===nothing ? nothing : f.models[i]
    end
    return nothing
end

function batch_mapvalue14(f,U)
    f isa Flux.Chain && return Matrix(f(U))
    f isa STLSQresult && return permutedims(reshape(f(permutedims(U)),size(U,2),3))
    return hcat([mapvalue14(f,u) for u in eachcol(U)]...)
end

# Each column has its own reservoir state; countries never share memory.
function basin_step14(f,U,R)
    memory=rc_memory14(f)
    nextR=if memory===nothing
        R
    else
        p=memory.payload;m=p.model
        (1-m.α).*R.+m.α.*tanh.(m.A*R+m.Win*((U.-p.input_mean)./p.input_scale).+m.ξ)
    end
    members=f isa EnsembleMap14 ? f.models : (f,)
    nextU=zeros(3,size(U,2))
    for model in members
        nextU .+= model isa RCmap14 ? model.payload.model.Wout*nextR.+model.payload.model.c : batch_mapvalue14(model,U)
    end
    return nextU./length(members),nextR
end

"""
    country_basins(f, data, ids, points; year=2024, warmup_start=1950,
                   maxsteps=10000, tol=1e-5, settle=20, vars=[:x,:y,:z])

Classify countries by their limiting stable fixed point, numbered by increasing z.
Use the observed state in `year`; RC memory is warmed with actual observations
from warmup_start through year-1. The first recursive step then consumes `year`
once. Future observations are never used. Require both observed and reservoir
states to stay within tol of the same point for settle consecutive steps.
Unresolved/nonfinite trajectories retain a missing point label.
"""
function country_basins(f,data,ids,points;year=2024,warmup_start=1950,
    maxsteps=10000,tol=1e-5,settle=20,vars=[:x,:y,:z])
    maxsteps>=settle>=1 && tol>0 && warmup_start<=year||throw(ArgumentError("Invalid convergence settings."))
    ids=sort(string.(collect(ids)))
    !isempty(ids) && length(unique(ids))==length(ids)||throw(ArgumentError("Country IDs must be nonempty and unique."))
    refs=sort([Float64.(collect(p)) for p in points];by=last)
    !isempty(refs) && all(p->length(p)==3 && all(isfinite,p),refs)||throw(ArgumentError("Supply three-dimensional stable points."))
    memory=rc_memory14(f);N=memory===nothing ? 0 : memory.payload.model.N
    targets=zeros(N,length(refs))
    for (j,u) in enumerate(refs)
        if memory===nothing
            norm(mapvalue14(f,u)-u,Inf)<=1e-7||throw(ArgumentError("Reference point $j is not a fixed point."))
        else
            targets[:,j]=rc_equilibrium14(memory.payload,u)
            state=[u;targets[:,j]]
            norm(fp_joint_step14(f,state)-state,Inf)<=1e-7||throw(ArgumentError("Reference point $j is not a joint RC fixed point."))
        end
    end
    count=length(ids);initial=zeros(3,count);R=zeros(N,count)
    for (i,id) in enumerate(ids)
        raw=sort(data[id],:t)
        history=raw[(warmup_start .<= raw.t) .& (raw.t .<= year),[:t;vars]]
        @assert history.t==collect(warmup_start:year) "Missing/duplicate history year: $id"
        H=Matrix(history[:,vars]);@assert all(isfinite,H)
        initial[:,i]=H[end,:]
        if memory!==nothing
            r=zeros(N)
            for j in 1:size(H,1)-1;r=rc_state14(memory.payload,r,H[j,:]);end
            R[:,i]=r
        end
    end
    final=copy(initial);U=copy(initial);active=collect(1:count)
    labels=zeros(Int,count);nearest=zeros(Int,count);streak=zeros(Int,count);lastpoint=zeros(Int,count)
    steps=zeros(Int,count);distance=fill(Inf,count);memory_distance=fill(Inf,count)
    status=fill("unresolved",count)
    for step in 1:maxsteps
        nextU,nextR=basin_step14(f,U,R);keep=Int[]
        for (column,i) in enumerate(active)
            final[:,i]=nextU[:,column];steps[i]=step
            if !all(isfinite,nextU[:,column]) || !all(isfinite,nextR[:,column])
                status[i]="nonfinite";distance[i]=Inf;memory_distance[i]=Inf;nearest[i]=0;continue
            end
            distances=[norm(nextU[:,column]-u) for u in refs]
            j=argmin(distances);nearest[i]=j;distance[i]=distances[j]
            memory_distance[i]=N==0 ? 0.0 : norm(nextR[:,column]-targets[:,j],Inf)
            close=distance[i]<=tol && memory_distance[i]<=tol
            streak[i]=close ? (lastpoint[i]==j ? streak[i]+1 : 1) : 0
            lastpoint[i]=close ? j : 0
            if streak[i]>=settle
                labels[i]=j;status[i]="converged"
            else
                push!(keep,column)
            end
        end
        active=active[keep];U=nextU[:,keep];R=nextR[:,keep]
        isempty(active)&&break
        step%500==0 && @info "Country convergence" step remaining=length(active)
    end
    result=DataFrame(ISO3=ids,start_year=fill(year,count),point=Union{Missing,Int}[j==0 ? missing : j for j in labels],
        point_name=[j==0 ? "미수렴" : "제$(j)점" for j in labels],status=status,steps=steps,
        nearest_point=Union{Missing,Int}[j==0 ? missing : j for j in nearest],distance=distance,reservoir_distance=memory_distance)
    insertcols!(result,3,(Symbol(string(v)*string(year))=>vec(initial[j,:]) for (j,v) in enumerate(vars))...)
    return hcat(result,DataFrame(permutedims(final),Symbol.(string.(vars).*"_end")))
end

function fit_rc14(train,X,Y,spec;vars=[:x,:y,:z])
    N,α=spec.n_nodes,spec.alpha
    β,ρ,ξ=1e-2,0.9,1.0 # Fixed settings; only N and α are searched.
    rng=MersenneTwister(0)
    A=sprand(rng,N,N,2/N);A.nzval .= 2rand(rng,nnz(A)).-1
    A .*= ρ/maximum(abs.(eigvals(Matrix(A))))
    Win=zeros(N,3)
    for i in 1:N;Win[i,rand(rng,1:3)]=2rand(rng)-1;end
    μ,σ=vec(mean(X,dims=2)),max.(vec(std(X,dims=2)),1e-8)
    m=Reservoir(3,3,N,α,β,ρ,A,ξ,Win,zeros(3,N),zeros(3),zeros(N))
    initial=(model=m,input_mean=μ,input_scale=σ)
    R=zeros(N,size(X,2));column=0
    for id in sort(collect(keys(train)))
        U=Matrix(train[id][:,vars]);r=zeros(N)
        for i in 1:size(U,1)-1
            r=rc_state14(initial,r,U[i,:]);column+=1;R[:,column]=r
        end
    end
    @assert column==size(Y,2)
    μR,μY=mean(R,dims=2),mean(Y,dims=2)
    Rc,Yc=R.-μR,Y.-μY
    Wout=Matrix((Symmetric(Rc*Rc'+β*I)\(Rc*Yc'))')
    c=vec(μY-Wout*μR)
    model=Reservoir(3,3,N,α,β,ρ,A,ξ,Win,Wout,c,zeros(N))
    return (model=model,input_mean=μ,input_scale=σ)
end

# Reservoir equilibrium under a constant observed state, used only for root search.
function rc_equilibrium14(payload,u)
    m=payload.model
    all(isfinite,u)||return fill(NaN,m.N)
    drive=m.Win*((u.-payload.input_mean)./payload.input_scale).+m.ξ
    r=tanh.(drive)
    for _ in 1:50
        t=tanh.(m.A*r+drive);g=t-r
        norm(g,Inf)<=1e-12 && return r
        step=try (I-Diagonal(1 .-t.^2)*m.A)\g catch; return fill(NaN,m.N);end
        accepted=false
        for j in 0:20
            candidate=r+0.5^j*step
            if norm(tanh.(m.A*candidate+drive)-candidate,Inf)<norm(g,Inf)
                r=candidate;accepted=true;break
            end
        end
        accepted||break
    end
    return r
end

function rc_stability14(payload,u,r)
    m=payload.model
    t=tanh.(m.A*r+m.Win*((u.-payload.input_mean)./payload.input_scale).+m.ξ)
    D=Diagonal(1 .-t.^2)
    B=m.α.*D*m.Win*Diagonal(1 ./payload.input_scale)
    C=Matrix((1-m.α)*I+m.α.*D*m.A)
    J=[m.Wout*B m.Wout*C; B C]
    step(v)=begin
        next=rc_state14(payload,v[4:end],v[1:3])
        [m.Wout*next+m.c;next]
    end
    values=eigvals(J);growth=log(maximum(abs.(values)))
    label=growth < -1e-6 ? "stable" : growth > 1e-6 ? "unstable" : "inconclusive"
    difference=maximum(abs.(J-numerical_jacobian(step,[u;r])))
    return values,growth,label,difference
end

function fit_sigmoid14(X,Y,spec)
    μX,σX=vec(mean(X,dims=2)),max.(vec(std(X,dims=2)),1e-8)
    μY,σY=vec(mean(Y,dims=2)),max.(vec(std(Y,dims=2)),1e-8)
    Xn,Yn=Float32.((X.-μX)./σX),Float32.((Y.-μY)./σY)
    Random.seed!(spec.seed)
    model=Flux.Chain(Flux.Dense(size(X,1)=>spec.width,Flux.sigmoid),
        [Flux.Dense(spec.width=>spec.width,Flux.sigmoid) for _ in 2:spec.hidden_layers]...,Flux.Dense(spec.width=>size(Y,1)))
    optimizer=Flux.setup(Flux.Adam(spec.rate),model)
    best=deepcopy(model);losses=[mlp_loss(model,Xn,Yn)];best_loss=first(losses)
    threads=BLAS.get_num_threads();BLAS.set_num_threads(1)
    try
        for epoch in 1:spec.epochs
            Flux.train!(mlp_loss,model,[(Xn,Yn)],optimizer)
            loss=mlp_loss(model,Xn,Yn)
            if loss<best_loss;best,best_loss=deepcopy(model),loss;end
            push!(losses,loss)
            epoch%1000==0 && @info "Sigmoid MLP" ID=spec.ID epoch best_loss
        end
    finally
        BLAS.set_num_threads(threads)
    end
    return (model=Flux.f64(best),input_mean=μX,input_scale=σX,output_mean=μY,output_scale=σY,losses=losses)
end

function fit_model14(kind,spec,train,X,Y;vars=[:x,:y,:z],tree=nothing)
    if kind==:RC
        payload=fit_rc14(train,X,Y,spec;vars)
        return restore_model14(kind,payload),payload
    elseif kind==:MLP
        run=fit_sigmoid14(X,Y,spec)
        payload=(run=run,)
        return restore_model14(kind,payload),payload
    elseif kind==:RBF
        balls=treeregion.(ball_leaves(tree,spec.n_basis))
        centers=stack(b.center for b in balls);widths=spec.width_scale.*max.([b.r for b in balls],1e-6)
        sysms=("next_".*string.(vars),string.(vars))
        recipe=cook_rbf(sysms,centers,widths;poly=0:0)[2:end,:]
        recipe.index=1:nrow(recipe)
        @assert nrow(recipe)==spec.n_basis && all(g->g isa GaussianRBF,recipe.func)
        data=hcat(DataFrame(X',last(sysms)),DataFrame(Y',first(sysms)))
        f=SINDy(data,sysms,recipe) # Default λ=0: least squares without thresholding.
        payload=(model=f,centers=centers,widths=widths)
        return restore_model14(kind,payload),payload
    end
    error("Unknown model: $kind")
end

function restore_model14(kind,payload)
    kind==:RC && return RCmap14(payload,zeros(payload.model.N))
    return kind==:MLP ? mlp_map(payload.run) : payload.model
end
function best_model14(kind,table,output)
    valid=filter(:status=>==("ok"),table)
    isempty(valid)&&error("No valid $kind models")
    id=valid.ID[argmin(valid.test_rmse)]
    return restore_model14(kind,JLD2.load(joinpath(output,"weights",id*".jld2"))["payload"])
end

function value_jacobian14(kind,payload)
    if kind==:RC
        m=payload.model;V=m.Win*Diagonal(1 ./payload.input_scale)
        previous=fill(NaN,3);cached=(zeros(3),zeros(3,3))
        return u->begin
            u==previous && return cached
            r=rc_equilibrium14(payload,u)
            if !all(isfinite,r) || norm(tanh.(m.A*r+V*u.-m.Win*(payload.input_mean./payload.input_scale).+m.ξ)-r,Inf)>1e-9
                return fill(NaN,3),fill(NaN,3,3)
            end
            D=Diagonal(1 .-r.^2)
            J=try m.Wout*((I-D*m.A)\(D*V)) catch; fill(NaN,3,3);end
            previous=copy(u);cached=(vec(m.Wout*r+m.c),J)
            cached
        end
    elseif kind==:MLP
        run=payload.run
        return u->begin
            value=(u.-run.input_mean)./run.input_scale;J=Matrix(Diagonal(1 ./ run.input_scale))
            layers=run.model.layers
            for layer in layers[1:end-1]
                value=Flux.sigmoid.(layer.weight*value+layer.bias)
                J=(value.*(1 .-value)).*(layer.weight*J)
            end
            value=layers[end].weight*value+layers[end].bias
            J=Diagonal(run.output_scale)*layers[end].weight*J
            value.*run.output_scale.+run.output_mean,J
        end
    else
        W=Matrix(payload.model.matrix);centers=payload.centers;widths=payload.widths
        return u->begin
            value=zeros(3);J=zeros(3,3)
            for k in axes(centers,2)
                difference=u-centers[:,k];phi=exp(-sum(abs2,difference)/(2widths[k]^2))
                iszero(phi)&&continue
                for i in 1:3
                    weight=W[k,i]*phi;value[i]+=weight
                    for j in 1:3;J[i,j]-=weight*difference[j]/widths[k]^2;end
                end
            end
            value,J
        end
    end
end

function fixedpoints14(kind,f,payload,states;vars=[:x,:y,:z])
    bounds=[extrema(states[!,v]) for v in vars];rng=MersenneTwister(0)
    starts=[vec([collect(p) for p in Iterators.product([range(lo,hi,length=7) for (lo,hi) in bounds]...)]);
        [collect(states[i,vars]) for i in unique(round.(Int,range(1,nrow(states),length=200)))];
        [collect(center).+collect(offset) for center in ((5.3,5.9,5.4),(7.,8.,7.)) for offset in Iterators.product(fill(-0.25:0.25:0.25,3)...)];
        [[lo+rand(rng)*(hi-lo) for (lo,hi) in bounds] for _ in 1:128]]
    vj=value_jacobian14(kind,payload)
    fixed(u)=kind==:RC ? vec(payload.model.Wout*rc_equilibrium14(payload,u)+payload.model.c) : vec(f(u))
    g(u)=first(vj(u))-u;Jg(u)=last(vj(u))-I
    result=DataFrame(x=Float64[],y=Float64[],z=Float64[],residual=Float64[],in_data_bounds=Bool[],spectral_radius=Float64[],growth_rate=Float64[],stability=String[],jacobian_difference=Float64[],state_dimension=Int[],reservoir=Vector{Float64}[])
    for start in starts
        u=fixedpoint_newton(g,start;jacobian=Jg)
        all(isfinite,u) && norm(g(u),Inf)<=1e-8 || continue
        residual=norm(fixed(u)-u,Inf);isfinite(residual)&&residual<=1e-8||continue
        all(norm(u-collect(r[vars]),Inf)>1e-5 for r in eachrow(result))||continue
        r=kind==:RC ? rc_equilibrium14(payload,u) : Float64[]
        if kind==:RC
            residual=max(residual,norm(rc_state14(payload,r,u)-r,Inf))
            residual<=1e-8||continue
            values,growth,label,difference=rc_stability14(payload,u,r)
        else
            values,growth,label=stability(f,u;discrete=true,lags=1)
            difference=maximum(abs.(last(vj(u))-numerical_jacobian(fixed,u)))
        end
        inside=all(lo<=v<=hi for (v,(lo,hi)) in zip(u,bounds))
        push!(result,(u...,residual,inside,maximum(abs.(values)),growth,label,difference,length(values),r))
    end
    return sort!(result,vars)
end

function rollout14(kind,f,U,n)
    u=vec(U[end,:]);pred=fill(NaN,n,3)
    memory=rc_memory14(f)
    if memory!==nothing
        reset_rc14!(memory)
        for i in 1:size(U,1)-1;memory(vec(U[i,:]));end
    end
    for i in 1:n
        u=vec(f(u));pred[i,:]=u
        all(isfinite,u)||break
    end
    return pred
end

# Pure maps use their ordinary Jacobian; RC uses the constant-input equilibrium.
function fp_value_jacobian14(f)
    f isa RCmap14 && return value_jacobian14(:RC,f.payload)
    if f isa EnsembleMap14
        derivatives=map(fp_value_jacobian14,f.models)
        return u->begin
            values=map(vj->vj(u),derivatives)
            sum(first.(values))/length(values),sum(last.(values))/length(values)
        end
    end
    value=u->mapvalue14(f,u)
    return u->(value(u),numerical_jacobian(value,u))
end

# One autonomous ensemble step, with reservoir state made explicit for analysis.
function fp_joint_step14(f,v)
    memory=rc_memory14(f);p=memory.payload
    u=v[1:3];r=rc_state14(p,v[4:end],u)
    members=f isa EnsembleMap14 ? f.models : (f,)
    values=[model isa RCmap14 ? vec(p.model.Wout*r+p.model.c) : mapvalue14(model,u) for model in members]
    return [sum(values)/length(values);r]
end

"""
    findfp(f; bounds=((0.,10.),(0.,10.),(0.,10.)), ngrid=7, starts=nothing,
           discrete=true, atol=1e-8, merge_tol=1e-5, stability_tol=1e-6, maxiters=100)

Return discovered stable fixed points as Vector{Vector{Float64}}.
For an annual map solve f(u)=u; discrete=false instead solves a vector field f(u)=0.
Pure callables may accept a vector or three scalars. Combine a stateful RC with
ensemble(...) so its equilibrium and full-state stability can be handled explicitly.
Search bounds generate starting points, rather than excluding roots outside them.
An RC ensemble is discrete-time, and its reservoir memory is preserved by this search.
"""
function findfp(f;bounds=((0.,10.),(0.,10.),(0.,10.)),ngrid=7,starts=nothing,
    discrete=true,atol=1e-8,merge_tol=1e-5,stability_tol=1e-6,maxiters=100)
    length(bounds)==3 && all(b->length(b)==2 && all(isfinite,b) && b[1]<b[2],bounds)||
        throw(ArgumentError("bounds must contain three finite (lower,upper) pairs."))
    ngrid>=2 && atol>0 && merge_tol>0 && stability_tol>=0 && maxiters>=1||
        throw(ArgumentError("Invalid search tolerances or iteration count."))
    memory=rc_memory14(f)
    memory!==nothing && !discrete && throw(ArgumentError("RC ensembles are discrete maps; use discrete=true."))
    seeds=starts===nothing ? Iterators.product([range(lo,hi,length=ngrid) for (lo,hi) in bounds]...) : starts
    vj=fp_value_jacobian14(f)
    g(u)=first(vj(u))-(discrete ? u : zeros(3))
    Jg(u)=last(vj(u))-(discrete ? I : zeros(3,3))
    roots=Vector{Vector{Float64}}();states=Vector{Vector{Float64}}()
    for start in seeds
        length(start)==3||throw(DimensionMismatch("Each starting point must have three coordinates."))
        u=fixedpoint_newton(g,start;tol=min(1e-12,atol),maxiters,jacobian=Jg)
        all(isfinite,u) && norm(g(u),Inf)<=atol||continue
        all(norm(u-v,Inf)>merge_tol for v in roots)||continue
        state=memory===nothing ? copy(u) : [u;rc_equilibrium14(memory.payload,u)]
        if memory!==nothing
            all(isfinite,state) && norm(fp_joint_step14(f,state)-state,Inf)<=atol||continue
        end
        push!(roots,u);push!(states,state)
    end
    stable=Vector{Vector{Float64}}()
    for (u,state) in zip(roots,states)
        J=memory===nothing ? last(vj(u)) : numerical_jacobian(v->fp_joint_step14(f,v),state)
        all(isfinite,J)||continue
        values=eigvals(J)
        growth=discrete ? log(maximum(abs.(values))) : maximum(real.(values))
        growth < -stability_tol && push!(stable,u)
    end
    return sort!(stable;by=Tuple)
end

function model_scores14(kind,f,train,test;vars=[:x,:y,:z])
    totals=zeros(2);counts=zeros(Int,2);rows=DataFrame(ISO3=String[],split=String[],rmse=Float64[],n_predictions=Int[])
    for id in sort(collect(keys(train)))
        U=Matrix(train[id][:,vars]);V=Matrix(test[id][:,vars])
        if kind==:RC
            reset_rc14!(f)
            fitted=permutedims(hcat([f(vec(U[i,:])) for i in 1:size(U,1)-1]...))
        else
            fitted=kind==:MLP ? permutedims(f(permutedims(U[1:end-1,:]))) : reshape(f(U[1:end-1,:]),size(U,1)-1,3)
        end
        forecast=rollout14(kind,f,U,size(V,1))
        for (i,(actual,predicted)) in enumerate(((U[2:end,:],fitted),(V,forecast)))
            N=size(actual,1)
            @assert size(predicted)==size(actual)
            error=sum(abs2,predicted-actual);isfinite(error)||(error=Inf)
            totals[i]+=error;counts[i]+=length(actual)
            push!(rows,(id,i==1 ? "train" : "test",sqrt(error/length(actual)),N))
        end
    end
    @assert counts==length(vars).*[sum(nrow(v)-1 for v in values(train)),sum(nrow(v) for v in values(test))]
    return sqrt.(totals./counts)...,rows
end

function model_signature14(train,test,X,Y;vars=[:x,:y,:z])
    io=IOBuffer();write(io,"models14-recursive-test-v1");write(io,join(sort(collect(keys(train))),"\n"));write(io,X);write(io,Y)
    for id in sort(collect(keys(train)))
        write(io,Matrix(train[id][:,vars]));write(io,Matrix(test[id][:,vars]))
    end
    return bytes2hex(sha256(take!(io)))
end

function training_signature14(train,X,Y;vars=[:x,:y,:z])
    io=IOBuffer();write(io,"models14-rc-training-v1");write(io,join(string.(vars),"\n"));write(io,X);write(io,Y)
    for id in sort(collect(keys(train)))
        write(io,id);write(io,Float64.(train[id].t));write(io,Matrix(train[id][:,vars]))
    end
    return bytes2hex(sha256(take!(io)))
end

function classify_fixedpoints14!(points;vars=[:x,:y,:z])
    points.dist_fp5=[norm(collect(r[vars]).-[5.,5.,5.]) for r in eachrow(points)]
    points.dist_fp7=[norm(collect(r[vars]).-[7.,7.,7.]) for r in eachrow(points)]
    # Equal distances belong to fp5, so each stable point is counted exactly once.
    points.fp_group=ifelse.(points.dist_fp5 .<= points.dist_fp7,"fp5","fp7")
    return points
end

function run_model_grid14(kind,specs,train,test,states,output;filename,resume=true,vars=[:x,:y,:z])
    mkpath(output);mkpath(joinpath(output,"weights"))
    mkpath(dirname(joinpath(output,filename)))
    X,Y=model14_pairs(train;vars)
    tree=kind==:RBF ? BallTree(X;leafsize=1,reorder=false) : nothing
    signature=model_signature14(train,test,X,Y;vars)
    train_signature=training_signature14(train,X,Y;vars)
    train_years=first(values(train)).t;test_years=first(values(test)).t
    test_samples=sum(nrow(v) for v in values(test))
    table=DataFrame(ID=String[],train_rmse=Union{Missing,Float64}[],test_rmse=Union{Missing,Float64}[],fp5=Union{Missing,Int}[],fp7=Union{Missing,Int}[])
    for spec in specs
        path=joinpath(output,"weights",spec.ID*".jld2")
        @info "Model grid" ID=spec.ID kind
        try
            saved=resume && isfile(path) ? JLD2.load(path) : nothing
            reusable=saved!==nothing && saved["spec"]==spec &&
                get(saved,"train_signature",nothing)==train_signature
            fit_seconds,fp_seconds=0.0,0.0
            if reusable
                payload,points=saved["payload"],saved["points"]
                f=restore_model14(kind,payload)
            else
                fit_seconds=@elapsed begin
                    f,payload=fit_model14(kind,spec,train,X,Y;vars,tree)
                end
                @info "Training complete; finding fixed points" ID=spec.ID fit_seconds
                fp_seconds=@elapsed points=fixedpoints14(kind,f,payload,states;vars)
            end
            score_seconds=@elapsed begin
                train_rmse,test_rmse,scores=model_scores14(kind,f,train,test;vars)
            end
            @info "Model timings (seconds)" ID=spec.ID cached=reusable fit_seconds fp_seconds score_seconds
            classify_fixedpoints14!(points;vars)
            JLD2.jldsave(path;kind,spec,payload,points,train_rmse,test_rmse,scores,signature,train_signature,train_years,test_years,test_mode="recursive")
            stable=points[points.stability .== "stable",:]
            fp5=count(==("fp5"),stable.fp_group);fp7=count(==("fp7"),stable.fp_group)
            @assert fp5+fp7==nrow(stable)
            row=merge((ID=spec.ID,train_rmse=train_rmse,test_rmse=test_rmse,fp5=fp5,fp7=fp7),spec,
                (fp_total=nrow(stable),train_samples=size(X,2),test_samples=test_samples,test_mode="recursive",status=all(isfinite,(train_rmse,test_rmse)) ? "ok" : "nonfinite_prediction",error=""))
            push!(table,row;cols=:union)
        catch error
            @error "Model failed" ID=spec.ID exception=(error,catch_backtrace())
            push!(table,merge((ID=spec.ID,train_rmse=missing,test_rmse=missing,fp5=missing,fp7=missing),spec,
                (fp_total=missing,train_samples=size(X,2),test_samples=test_samples,test_mode="recursive",status="error",error=sprint(showerror,error)));cols=:union)
        end
        @assert names(table)[1:5]==["ID","train_rmse","test_rmse","fp5","fp7"]
        CSV.write(joinpath(output,filename),table)
    end
    return table
end
