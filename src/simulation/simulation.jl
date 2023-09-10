function initializer(POPULATION, y)
    Npop = sum(POPULATION[:,y])
    gdf = groupby(POPULATION, :location)
    
    location_ = zeros(Int8, Npop)
    for cnpop = cumsum(combine(gdf, y => sum => :npop).npop)
        location_[1:cnpop] .+= 1
    end
    location_ .-= 18
    
    gender_ = Bool[]
    for k ∈ 1:17
        male, female = combine(groupby(gdf[k], :gender), y => sum => :npop).npop
        append!(gender_, zeros(Bool, male))
        append!(gender_, ones(Bool, female))
    end
        
    age_ = Int8[]
    for (age, npop) ∈ zip(mod.((1:3400) .- 1, 100), POPULATION[:,y])
        append!(age_, repeat([age], npop))
    end

    return location_, gender_, age_
end

function validizer(p)
    if isnan(p) return rand() end
    return min(1.0, max(0.0, p))
end

function simulation(seed, tbgn, tend)
    location_, gender_, age_ = initializer(POPULATION, "y$tbgn")
    traj = deepcopy(POPULATION)
    dead = deepcopy(MORTALITY)
    mgrn = deepcopy(MOBILITY)

    pop5 = sum(reshape(POPULATION[:,"y$tbgn"], 5, :), dims = 1)
    pop5 = reshape(pop5, 20, :)
    pop5[17,:] = sum(pop5[17:end,:], dims = 1)  # 가장 아래에 있는 4행은 80세 이상
    pop5 = pop5[1:17, :] # mobility랑 칸수를 맞추기 위해 제거, 34열 = 2성별 * 17시도
    pop5 = reshape(pop5, 34, :) # 성별이 먼저 나오기 때문에 같이 복제되어야함
    pop5 = vec(repeat(pop5, 17)) # 17개 전출지별로 복제
    argyear = findfirst(names(FERTILITY) .== string(tbgn))

    tensor_population = reshape(POPULATION[:,"y$tbgn"], 100, 2, 17) # 연령(100) * 성별(2) * 시도(17)
    tensor_mortality = reshape(MORTALITY[:,"y$tbgn"], 100, 2, 17) ./ tensor_population
    tensor_mobility = reshape(MOBILITY[:,"y$tbgn"] ./ pop5, 17, 2, 17, 17) # 연령(17) * 성별(2) * 전입(17) * 전출(17)
    tensor_fertility = parse.(Float64, Matrix(FERTILITY[Not(1), argyear:(argyear + 6)])) ./ 1000

    Random.seed!(seed)

for t in tbgn:tend
    # 죽음 시작
    print('|')
    bit_age_ = Dict([a => (age_ .== a) for a ∈ 0:100])
    bit_gender_ = Dict([gen => (gender_ .== gen) for gen ∈ [false, true]])
    bit_location_ = Dict([loc => (location_ .== loc) for loc ∈ -(1:17)])
    population = Int64[]
    death = Int64[]
    for loc ∈ -(1:17), gen ∈ [false, true]
        bit_double = bit_location_[loc] .&& bit_gender_[gen]
        for a ∈ 0:99 # 반복문의 순서 자체는 기록을 위해서 바뀌어선 안 됨
            μ = validizer(tensor_mortality[a+1, gen+1, -loc])

            bit_triple = (bit_age_[a] .&& bit_double)
            pidx = findall(bit_triple)
            push!(population, length(pidx))
            died = rand(Bernoulli(μ), length(pidx))
            push!(death, count(died))
            location_[pidx[died]] .= -18
        end
    end

    traj[!, string('y', t)] = population
    dead[!, string('y', t)] = death
    println(Dates.now())
    CSV.write("G:/rslt $(lpad(seed, 4, '0')).csv", traj, encoding = "UTF-8", bom = true)
    CSV.write("G:/dead $(lpad(seed, 4, '0')).csv", dead, encoding = "UTF-8", bom = true)

    location_[age_ .≥ 100] .= -18
    bit_location_ = Dict([loc => (location_ .== loc) for loc ∈ -(1:17)])
    # 죽음 끝

    # 출산 시작
    bit_age5_ = Dict([a => (a .≤ age_ .< (a + 5)) for a ∈ 0:5:80])
    bit_age5_[80] = (80 .≤ age_)
    birth_location = zeros(Int64, 17)
    flow = Int64[]
    for loc ∈ -(1:17), gen ∈ [false, true]
        bit_double = bit_location_[loc] .&& bit_gender_[gen]
        for a ∈ (0:5:80)
            aidx = (a ÷ 5) + 1
            bit_triple = bit_age5_[a] .&& bit_double
            pidx = shuffle(findall(bit_triple))
            
            # 이동 시작
            npop = length(pidx)
            σ = tensor_mobility[aidx, gen+1, :, -loc]
            if (( 1000 ≤ seed ≤ 1017 ) && (t ≤ 2030)) σ[seed - 1000] = 2σ[seed - 1000] end
            moved = trunc.(Int64, σ .* npop)
            append!(flow, moved)
            moved = reverse(cumsum(moved))
            for to in 1:17
                location_[pidx[1:moved[to]]] .= (to-18) # 인구수 이상으로는 이동 불가
                # location_[pidx[1:min(npop, moved[to])]] .= (to-18) # 인구수 이상으로는 이동 불가
                # 자연스럽게 만약 다 못간다면 서울이 가장 우선시 됨
            end
            # 이동 끝

            if gen .&& (15 ≤ a ≤ 45)
                β = validizer(tensor_fertility[-loc, aidx - 3])
                birth_location[-loc] += sum(rand(Bernoulli(β), npop))
            end
        end
    end
    flow = vcat([vec(reshape(flow, 17, 34, :)[:,:,k]') for k ∈ 1:17]...) # 원본데이터와 일치하도록 수정
    mgrn[!, string('y', t)] = flow
    CSV.write("G:/mgrn $(lpad(seed, 4, '0')).csv", mgrn, encoding = "UTF-8", bom = true)
    age_[location_ .!= (-18)] .+= 1
    append!(age_, zeros(Int8, sum(birth_location)))
    append!(gender_, rand(Bernoulli(female_ratio), sum(birth_location))) # 성비는 남:여 = 105:100
    append!(location_, vcat([repeat([loc], birth_location[-loc]) for loc ∈ -(1:17)]...))
    # 출산 끝
end

end