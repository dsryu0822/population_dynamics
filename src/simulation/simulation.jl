function initializer(POPULATION, y)
    Npop = sum(POPULATION[:,end])
    gdf = groupby(POPULATION, :location)
    
    location_ = zeros(Int8, Npop)
    for cnpop = cumsum(combine(gdf, y => sum => :npop).npop)
        location_[1:cnpop] .+= 1
    end
    location_ .-= 18
    # location_ = CuArray(location_)
    
    gender_ = Bool[]
    for k ∈ 1:17
        male, female = combine(groupby(gdf[k], :gender), y => sum => :npop).npop
        append!(gender_, zeros(Bool, male))
        append!(gender_, ones(Bool, female))
    end
    # gender_ = CuArray(gender_)
        
    age_ = Int8[]
    for (age, npop) ∈ zip(mod.((1:3400) .- 1, 100), POPULATION[:,end])
        append!(age_, repeat([age], npop))
    end
    # age_ = CuArray(age_)

    return location_, gender_, age_
end

function simulation(seed, POPULATION)
    location_, gender_, age_ = initializer(POPULATION, :y2021)
    Random.seed!(seed)

    for t in 2021:2050

    # 죽음 시작
    print('|')
    bit_location_ = Dict([loc => (location_ .== loc) for loc ∈ -(1:18)])
    bit_age_ = Dict([a => (age_ .== a) for a ∈ 0:100])
    bit_gender_ = Dict([gen => (gender_ .== gen) for gen ∈ [false, true]])
    population = Int64[]
    for loc ∈ -(1:17), gen ∈ [false, true], a ∈ 0:99
        # println("loc: $loc, gen: $gen, a: $a")
        pidx = findall(bit_location_[loc] .&& bit_gender_[gen] .&& bit_age_[a])
        push!(population, length(pidx))

        μ = tensor_mortality[maximum(a)+1,-loc,gen+1]
        died = rand(Bernoulli(μ), length(pidx)) # rand() .< μ 랑 속도는 같음
        location_[pidx[died]] .= -18
    end
    POPULATION[!, string('y', t)] = population
    location_[age_ .≥ 100] .= -18
    bit_location_ = Dict([loc => (location_ .== loc) for loc ∈ -(1:18)])
    # 죽음 끝

    # 출산 시작
    print('.')

    birth_location = zeros(Int64, 17)
    bit_age5_ = Dict([a => (a .≤ age_ .< (a + 5)) for a ∈ 0:5:80])
    for gen ∈ [false, true], loc ∈ -(1:17), a ∈ (0:5:80)
        aidx = (a ÷ 5) + 1
        # 이동 시작
        pidx = shuffle(findall(bit_location_[loc] .&& bit_gender_[gen] .&& bit_age5_[a]))
        σ = tensor_mobility[aidx, gen+1, :, -loc]
        moved = reverse(cumsum(σ))
        for j in 1:17
            location_[pidx[1:moved[j]]] .= (j-18)
        end
        # 이동 끝

        if gen .&& (15 ≤ a ≤ 45)
            npop = length(pidx)
            β = tensor_fertility[-loc, aidx - 3]
            birth_location[-loc] += trunc(Int64, npop*β/1000)
        end
    end
    age_ .+= 1
    append!(gender_, rand(Bernoulli(female_ratio), sum(birth_location)))
    append!(age_, zeros(Int8, sum(birth_location)))
    append!(location_, vcat([repeat([loc], birth_location[-loc]) for loc ∈ -(1:17)]...))
    # 출산 끝

    CSV.write("D:/rslt $(lpad(seed, 4, '0')).csv", POPULATION, encoding = "UTF-8", bom = true)
    end
end