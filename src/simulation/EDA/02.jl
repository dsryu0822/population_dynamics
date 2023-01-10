include("basic.jl")

t = 2021
mtrx_mgrn_ = Dict()
vctr_mgrn = mgrn[:, "y$t"]
mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
mtrx_mgrn .*= (1 .- I(17))
