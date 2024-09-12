using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays


function rand_pm1_matrix(M, N)
     return rand(Bool, M, N) .* 2 .- 1.0
end
 
 


M = 100
S1 = rand(M, M)
S2 = rand(M, M)

S1 = S1/norm(S1)
S2 = S2/norm(S2)

@show cond(S1)
@show cond(S2)

@assert norm(S1) ≈ norm(S2) # Vergleichbarkeit...

b = rand(M)



# exact solution

x1 = S1 \ b 
x2 = S2 \ b 
""

##

##


N = 500
Δx1_list = Vector{Float64}(undef,N)
Δx2_list = Vector{Float64}(undef,N)
max_rel_integral_error = 1.0e-4

for i in 1:N
     
     ΔS1 = S1 .* rand_pm1_matrix(M, M) .* rand(M, M) * max_rel_integral_error
     ΔS2 = S2 .* rand_pm1_matrix(M, M) .* rand(M, M) * max_rel_integral_error
     Δb = b .* rand_pm1_matrix(M, 1) .* rand(M, 1) * max_rel_integral_error

     x1e = (S1 + ΔS1) \ (b + Δb) 
     #norm((S1 + ΔS1)*x1e - (b + Δb))

     x2e = (S2 + ΔS2) \ (b + Δb) 
     #norm((S2 + ΔS2)*x2e - (b + Δb))

     Δx1_list[i] = norm(x1e - x1)/norm(x1)
     Δx2_list[i] = norm(x2e - x2)/norm(x2)
end


@show sum(Δx1_list)/N
@show sum(Δx2_list)/N
@show maximum(Δx1_list)
@show maximum(Δx2_list)

""
