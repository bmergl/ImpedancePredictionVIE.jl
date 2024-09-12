using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays


function rand_pm1_matrix(M, N)
     return rand(Bool, M, N) .* 2 .- 1.0
end
 
 
 
# cond völlig verschieden, aber Fehler etwa gleich!
# S1 = [100.0 1.0; 
#      10.0 20000.0]
# S2 = [1.0 100.0; 
#      10.0 20000.0]

# cond völlig verschieden, Fehler verschieden
# S1 = [100.0 1.0; 
#      10.0 2000.0]
# S2 = [1.0 100.0; 
#      10.0 2000.0]


S1 = [10000.0 0.01; 
     1.0 20000.0]
S2 = [0.01 10000.0; 
     1.0 20000.0]



@show cond(S1)
@show cond(S2)

@assert norm(S1) ≈ norm(S2) # Vergleichbarkeit...

b = [2.0, 1.0]



# exact solution

x1 = S1 \ b 
norm(S1*x1 - b)

x2 = S2 \ b 
norm(S2*x2 - b)


# solution if S has errors (not b)

#ΔS1 = 0.05*S1 #  Achtung! hier gibt es das Problem noch nicht! 
#ΔS2 = 0.05*S2

N = 10000
Δx1_list = Vector{Float64}(undef,N)
Δx2_list = Vector{Float64}(undef,N)
max_rel_integral_error = 0.005

for i in 1:N
     


     ΔS1 = S1.*rand(2,2)*max_rel_integral_error.*rand_pm1_matrix(2, 2)
     ΔS2 = S2.*rand(2,2)*max_rel_integral_error.*rand_pm1_matrix(2, 2)
     Δb = b.*rand(2,1)*max_rel_integral_error.*rand_pm1_matrix(2, 1) 

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
