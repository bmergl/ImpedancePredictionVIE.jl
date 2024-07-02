using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using Plots

using SauterSchwab3D



## # 1D


g = x1 -> 1.0
Ig_exact = 1.0 # ∫[0...1]g(x1)dx1


n_list = [i for i in 1:10000] #collect(5000:8000)
err_I_sum1 = []
err_I_sum2 = []
err_I_sum3 = []

for n in n_list
    qps = SauterSchwab3D._legendre(n,0,1)

    # V1
    I_sum1 = sum(w1*g(x1) for (x1, w1) in qps)
          
    # V2
    el_vec = [w1*g(x1) for (x1, w1) in qps]
    
    n == n_list[end] && @show minimum(el_vec)
    I_sum2 =  1.0
    for el in el_vec
        I_sum2 += el
    end
    #I_sum2 = I_sum2 -1.0

    # V3
    I_sum3 = sum(el_vec) # = sum([w1*g(x1) for (x1, w1) in qps])

    push!(err_I_sum1, abs(I_sum1-Ig_exact)/Ig_exact)
    push!(err_I_sum2, abs(I_sum2-Ig_exact)/Ig_exact)
    push!(err_I_sum3, abs(I_sum3-Ig_exact)/Ig_exact)

end

## 

#plt = plot(xscale=:log2)
plt = plot(n_list, err_I_sum1, label="sum(... in qps) - SauterSchwab3D", linewidth=4.0, size=(1000,800))
plot!(plt, n_list, err_I_sum2, label="add elements in sequence", linewidth=2.0)
plot!(plt, n_list, err_I_sum3, label="sum([... in qps])", linewidth=2.0)

##



sum(w1*w2*w3*w4*SauterSchwab3D.k5p_ce(INTEGRAND1, x1, x2, x3, y1, y2)
for (x1, w1) in qps for (x2, w2) in qps for (x3, w3) in qps for (y1, w4) in qps for (y2, w5) in qps)

##


# n 





qps = ce_ref.qps
A = [w1*w2*w3*w4*w5*SauterSchwab3D.k5p_ce(INTEGRAND1, x1, x2, x3, y1, y2)
		for (x1, w1) in qps for (x2, w2) in qps for (x3, w3) in qps for (y1, w4) in qps for (y2, w5) in qps]
sum(A)-exact1


@show s
s-exact1



