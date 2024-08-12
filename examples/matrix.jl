using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays
using BEAST
using FastBEAST
using ImpedancePredictionVIE
using CompScienceMeshes
using IterativeSolvers

using JLD2

using Plots
using Ploty

md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumSauterQstrat(3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)
# qs3D = BEAST.SingleNumQStrat(6)
# qs4D = BEAST.DoubleNumWiltonSauterQStrat(6,6,6,6,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
# qs5D6D = BEAST.SauterSchwab3DQStrat(6,6,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


material = IP.constantmaterial(10.0, nothing)
κ0 = 1.0
ϵ0 = nothing
ω = nothing
potential_top = 0.5
potential_bottom = -0.5

# Material
κ, ϵ = material()
τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)

# Excitation
v_top = ones(length(md.topnodes)) * potential_top
v_bottom = ones(length(md.bottomnodes)) * potential_bottom
v = vcat(v_top, v_bottom)

# Operators row 1
B11_Γ = BEAST.Identity() #Verschwindet!!! <---- klären....
B11_ΓΓ = Helmholtz3D.doublelayer(gamma = T(0.0), alpha = 1.0)# + (-1/2)*Identity()

B12_ΓΓ =  IP.MaterialSL(T(0.0), 1.0, inv_τ)

B13_ΓΓ = IP.MaterialSL(T(0.0), 1.0, χ)
B13_ΓΩ = IP.gradG_ΓΩ(T(0.0), -1.0, χ)


# Operators row 2
B21_ΓΓ = Helmholtz3D.hypersingular(gamma = T(0.0), beta = -1.0) # das versteckt alpha std. Null für gamma =0.0!!!  # -1.0 siehe BEAST & Steinbach 6.5, Ja, Float64 statt T

B22_Γ = IP.MatId(T(-1.0), inv_τ)
B22_ΓΓ = IP.MaterialADL(T(0.0), 1.0, inv_τ)

B23_ΓΓ = IP.MaterialADL(T(0.0), 1.0, χ)
B23_ΓΩ = IP.ncgrad_gradGc_ΓΩ(T(0.0), 1.0, χ)


# Operators row 3
B31_ΓΓ = Helmholtz3D.doublelayer(gamma = T(0.0), alpha = 1.0)# + (-1/2)*Identity()
B31_ΩΓ = IP.div_ngradG_ΩΓ(T(0.0), -1.0, x->1.0)

B32_ΓΓ = IP.MaterialSL(T(0.0), 1.0, inv_τ)
B32_ΩΓ = IP.div_G_ΩΓ(T(0.0), -1.0, inv_τ)

B33_Ω =  IP.MatId(T(-1.0), inv_τ) # different_tau scalartype(MaterialIdentity) = alpha ...
B33_ΓΓ = IP.MaterialSL(T(0.0), 1.0, χ)
B33_ΓΩ = IP.n_gradG_ΓΩ(T(0.0), -1.0, χ)
B33_ΩΓ = IP.div_G_ΩΓ(T(0.0), -1.0, χ)
B33_ΩΩ = IP.div_gradG_ΩΩ(T(0.0), 1.0, χ)

# fns spaces
y_d = md.y_d
y = md.y
w = md.w
X = md.X
ntrcX = md.ntrcX 
#wdual = duallagrangec0d1(md.Γ, md.Γ_c) # DoubleNumSauterQstrat ist ohne wilton!
#ydual = duallagrangecxd0(md.Γ, md.Γ_nc)
t_eq1 = w
t_eq2 = y

# assemble
B11 = assemble(B11_Γ, t_eq1, y) + assemble(B11_ΓΓ, t_eq1, y) -(1/2)*assemble(BEAST.Identity(), t_eq1, y)
B12 = assemble(B12_ΓΓ, t_eq1, w)
B13 = assemble(B13_ΓΓ, t_eq1, ntrcX) + assemble(B13_ΓΩ, t_eq1, X) 

B21 = assemble(B21_ΓΓ, t_eq2, y)
B22 = assemble(B22_Γ, t_eq2, w) + assemble(B22_ΓΓ, t_eq2, w)
B23 = assemble(B23_ΓΓ, t_eq2, ntrcX) + assemble(B23_ΓΩ, t_eq2, X)

B31 = assemble(B31_ΓΓ, ntrcX, y) -(1/2)*assemble(BEAST.Identity(), ntrcX, y) + assemble(B31_ΩΓ, X, y) 
B32 = assemble(B32_ΓΓ, ntrcX, w) + assemble(B32_ΩΓ, X, w)
B33 = assemble(B33_Ω, X, X) +
    assemble(B33_ΓΓ, ntrcX, ntrcX) +
    assemble(B33_ΓΩ, ntrcX, X) +
    assemble(B33_ΩΓ, X, ntrcX) + 
    assemble(B33_ΩΩ, X, X)

R11 = -assemble(B11_Γ, t_eq1, y_d) -assemble(B11_ΓΓ, t_eq1, y_d) +(1/2)*assemble(BEAST.Identity(), t_eq1, y_d)
R21 = -assemble(B21_ΓΓ, t_eq2, y_d)
R31 = -assemble(B31_ΓΓ, ntrcX, y_d) +(1/2)*assemble(BEAST.Identity(), ntrcX, y_d) -assemble(B31_ΩΓ, X, y_d)

#error("STOP HERE!")

#ROW1 = hcat(B21,B22,B23) #swap 1&2
#ROW2 = hcat(B11,B12,B13) # "-"
#ROW3 = hcat(B31,B32,B33)
#S = vcat(ROW1,ROW2,ROW3)
R = vcat(R21,R11,R31) # "-"
b = R*v


rownrmB21 = norm(B21) #[norm(B21[i, :]) for i in 1:size(B21, 1)] #
rownrmB12 = norm(B12) #[norm(B12[i, :]) for i in 1:size(B12, 1)] #
rownrmB33 = norm(B33) #[norm(B33[i, :]) for i in 1:size(B33, 1)] #

α_T = 1.0 #./norm(b[1:length(w.fns)]) 
α_B = 1.0 #./(α_T.*rownrmB21)

β_T = 1.0 #./norm(b[length(w.fns)+1:length(w.fns)+length(y.fns)]) 
β_B = 1.0 #./(β_T.*rownrmB12)

γ_T = 1.0 #./norm(b[length(w.fns)+length(y.fns)+1:end]) 
γ_B = 1.0 #./(γ_T.*rownrmB33)

@show α_T
@show α_B
@show β_T
@show β_B
@show γ_T
@show γ_B

ROW1 = α_T.*hcat(α_B.*B21, β_B.*B22, γ_B.*B23) #swap 1&2
ROW2 = β_T.*hcat(α_B.*B11, β_B.*B12, γ_B.*B13) # "-"
ROW3 = γ_T.*hcat(α_B.*B31, β_B.*B32, γ_B.*B33)
S = vcat(ROW1,ROW2,ROW3)
R = vcat(α_T.*R21, β_T.*R11, γ_T.*R31) # new
b = R*v

@show opnorm(B11)
@show opnorm(B12)
@show opnorm(B13)

@show opnorm(B21)
@show opnorm(B22)
@show opnorm(B23)

@show opnorm(B31)
@show opnorm(B32)
@show opnorm(B33)



@show cond(B11)
@show cond(B12)
@show cond(B13)

@show cond(B21)
@show cond(B22)
@show cond(B23)

@show cond(B31)
@show cond(B32)
@show cond(B33)


~,s,~ = svd(B11)

Plots.plot(s, marker=:x)

