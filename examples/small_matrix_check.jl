using MKL

using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using Plots
using SphericalScattering

using ImpedancePredictionVIE


BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(6,6,6,6,6,6)
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,5,5,5,5)

# Mesh, Basis
v_list = [point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0), point(0.0,-1.0,0.0), point(-1.0,0.0,0.0), point(1.0,1.0,1.0), point(1.0,1.0,0.0)]
indΩ1 = SVector(1,2,4,3)
indΩ2 = SVector(1,4,2,6)
#indΩ3 = SVector(1,6,7,4)
#indΩ4 = SVector(1,7,3,4)

#indΩ5 = SVector(1,2,3,5)
#indΩ6 = SVector(1,2,5,6)
#indΩ7 = SVector(1,6,5,7)
#indΩ8 = SVector(1,5,3,7)

#indΩ9 = SVector(2,4,3,8)
#indΩ10 = SVector(2,3,9,8)

#f_list = [indΩ1, indΩ2, indΩ3, indΩ4, indΩ5, indΩ6, indΩ7, indΩ8, indΩ9, indΩ10] # <-------------------------
f_list = [indΩ1, indΩ2]#, indΩ5, indΩ6]
mesh = Mesh(v_list, f_list)


X = nedelecd3d(mesh)

bnd = boundary(mesh)
ntrc = fnsspace -> ntrace(fnsspace, bnd) 

y = lagrangec0d1(bnd; dirichlet = false)


@show numfunctions(X)
@show numfunctions(y)


#Visu.mesh(mesh)

##



# VIE Operators
χ = x -> 0.5 # ACHTUNG MUSS KONSTANT SEIN OHNE SPRÜNGE
B23_ΓΓ = IPVIE2.B23_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
B23_ΓΩ = IPVIE2.B23_ΓΩ(alpha = 1.0, gammatype = Float64, chi = χ)
B23_alternativ = IPVIE2.B23_alternativ(alpha = 1.0, gammatype = Float64, chi = χ)


# Assembly
Z_B = assemble(B23_ΓΓ, y, ntrc(X))
Z_V = assemble(B23_ΓΩ, y, X)
Z_0 = Z_B + Z_V

Z_A = assemble(B23_alternativ, y, X)

relErrM = abs.(Z_0 - Z_A)./abs.(Z_A)
@show maximum(relErrM)
display(relErrM)

## # Stop here

Z_B
Z_V
Z_0

Z_A


j_swg = 1

relErrM[:,j_swg]

Z_A[:,j_swg]
Z_0[:,j_swg]




X.pos[j_swg] # hat bei bnd keinen Beitrag zur matrix...
Z_B[:,j_swg]
Z_V[:,j_swg]

Z_A[:,j_swg]

##
Visu.points([y.pos[3],X.pos[j_swg]],Visu.mesh(mesh))

##



ind2Mind = CartesianIndices(relErrM)
pair_list = []
for (index, element) in enumerate(relErrM) 
    i = ind2Mind[index][1]
    j = ind2Mind[index][2]
    z0 = Z_0[i,j]
    za = Z_A[i,j]

    lim = 1.0
    if element > lim
        println("relErr: $element, Z_0[$i,$j]=$z0, Z_A[$i,$j]=$za")
        #abs(z0) > 0.00000001 && println("relErr: $element, Z_0[$i,$j]=$z0, Z_A[$i,$j]=$za")
        push!(pair_list,[i,j])
    end

end

display(pair_list)

##




plt = Visu.iplot()

Visu.mesh(mesh,plt)

pair = pair_list[2]

i = pair[1]
j = pair[2]  
Visu.points([y.pos[i],X.pos[j]],plt)





## #############################################







#mittlerer Radius

@show sum(norm.(X.pos))/length(X.pos)
r_list
@show sum(r_list)/length(r_list)







##

cond(Z_version1)

# abs1 = abs.(Z_Y)
# abs.(Z_0 - Z_Y)./abs1

# ~,sv,~ = svd(Z_0)
# plot(sv)

# ~,sv2,~ = svd(Z_Y)
# plot!(sv2)


# System matrix
Z_version1 = Z_I + Z_V + Z_B
Z_version2 = Z_I + Z_Y

cond(Z_version1)
cond(Z_version2)

# Solution of the linear system of equations
u_version1 = Z_version1 \ Vector(b)
u_version2 = Z_version2 \ Vector(b)


