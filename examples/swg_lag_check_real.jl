using MKL

using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using Plots
using SphericalScattering

using ImpedancePredictionVIE



# Mesh, Basis
v_list = [point(-0.5, -0.5, -0.5), point(-0.5, -0.5, -0.3333333333328717), point(-0.3333333333328717, -0.5, -0.5), point(-0.3779915320715162, -0.5, -0.3779915320715162), point(-0.2528831044819755, -0.5, -0.3630506580181511), point(-0.3611620212023955, -0.5, -0.2522937250014939),     point(-0.5, -0.3779915320723491, -0.3779915320721015), point(-0.3779915320721015, -0.3779915320723491, -0.5)]

indΓ1 = SVector(4,5,6)
indΓ2 = SVector(2,4,6)
indΓ3 = SVector(1,4,2)
indΓ4 = SVector(1,3,4)
indΓ5 = SVector(3,5,4)
bf_list = [indΓ1, indΓ2, indΓ3, indΓ4, indΓ5]
bmesh = Mesh(v_list, bf_list)
Visu.mesh(bmesh)
Visu.points(v_list)

indΩ1 = SVector(7,1,4,8)
indΩ2 = SVector(4,7,2,1)
f_list = [indΩ1,indΩ2]# <-------------------------
mesh = Mesh(v_list, f_list)


# grid points (lag)
[-0.5, -0.5, -0.5] # 1
[-0.5, -0.5, -0.3333333333328717] # 2
[-0.3333333333328717, -0.5, -0.5] # 3
[-0.3779915320715162, -0.5, -0.3779915320715162] # 4
[-0.2528831044819755, -0.5, -0.3630506580181511] # 5
[-0.3611620212023955, -0.5, -0.2522937250014939] # 6



# swg points
[-0.5, -0.3779915320723491, -0.3779915320721015] # ist neu: 7
[-0.5, -0.5, -0.5]                               # ist 1
[-0.3779915320715162, -0.5, -0.3779915320715162] # ist 4
[-0.3779915320721015, -0.3779915320723491, -0.5] # ist neu: 8

[-0.3779915320715162, -0.5, -0.3779915320715162] # ist 4
[-0.5, -0.3779915320723491, -0.3779915320721015] # ist 7
[-0.5, -0.5, -0.3333333333328717]                # ist 2
[-0.5, -0.5, -0.5]                               # ist 1







# ntrc auf betreachtetem Abschnitt = 0 => BNDTERM = 0 => Volterm1=Volterm2 (ideal)
X = nedelecd3d(mesh)
y = lagrangec0d1(bmesh, dirichlet = true)

i_swg = -1
for (i,el) in enumerate(X.fns)
    length(el) == 2 && (i_swg = i)
end
pnts = Vector{SVector{3,Float64}}()
push!(pnts,X.pos[i_swg])
X.pos = pnts
X.fns = [X.fns[i_swg],]
bnd = boundary(mesh)
ntrc = fnsspace -> ntrace(fnsspace, bnd)


@show numfunctions(X)
@show numfunctions(y)

#Visu.fnspos(X,Visu.mesh(mesh))
#Visu.mesh(mesh,Visu.mesh(bmesh))


##

q = 6
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(q,q,q,6,6,6)
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,5,5,5,5)

# VIE Operators
χ = x -> 0.5 # ACHTUNG MUSS KONSTANT SEIN OHNE SPRÜNGE
B23_ΓΓ = IPVIE2.B23_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
B23_ΓΩ = IPVIE2.B23_ΓΩ(alpha = 1.0, gammatype = Float64, chi = χ)
B23_alternativ = IPVIE2.B23_alternativ(alpha = 1.0, gammatype = Float64, chi = χ)


# Assembly
Z_B = assemble(B23_ΓΓ, y, ntrc(X))
Z_V = assemble(B23_ΓΩ, y, X)

Z_A = assemble(B23_alternativ, y, X)

@assert Z_B[1,1] == 0.0

@show (Z_V[1,1] - Z_A[1,1])/Z_A[1,1]


##

Z_B[1,1]
Z_V[1,1]

Z_A[1,1]


##

j=1
c1 = X.fns[j][1].cellid
c2 = X.fns[j][2].cellid
vertc1 = mesh.vertices[mesh.faces[c1]]
vertc2 = mesh.vertices[mesh.faces[c2]]

plt = Visu.mesh(bmesh)
Visu.simplex(plt,vertc1)
Visu.simplex(plt,vertc2)
Visu.points([y.pos[1],X.pos[j]],plt)