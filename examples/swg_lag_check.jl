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
v_list = [point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(-1.0, 0.0, 0.0), point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0),     point(0.0, 1.0, -0.3), point(-0.5, 0.5, 0.5)]

indΓ1 = SVector(1,2,4)
indΓ2 = SVector(1,5,2)
indΓ3 = SVector(1,4,3)
indΓ4 = SVector(1,3,5)
bf_list = [indΓ1, indΓ2, indΓ3, indΓ4]#, indΩ5, indΩ6]
bmesh = Mesh(v_list, bf_list)


indΩ1 = SVector(1,5,2,6)
indΩ2 = SVector(1,6,2,7)
f_list = [indΩ1,indΩ2]# <-------------------------
mesh = Mesh(v_list, f_list)


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