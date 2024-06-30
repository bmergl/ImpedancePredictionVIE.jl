using MKL

using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays
using Plotly


# hole Daten von main_example:

# vertices von Lag
[-0.5, -0.3779915320723491, -0.3779915320721015] #1
[-0.5, -0.2522937250027503, -0.3611620212027685] #2
[-0.5, -0.3333333333337957, -0.5]                #3
[-0.5, -0.1666666666675912, -0.5]                #4
[-0.5, -0.07888408893686816, -0.3316888586079301]#5
[-0.5, -0.2179091226067498, -0.2120159735025309] #6
[-0.5, -0.3630506580185957, -0.2528831044818892] #7

# swg
[-0.3503355194295261, -0.09839474790376397, -0.5]               # 8
[-0.5, -0.1666666666675912, -0.5]                               # doppelt 4
[-0.5, -0.07888408893686816, -0.3316888586079301]               # doppelt 5
[-0.3011349771082896, -0.3004935351285941, -0.2989213866644215] # 9
[-0.5, -0.2522937250027503, -0.3611620212027685]                # doppelt 2

# 1tet
[-0.5, -0.07888408893686816, -0.3316888586079301]
[-0.5, -0.2522937250027503, -0.3611620212027685]
[-0.3503355194295261, -0.09839474790376397, -0.5]
[-0.5, -0.1666666666675912, -0.5]

# 2tet
[-0.5, -0.2522937250027503, -0.3611620212027685]
[-0.3503355194295261, -0.09839474790376397, -0.5]
[-0.3011349771082896, -0.3004935351285941, -0.2989213866644215]
[-0.5, -0.07888408893686816, -0.3316888586079301]






# Mesh, Basis
v_list = [point(-0.5, -0.3779915320723491, -0.3779915320721015),point(-0.5, -0.2522937250027503, -0.3611620212027685),point(-0.5, -0.3333333333337957, -0.5),point(-0.5, -0.1666666666675912, -0.5),point(-0.5, -0.07888408893686816, -0.3316888586079301),point(-0.5, -0.2179091226067498, -0.2120159735025309),point(-0.5, -0.3630506580185957, -0.2528831044818892),       point(-0.3503355194295261, -0.09839474790376397, -0.5), point(-0.3011349771082896, -0.3004935351285941, -0.2989213866644215)]#Ω.vertices

# manip: gut um Near Singularities zu zeigen!
#v_list = [point(-0.5, -0.3779915320723491, -0.3779915320721015),point(-0.5, -0.2522937250027503, -0.3611620212027685),point(-0.5, -0.3333333333337957, -0.5),point(-0.5, -0.1666666666675912, -0.5),point(-0.5, -0.07888408893686816, -0.3316888586079301),point(-0.5, -0.2179091226067498, -0.2120159735025309),point(-0.5, -0.3630506580185957, -0.2528831044818892),       point(-0.35, -0.1, -0.4), point(-0.4511349771082896, -0.3504935351285941, -0.213866644215)]

# damit soll CommonFace5D näher untersucht werden!
v_list = [point(-0.5, -0.3779915320723491, -0.3779915320721015),point(-0.5, -0.2522937250027503, -0.3611620212027685),point(-0.5, -0.3333333333337957, -0.5),point(-0.5, -0.1666666666675912, -0.5),point(-0.5, -0.07888408893686816, -0.3316888586079301),point(-0.5, -0.2179091226067498, -0.2120159735025309),point(-0.5, -0.3630506580185957, -0.2528831044818892),       point(-0.3503355194295261, -0.09839474790376397, -0.3), point(-0.45011183980984204, -0.14319085394779416, -0.325)]#Ω.vertices



indΓ1 = SVector(2,3,1) #Γ_nc.faces[c_1]
indΓ2 = SVector(2,1,7) #Γ_nc.faces[c_2]
indΓ3 = SVector(2,4,3) #Γ_nc.faces[c_3]
indΓ4 = SVector(2,5,4) #Γ_nc.faces[c_4]
indΓ5 = SVector(2,6,5) #Γ_nc.faces[c_5]
indΓ6 = SVector(2,7,6) #Γ_nc.faces[c_6]
simplex(v_list[indΓ6]).normals[1]


bf_list = [indΓ1, indΓ2, indΓ3, indΓ4, indΓ5, indΓ6]
bmesh = Mesh(v_list, bf_list) 

indΩ1 = SVector(5,2,8,4) #Ω.faces[c1]
indΩ2 = SVector(2,8,9,5) #Ω.faces[c2]
f_list = [indΩ1,indΩ2]
mesh_ = Mesh(v_list, f_list)


# ntrc auf betreachtetem Abschnitt = 0 => BNDTERM = 0 => Volterm1=Volterm2 (ideal)
X_ = nedelecd3d(mesh_)
y_ = lagrangec0d1(bmesh, dirichlet = true)

i_swg = -1
for (i,el) in enumerate(X_.fns)
    length(el) == 2 && (i_swg = i)
end
pnts = Vector{SVector{3,Float64}}()
push!(pnts,X_.pos[i_swg])
X_.pos = pnts
X_.fns = [X_.fns[i_swg],]
bnd_ = boundary(mesh_)
ntrc_ = fnsspace -> ntrace(fnsspace, bnd_)


# manipuliere swg - sodass nur die halbe... STOPP  VERBOTEN siehe BND1
#X_.fns = [[X_.fns[1][1]],]
# manipuliere lagrange:
k = 4 # <------------------------ 4 ist problematisch!
y_.fns = [[y_.fns[1][k]],] 


display(X_.fns)
display(y_.fns)
@show numfunctions(X_)
@show numfunctions(y_)

#Visu.fnspos(X,Visu.mesh(mesh))
#Visu.mesh(mesh_,Visu.mesh(bmesh))
#Visu.mesh(mesh_)

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
Z_B = assemble(B23_ΓΓ, y_, ntrc_(X_))
Z_V = assemble(B23_ΓΩ, y_, X_)

Z_A = assemble(B23_alternativ, y_, X_)

@assert Z_B[1,1] == 0.0

@show (Z_V[1,1] - Z_A[1,1])/Z_A[1,1]


##

Z_B[1,1]
Z_V[1,1]

Z_A[1,1]

## 

j=1
c1 = X_.fns[j][1].cellid
c2 = X_.fns[j][2].cellid
vertc1 = mesh_.vertices[mesh_.faces[c1]]
vertc2 = mesh_.vertices[mesh_.faces[c2]]
plt = Visu.mesh(bmesh)
Visu.simplex(plt,vertc1)
Visu.simplex(plt,vertc2)

@assert length(y_.fns[1]) == 1 #sonst ...
cy = y_.fns[1][1].cellid
vertcy = bmesh.vertices[bmesh.faces[cy]]
lagcenter = cartesian(CompScienceMeshes.center(simplex(vertcy)))


Visu.points([lagcenter,X_.pos[j]],plt)




#Visu.points([y_.pos[1],X_.pos[j]],plt)


# Ja near Singularities durch CommonEdge5D oder CommonVertex5D konvergieren langsamer (je nach Winkel...), aber sie konvergieren
# "Grundverschiedene" Elemente nur im CommonFace5D Fall beobachtet - Frage ist der Alternativterm auch falsch?
# ...würde erklären warum es noch keine Konvergenz gibt... CommonFace5D unter GENERALVERDACHT!
# Kann es sein dass doch irgendeine identität Fehlt oder gibt es klare Gegenbeweise???