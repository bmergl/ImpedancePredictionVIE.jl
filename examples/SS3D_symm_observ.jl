using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays


function tet_circ_LHS(s) # simplex als input
    is_circ_lhs = false
    #s = simplex(mesh.vertices[face])
    t41=s.vertices[1]-s.vertices[4]
    t42=s.vertices[2]-s.vertices[4]
    t43=s.vertices[3]-s.vertices[4]

    @assert t41 == s.tangents[1]
    @assert t42 == s.tangents[2]
    @assert t43 == s.tangents[3]

    c = cross(t41,t42)
    dot(c,t43) > 0.0 && (is_circ_lhs = true)


    # t12=s.vertices[1]-s.vertices[2]
    # t13=s.vertices[1]-s.vertices[3]
    # t14=s.vertices[4]-s.vertices[1]
    # c2=cross(t12,t13)
    # @assert dot(c2,t14)<0


    """
    Ergebnis:

    1 2 3 CircRechteHand => n1
    => 4 ist in -n1  Richtung
    
    Alternativ direkt Linke hand
    """

    return is_circ_lhs
end

function findpoint(pointlist, point)
    for (j,p) in enumerate(pointlist)
        norm(p-point) < 1e-13 && return j
    end
    return error("point not in pointlist!")
end

#Mesh
v_list = [point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0)]

# indΩ1 = SVector(2,3,4,1)
# indΩ2 = SVector(3,2,5,1)

# STANDARD ########################
indΩ1 = SVector(1,2,4,3) #eigentlich erlaubt, macht aber massive Probleme bei evie bzgl Symm
indΩ2 = SVector(1,5,2,3)
###################################

# indΩ1 = SVector(1,4,3,2) #eigentlich erlaubt, macht aber massive Probleme bei evie bzgl Symm
# indΩ2 = SVector(1,3,5,2)


# indΩ1 = SVector(1,2,3,4) #RHR, ggf. verboten laut anderen meshes
# indΩ2 = SVector(1,3,2,5) #...

# indΩ1 = SVector(1,2,3,4) #RHR, ggf. verboten laut anderen meshes
# indΩ2 = SVector(1,2,5,3) #...

# indΩ1 = SVector(2,3,4,1) #LHR
# indΩ2 = SVector(2,1,5,3) #...

# indΩ1 = SVector(2,3,4,1) #LHR
# indΩ2 = SVector(2,1,5,3) #...   


s = simplex(v_list[indΩ1])
tet_circ_LHS(s)


f_list = [indΩ1, indΩ2]
mesh=Mesh(v_list, f_list)
bnd=boundary(mesh)

y = lagrangec0d1(bnd)
w = lagrangecxd0(bnd)
X = nedelecd3d(mesh)
ntrc = X -> ntrace(X,bnd)

Y = lagrangec0d1(mesh, dirichlet = false)
strc = Y -> strace(Y,bnd)

invtau = x -> 1.0
chi = x -> 1.0

Op1 = IPVIE2.B31_ΩΓ(gammatype=Float64, alpha= 1.0)
Op2 = IPVIE2.B32_ΩΓ(gammatype=Float64, alpha = 1.0, invtau = invtau)
# IPVIE2.B33_ΩΩ(gammatype=Float64, alpha = 1.0, chi=chi) hier auf diag alles bestens!
Op3 = VIE.hhboundary(gamma=0.0, alpha =1.0, tau=chi)

BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(6,6,6,6,6,6)




# Das hier entscheidend das es bei IPVIE Probleme macht PWC und SWG Kombi...
Z2 = assemble(Op2, X, w)    # SWG(3D) und PWC(2D)

#Symmetrie: (SWG und PWC derselben Fläche) * 4
p_a = SVector(0, 1/3, 1/3)
p_b = SVector(1/3, 0, 1/3)
p_c = SVector(0, 1/3, -1/3)
p_d = SVector(1/3, 0, -1/3)

i_Xa = findpoint(X.pos, p_a)
i_wa = findpoint(w.pos, p_a)

i_Xb = findpoint(X.pos, p_b)
i_wb = findpoint(w.pos, p_b)

i_Xc = findpoint(X.pos, p_c)
i_wc = findpoint(w.pos, p_c)

i_Xd = findpoint(X.pos, p_d)
i_wd = findpoint(w.pos, p_d)

@assert X.pos[i_Xa] == w.pos[i_wa]
@assert X.pos[i_Xb] == w.pos[i_wb]
@assert X.pos[i_Xc] == w.pos[i_wc]
@assert X.pos[i_Xd] == w.pos[i_wd]

Z2[i_Xa,i_wa]
Z2[i_Xb,i_wb]
Z2[i_Xc,i_wc]
Z2[i_Xd,i_wd]

@show abs(abs(Z2[i_Xb,i_wb])-abs(Z2[i_Xa,i_wa])) / abs(Z2[i_Xa,i_wa])   # ACHTUNG DAS IST NICHT DER RELATIVE
@show abs(abs(Z2[i_Xc,i_wc])-abs(Z2[i_Xa,i_wa])) / abs(Z2[i_Xa,i_wa])   # Fehler sondern der relative
@show abs(abs(Z2[i_Xd,i_wd])-abs(Z2[i_Xa,i_wa])) / abs(Z2[i_Xa,i_wa])   # UNTERSCHIED!!!

#@show abs(abs(Z2[i_Xb,i_wb])-abs(st)) / abs(st)

#Symmetrie: (SWG und PWC derselben Fläche) * 2 das sind die beiden besonderen FLächen
s_a = SVector(1/3, 1/3, 1/3)
s_b = SVector(1/3, 1/3, -1/3)

j_Xa = findpoint(X.pos, s_a)
j_wa = findpoint(w.pos, s_a)

j_Xb = findpoint(X.pos, s_b)
j_wb = findpoint(w.pos, s_b)

@assert X.pos[j_Xa] == w.pos[j_wa]
@assert X.pos[j_Xb] == w.pos[j_wb]

Z2[j_Xa,j_wa]
Z2[j_Xb,j_wb]

@show abs(abs(Z2[j_Xa,j_wa])-abs(Z2[j_Xb,j_wb])) / abs(Z2[j_Xb,j_wb])


""
##

#für nicht vektorielle Basis auch Probleme bei 5D Integralen?
Z = assemble(Op3, strc(Y), Y)

x_a = SVector(0,0,1.0)
x_b = SVector(0,0,-1.0)

i_sYa = findpoint(strc(Y).pos, x_a)
i_Ya = findpoint(strc(Y).pos, x_a)

i_sYb = findpoint(strc(Y).pos, x_b)
i_Yb = findpoint(Y.pos, x_b)

Z[i_sYa, i_Ya]
Z[i_sYb, i_Yb]

@show abs(abs(Z[i_sYb, i_Yb]) - abs(Z[i_sYa, i_Ya])) / abs(Z[i_sYa, i_Ya])
""

#
Z1 = assemble(Op1, X, y)    # SWG(3D) und LinLag(2D)      für uns gar nicht sooo wichtig
#SymmPaar1: X.pos[5] y.pos[4] symmetrisch zu X.pos[7] y.pos[4]
@show X.pos[5]
@show X.pos[7]
@show y.pos[4]
Z1[5,4]
Z1[7,4]