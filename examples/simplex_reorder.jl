using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays
using SauterSchwab3D




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

##

d = 1.0
tet_circ_LHS_var = false
tet = nothing
tri = nothing
p1 = point(rand(),rand(),rand())
p2 = point(rand(),rand(),rand())
p3 = point(rand(),rand(),rand())
p4 = point(rand(),rand(),rand())

p5 = point(rand(),rand(),rand())

# tet = simplex(p1,p2,p3,p4) # darf nicht geändert werden sonst passt das mit p4 nicht mehr!!!
# tri = simplex(p3,p4,p2) # beliebig...wird schon Fehler anzeigen wenn Punkte nicht passen
tet = simplex(p1,p2,p3,p4) # darf nicht geändert werden sonst passt das mit p4 nicht mehr!!!
tri = simplex(p5,p1,p2) # beliebig...wird schon Fehler anzeigen wenn Punkte nicht passen



tet_circ_LHS_var = tet_circ_LHS(tet) 
c_tri = cartesian(CompScienceMeshes.center(tri))
d = dot(tri.normals[1],p1+p2+p3+p4-c_tri)

@assert d < 0.0
@assert tet_circ_LHS_var == true


##


tau = x -> 1.0
OP = IPVIE2.B32_ΩΓ(alpha = 1.0, gammatype = Float64, invtau = tau)
v_list = [point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0), point(0.0, 0.0, 1.0), point(0.0, 0.0, -1.0)]
indΩ1 = SVector(1,2,4,3) #eigentlich erlaubt, macht aber massive Probleme bei evie bzgl Symm
indΩ2 = SVector(1,5,2,3)
f_list = [indΩ1, indΩ2]
mesh = Mesh(v_list, f_list)
w = lagrangecxd0(mesh)
qs = BEAST.quaddata(OP, refspace(w), refspace(w),[tet,], [tri,], BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3))
qr = BEAST.qr_boundary(OP, refspace(w), refspace(w), nothing, tet, nothing,  tri, qs, BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3))


#Common Face
function reorder_TEST(sing::SauterSchwab3D.Singularity5DPoint)

    # Find the permutation P of t and s that make
    # Pt = [P, A1, A2, A3]
    # Ps = [P, B1, B2]
  
    i = sing.T[1]-1
    I = circshift([1,2,3,4],-i)
    for k in 1:i
        reverse!(I,2,4)
    end 
    
    J = circshift([1,2,3],-sing.S[1]+1)

    return SVector{4}(I), SVector{3}(J)
end

# I1 = circshift([1,2,3,4],2)
# plist = tet.vertices[I1]
# tet2 = simplex(plist[1],plist[2],plist[3],plist[4])
# tet_circ_LHS(tet2) 
# tet_circ_LHS(tet)


qr.sing
I,J = reorder_TEST(qr.sing)

new_tet = simplex(tet.vertices[I])
new_tri = simplex(tet.vertices[J])
@assert tet_circ_LHS(new_tet) == true #nötig?

##

# Das ist die Darstellung die ANGEBLICH durch die reorder Funktion erreicht wird!
@show norm(new_tet[1] - new_tri[1]) #< 1.0e-14
@show norm(new_tet[2] - new_tri[2]) #< 1.0e-14
@show norm(new_tet[4] - new_tri[3]) #< 1.0e-14

# Das ist die Darstellung die laut example_cf_2.5d.jl nötig ist
# const P = simplex(pI,pII,pIV,pIII)    # immerhin ist der spezielle Punkt auch an 3. Stelle! 
# const Q = simplex(pI,pIII,pII)
# Pt = [P1, P2, A1, P3]
# Ps = [P1, P2, P3]

@show norm(new_tet[1] - new_tri[1]) #< 1.0e-14
@show norm(new_tet[2] - new_tri[3]) #< 1.0e-14
@show norm(new_tet[4] - new_tri[2]) #< 1.0e-14


##

allp = collect(BEAST.permutations([1,2,3,4]))
for el in allp
    @show el
    vert = tet.vertices[el]
    s1 = simplex(vert[1],vert[2],vert[3],vert[4])
    @show tet_circ_LHS(s1)
    display("--------------------")
end
# ! es gibt auch andere die den LHS tet erhalten würden, aber ist das überhaupt nötig wenn man eh zurücktauscht???

#  K = circshift([4,1,3,2],2) #circshift([1,2,3,4],1)
#  P_new2 = simplex(P[K[1]],P[K[2]],P[K[3]],P[K[4]])
#  @test is_CSM_tet(P_new2) == true