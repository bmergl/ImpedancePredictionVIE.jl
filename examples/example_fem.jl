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
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
print("tehrahedrons: ", length(md.Ω.faces))





Y = lagrangec0d1(md.Ω, dirichlet = false) # sinnvoll anpassen...sollte einfach sein!
md.Ω.vertices
X = md.X

OP1 = IP.MatLoc(1.0, x->1.0)

#import BEAST

OP1 = IP._grad_Ω(1.0, x->1.0)
BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.NDLCDRefSpace{T}, ::BEAST.LagrangeRefSpace{T,D1,4}) where {T,D1} = BEAST.SingleNumQStrat(6)
function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.NDLCDRefSpace{T},
    f::BEAST.LagrangeRefSpace{T,Deg,4}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
    # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

assemble(OP1, X, Y)


OP2 = IP._div_Ω(1.0, x->1.0)
BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.LagrangeRefSpace{T,D1,4}, ::BEAST.NDLCDRefSpace{T}) where {T,D1} = BEAST.SingleNumQStrat(6)
function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.LagrangeRefSpace{T,Deg,4},
    f::BEAST.NDLCDRefSpace{T}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
    # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

assemble(OP2, Y, X)




# X2 = nedelecd3d(md.Ω)

# divX2 = divergence(X2, X2.geo, X2.fns) 

# strcdivX2 = strace(divX2, boundary(divX2.geo))

# els, ad, nr = assemblydata(divX2)




# els, ad, nr = assemblydata(md.X)

# faces(els[1])[1]

# cell_faces_center = cartesian.(CompScienceMeshes.center.(faces(els[1])))

# function BEAST.strace(x::BEAST.LagrangeRefSpace{T, 0, 4, 4}, cell, localid, face) where {T}

#     t = ones(T, 1, 4)
#     return t
#     #T = scalartype(x)
#     t = zeros(T, 4, 4)
    
#     for i in 1:4
#         t[i,1] = T(1.0)
#     end

#     #cell_faces = faces(cell)
#     #cell_faces_center = cartesian.(CompScienceMeshes.center.(cell_faces))
#     #face_center = cartesian(CompScienceMeshes.center(face))

#     return t
# end
# divergence2(space::BEAST.NDLCDBasis, geo, fns) = BEAST.LagrangeBasis{0,-1,4}(geo, fns, space.pos)
# function BEAST.strace(X::BEAST.LagrangeBasis{0,-1}, geo, fns::Vector)
#     # dimension of the boundary
#     #n = dimension(geo)
#     # degree of the space
#     #d = 1
#     # number of local shape functions
#     #NF = binomial(n+d,n)

#     trpos = deepcopy(positions(X))

#     NF = 4
#     @show NF
#     return BEAST.LagrangeBasis{0,-1,NF}(geo, fns, trpos)
# end

# ##

# divX = divergence2(md.X, md.X.geo, md.X.fns) # damit 4 fns je cell
# typeof(divX)
# refspace(divX)
# BEAST.numfunctions(refspace(divX))
# # Funktionen
# strcdivX = strace(divX, md.Γ) 
# strcdivX.fns

# for el in strcdivX.fns
#     @show el
# end


# # HAUPTPROBLEME: X ist ja schon weniger... und 


# els, ad, nr = assemblydata(strcdivX)



# X = md.X
# grady = gradient(md.y)
# gradyd = gradient(md.y_d)


# assemble(BEAST.Identity(), md.X, md.X)
# OP1 = IP.MatLoc(1.0, x->1.0)
# assemble(OP1, md.y, strcdivX)




# X2 = nedelecd3d(md.Ω)
# divX = BEAST.divergence2(X2, X2.geo, X2.fns) # damit 4 fns je cell
# typeof(divX)
# refspace(divX)
# BEAST.numfunctions(refspace(divX))
# strcdivX = strace(divX, boundary(divX.geo))





# ntrace(X2, boundary(X2.geo))


# X1 = lagrangecxd0(md.Ω)

# refspace(X1)
# BEAST.numfunctions(refspace(X1))


# strcdivX = strace(divX, boundary(divX.geo))

# assemble(OP1, md.y, divX)


# rtest = refspace(md.y)
# rtrial = refspace(md.X)

# OP1 = IP._grad_Ω(1.0, x->1.0)

# BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.LagrangeRefSpace{T,D1,3}, ::BEAST.NDLCDRefSpace{T}) where {T,D1} = SingleNumQStrat(6)
# function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.LagrangeRefSpace{T,Deg,3} where {T,Deg},
#     f::BEAST.NDLCDRefSpace, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat)

#     u, w = trgauss(qs.quad_rule)
#     qd = [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
#     A = BEAST._alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end

# OP2 = IP._div_Γ(1.0, x->1.0)
# assemble(OP2, md.y, md.X)


# BEAST.refspace(md.y)

# TS = divergence(md.X)
# strcTS = strace(TS, boundary(md.Ω))

# assemble(OP1, md.y_d, md.y_d)




##

qs3D = BEAST.SingleNumQStrat(3)
BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D

# STANDARD-TESTMATERIAL: IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2])

sol, S, R = IP.solve(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.pwlinx([[1.0, 2.0],[4.0, 10.0],[20.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]),   #IP.general_material(κ, nothing),  #  IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing), ,#, # #, #
    κ0 = 1.0, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
    #matalloc = :center,
)

# save
dataname = "test" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, sol) 



