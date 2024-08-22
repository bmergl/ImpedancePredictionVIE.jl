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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0001)
print("tehrahedrons: ", length(md.Ω.faces))
# hlp = 1.0
# dataname = ""
# jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, hlp) 


# dataname = "md_h0.0001m"
# datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
# md = load(datapath, "md")

##



qs3D = BEAST.SingleNumQStrat(3)
BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D


# STANDARD-TESTMATERIAL: IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2])

sol, S, R = IP.solvefem(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), #IP.constantmaterial(1.0, nothing), #IP.pwlinx([[1.0, 2.0],[4.0, 10.0],[20.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]),   #IP.general_material(κ, nothing),  #  IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing), ,#, # #, #
    κ0 = 1.0,
    ϵ0 = nothing, #1.0*IP.ε0,
    ω = nothing, #2*pi*1000.0, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    #matalloc = :center,
)


##
#dataname = "test" # for JLD2 save
#jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md)#, sol) 


##




##


# Stomdichte
range_ = range(-0.0049,stop=0.0049,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = IP.grideval(points, sol.u_J, md.X)#, type=Float64)
J_ana = IP.solution_J_ana(md.body, sol.material, md, sol, points, J_MoM)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.005,stop=0.005,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = IP.grideval(points2, sol.u_J, md.X)
J_ana2 = IP.solution_J_ana(md.body, sol.material, md, sol, points2, J_MoM2)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.0049
points3 = [point(x,y,0.0049) for x in range_xy for y in range_xy]
J_MoM3 = IP.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solution_J_ana(md.body, sol.material, md, sol, points3, J_MoM3)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
println()
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
I_top2, I_bottom2 = IP.getcurrent2(md, sol)
println()
@show norm(I_top2-I_ana)/norm(I_ana)
@show norm(I_bottom2-I_ana)/norm(I_ana)
println()

# Potential: Randknoten vs. Analytisch ACHTUNG!!! DIE FUNKTION GEHT NICHT MEHR
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol; FEM = true)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)

##
#display(Visu.fieldplot(points, J_MoM, 0.0007, Visu.mesh(md.Γ_c)))


## facecurrents Tests

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
u_Φ_full = vcat(sol.u_Φ,sol.v)
Y_full = BEAST.LagrangeBasis{1,0,4}(md.Ω, vcat(md.Y.fns,md.Y_d.fns), vcat(md.Y.pos,md.Y_d.pos))
fcr1, geo1 = facecurrents(u_Φ_full, strace(Y_full,md.Γ))
Plotly.plot(patch(geo1, fcr1))

u_Φ_full = vcat(sol.u_Φ,sol.v)
Y_full = BEAST.LagrangeBasis{1,0,4}(md.Ω, vcat(md.Y.fns,md.Y_d.fns), vcat(md.Y.pos,md.Y_d.pos))
fcr1, geo1 = facecurrents(u_Φ_full, Y_full)
Plotly.plot(patch(geo1, fcr1))





# J_n auf Γ_c mittels u_J d.h. mittels ntrace der Volumenlösung
fcr3, geo3 = facecurrents(sol.u_J, md.ntrcX)
Plotly.plot(patch(geo3, fcr3))


## x-line at y0, z0 - J_z only inside the sphere mesh valid!
y0 = 0.0
z0 = 0.0
x_range = range(-md.body.L_x/2, stop = md.body.L_x/2, length = 400)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

J_MoM_x = IP.grideval(points_x, sol.u_J, md.X)
~, ~, J_z   = pointlist2xyzlist(J_MoM_x)
J_ana_x = IP.solution_J_ana(md.body, sol.material, md, sol, points_x, J_MoM_x)
~, ~, J_z_ana = pointlist2xyzlist(J_ana_x)

## Plot

Plots.plot(x, -J_z_ana, label = "J_z_ana", size=(700,600))
plot!(x, -J_z, label = "J_z")
#xlims!(0.0, 1.0)
#ylims!(-0.3, 0.0)
title!("J_z(x, y0, z0)")
xlabel!("x")




##
t

# nondirichletnodes = Vector{Int64}()
# for k in 1:length(md.Ω.vertices)
#     push_ = true
#     for n in md.dirichletnodes 
#         if k == n 
#             push_ = false
#             break
#         end
#     end
#     push_ && push!(nondirichletnodes, k)
# end
# @assert length(nondirichletnodes) + length(md.dirichletnodes) == length(md.Ω.vertices)

# # build continuous linear Lagrange elements on a 3D manifold
# function BEAST.lagrangec0d1(mesh, vertexlist::Vector, ::Type{Val{4}})

#     T = coordtype(mesh)
#     U = universedimension(mesh)

#     cellids, ncells = vertextocellmap(mesh)

#     Cells = cells(mesh)
#     Verts = vertices(mesh)

#     # create the local shapes
#     fns = Vector{BEAST.Shape{T}}[]
#     pos = Vector{vertextype(mesh)}()

#     sizehint!(fns, length(vertexlist))
#     sizehint!(pos, length(vertexlist))
#     for v in vertexlist

#         numshapes = ncells[v]
#         numshapes == 0 && continue

#         shapes = Vector{BEAST.Shape{T}}(undef,numshapes)
#         for s in 1: numshapes
#             c = cellids[v,s]
#             # cell = mesh.faces[c]
#             cell = Cells[c]

#             localid = something(findfirst(isequal(v), cell),0)
#             @assert localid != 0

#             shapes[s] = BEAST.Shape(c, localid, T(1.0))
#         end

#         push!(fns, shapes)
#         push!(pos, Verts[v])
#     end

#     NF = 4
#     BEAST.LagrangeBasis{1,0,NF}(mesh, fns, pos)
# end




# # Basis
# X = md.X
# Y_d = lagrangec0d1(md.Ω, md.dirichletnodes, Val{4})
# Y = lagrangec0d1(md.Ω, nondirichletnodes, Val{4})

# # Operators

# O = zeros(Float64,length(Y.fns),length(Y.fns))
# Õ = zeros(Float64,length(Y.fns),length(Y_d.fns)) 

# I = assemble(BEAST.Identity(), X, X)


# Op_A = IP._grad_Ω(1.0, x->2.0)
# A = assemble(Op_A, X, Y)
# Ã = -assemble(Op_A, X, Y_d)


# #els, ad, nr = assemblydata(Y)




# #Y_test = lagrangec0d1(md.Ω, dirichlet = false)
# #Visu.fnspos(Y)

# Op_B = IP._div_Ω(1.0, x->2.0)
# B = assemble(Op_B, Y, X)





















##








# Y = lagrangec0d1(md.Ω, dirichlet = false) # sinnvoll anpassen...sollte einfach sein!
# md.Ω.vertices
# X = md.X

# topnodes = realnodes(md.Γ_c_t)
# bottomnodes = realnodes(md.Γ_c_b)
# dirichletnodes = vcat(topnodes, bottomnodes)

# Y_d = lagrangec0d1(md.Ω, dirichletnodes, Val{3})
# Visu.fnspos(Y_d)



# nondirichletnodes = Vector{Int64}()
# for k in 1:length(md.Ω.vertices)
#     push_ = true
#     for n in dirichletnodes 
#         if k == n 
#             push_ = false
#             break
#         end
#     end
#     push_ && push!(nondirichletnodes, k)
# end
# @assert length(nondirichletnodes) + length(dirichletnodes) == length(md.Ω.vertices)

# Y = lagrangec0d1(md.Ω, nondirichletnodes, Val{3})
# Visu.fnspos(Y)

# BEAST.interior_and_junction_vertices(md.Ω, md.Γ)

# OP1 = IP.MatLoc(1.0, x->1.0)

# #import BEAST

# OP1 = IP._grad_Ω(1.0, x->1.0)
# BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.NDLCDRefSpace{T}, ::BEAST.LagrangeRefSpace{T,D1,4}) where {T,D1} = BEAST.SingleNumQStrat(6)
# function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.NDLCDRefSpace{T},
#     f::BEAST.LagrangeRefSpace{T,Deg,4}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
#     # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

#     o, x, y, z = CompScienceMeshes.euclidianbasis(3)
#     reftet = simplex(x,y,z,o)
#     qps = quadpoints(reftet, qs.quad_rule)
#     qd = [(w, parametric(p)) for (p,w) in qps]
#     A = BEAST._alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end

# assemble(OP1, X, Y)


# OP2 = IP._div_Ω(1.0, x->1.0)
# BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.LagrangeRefSpace{T,D1,4}, ::BEAST.NDLCDRefSpace{T}) where {T,D1} = BEAST.SingleNumQStrat(6)
# function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.LagrangeRefSpace{T,Deg,4},
#     f::BEAST.NDLCDRefSpace{T}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
#     # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

#     o, x, y, z = CompScienceMeshes.euclidianbasis(3)
#     reftet = simplex(x,y,z,o)
#     qps = quadpoints(reftet, qs.quad_rule)
#     qd = [(w, parametric(p)) for (p,w) in qps]
#     A = BEAST._alloc_workspace(qd, g, f, tels, bels)
#     return qd, A
# end

# assemble(OP2, Y, X)




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






