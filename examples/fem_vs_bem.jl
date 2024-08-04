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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0008)
print("tehrahedrons: ", length(md.Ω.faces))


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

mat = IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2])

sol1, S1, R1 = IP.solvefem(;
    md = md, 
    material = mat,
    κ0 = 1.0,
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D,
)

sol2, S2, R2 = IP.solve0(;
    md = md, 
    material = mat, 
    κ0 = 1.0,
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
    #matalloc = :avg#:center,
)


## x-line at y0, z0 - J_z only inside the sphere mesh valid!
y0 = 0.0
z0 = 0.0
x_range = range(-md.body.L_x/2, stop = md.body.L_x/2, length = 200)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

J_FEM_x = IP.grideval(points_x, sol1.u_J, md.X)
~, ~, J_z_FEM   = pointlist2xyzlist(J_FEM_x)
J_MoM_x = IP.grideval(points_x, sol2.u_J, md.X)
~, ~, J_z_MoM   = pointlist2xyzlist(J_MoM_x)

J_ana_x = IP.solution_J_ana(md.body, sol.material, md, sol, points_x, J_MoM_x)
~, ~, J_z_ana = pointlist2xyzlist(J_ana_x)

## Plot

Plots.plot(x, -J_z_ana, label = "J_z_ana", size=(700,600))
plot!(x, -J_z_FEM, label = "J_z_FEM")
plot!(x, -J_z_MoM, label = "J_z_MoM")
#plot!(x, J_z*100, label = "J_z_modified")
#xlims!(0.0, 1.0)
#ylims!(1600, 2000)
title!("J_z(x, y0, z0)")
xlabel!("x")





