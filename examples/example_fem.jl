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





grady = gradient(md.y)
gradyd = gradient(md.y_d)
divX = divergence(md.X, md.X.geo, md.X.fns) # weil auf volumen genormt...

assemble(BEAST.Identity(), md.X, md.X)

OP1 = IP.MatLoc(1.0, x->1.0)
strcdivX = strace(divX, boundary(divX.geo))
assemble(OP1, md.y, divX)


rtest = refspace(md.y)
rtrial = refspace(md.X)

OP1 = IP._grad_Ω(1.0, x->1.0)

BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.LagrangeRefSpace{T,D1,3}, ::BEAST.NDLCDRefSpace{T}) where {T,D1} = SingleNumQStrat(6)
function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.LagrangeRefSpace{T,Deg,3} where {T,Deg},
    f::BEAST.NDLCDRefSpace, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat)

    u, w = trgauss(qs.quad_rule)
    qd = [(w[i], SVector(u[1,i], u[2,i])) for i in 1:length(w)]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

OP2 = IP._div_Γ(1.0, x->1.0)
assemble(OP2, md.y, md.X)


BEAST.refspace(md.y)

TS = divergence(md.X)
strcTS = strace(TS, boundary(md.Ω))

assemble(OP1, md.y_d, md.y_d)




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



