using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.18)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Γ) 

##

elements, ad, cells = assemblydata(md.X)

ad[4]

cells[end]
ad[1][1]

for (n,b) in ad[1000][1]
    @show n, b
end
xy
md.X.fns
# erstelle X_mat basis für trial in den speziellen Fällen!



mutable struct testop
    mf::Vector{Float64}
    v::Float64
end

a = Vector([0.1, 0.9, 1.3])
v_ini = 0.0
op1 = testop(a,v_ini)

function change!(op)
    #new_op = testop(op.mf,50.8)
    #op = new_op
    op.v = 50.8
    return nothing
end

change!(op1)
op1







# Quadstrat
qs3D = BEAST.SingleNumQStrat(6)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D