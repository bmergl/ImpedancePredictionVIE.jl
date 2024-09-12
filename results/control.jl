using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays
using BEAST
using ImpedancePredictionVIE
using CompScienceMeshes

using JLD2

using Plots
using Plotly

#md = IP.setup(; meshname = "coarsecube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0033)
md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
println("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

f_list = [1000000.0] #[1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0] # in Hz
Z_mom_list = []
Z_momhq_list = []
Z_fem_list = []
#Z_ana_list = []
mat = IP.constantmaterial(0.2, 10000.0*IP.ε0) #<------------------------------------------------------- KLÄREN!!!!
#mat = IP.constantmaterial(0.2, nothing)
solmom = nothing
solmom0 = nothing
Φtop = 0.5
Φbottom = -0.5

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


f = 1000000.0

solmom, S, R = IP.solve(;
    md = md, 
    material = mat,
    κ0 = 0.1,#0.1,
    ϵ0 = 1.0*IP.ε0,
    ω = 2*pi*f,  
    potential_top = Φtop, 
    potential_bottom = Φbottom,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
)

solmom0, S0, R0 = IP.solve0(;
    md = md, 
    material = mat,
    κ0 = 0.1,#0.1,
    ϵ0 = 1.0*IP.ε0,
    ω = 2*pi*f,  
    potential_top = Φtop, 
    potential_bottom = Φbottom,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
)

##

DIFF = log.(abs.(S-S0) .+eps(1.0))
p = Plots.heatmap(DIFF, title = "2D Rasterplot der Differenzmatrix")
yflip!(p; size=(1200,1000))

##
norm(solmom0.u-solmom.u)/norm(solmom.u)


##

DIFF = log.(abs.(R-R0) .+eps(1.0))
p = Plots.heatmap(DIFF, title = "2D Rasterplot der Differenzmatrix")
yflip!(p; size=(1200,1000))






