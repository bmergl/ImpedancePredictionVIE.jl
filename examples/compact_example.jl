using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly




meshdata = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.18)
dataname = "test1" # for JLD2 save

κ = x -> 2.0
κ0 = 1.0


# Quadstrat
qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


solution = IP.solve(;
    meshdata = meshdata, 
    κ = κ, 
    κ0 = κ0, 
    ϵ = nothing, 
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D 
)


# ACHTUNG !!!! ANPASSEN VON L_x usw von cuboid auf den gmsh command!!! für cyl extra command...
# NEUE BERECHUNG DER STROMSTÄRKE norm.(J_vec-J_ana_vec)
# Postproc.jl wo die auswertung zusammengefasst ist, integrals.jl auflösen NE Doch einfach in ip.jl





##















jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; meshdata) 












## this could be in a seperate script
using JLD2
dataname = "test1"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
all_var=load(datapath)

meshdata = all_var["meshdata"]
meshdata = load(datapath, "meshdata")






