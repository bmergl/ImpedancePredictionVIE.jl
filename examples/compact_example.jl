using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


meshdata = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.09)
print("tehrahedrons:", length(meshdata.Ω.faces))

#Visu.mesh(meshdata.Γ) 



##



# für die materialtypen möglicherweise anonyme funktion des structs: constantmaterial("ϵ") liefert x -> ϵ
# aber wie konkret hier init - allg fall muss auch geplant sein

# analytische_Lösungen.jl erzeugen für die wichtigsten Fälle!

# Postproc.jl wo die auswertung zusammengefasst ist?



κ = x -> 2.0
κ0 = 1.0


# Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(4,4,4,4,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(4,4,6,6,6,6)

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
dataname = "solution_test1" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; solution) 


## # Constant Material Function ##########################################

# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, solution.u_J, meshdata.X)#, type=Float64)
display("Stomdichte Gesamtvolumen")
Jz_ana = -κ(point(0,0,0))*(solution.potential_top-solution.potential_bottom)
J_ana = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM))
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)


# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, solution.u_J, meshdata.X)
display("Stromdichte Mitte: Ebene z=0.0")
J_ana2 = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM2))
@show norm(J_MoM2-J_ana2)/norm(J_ana2)
#display(Visu.fieldplot(points2, J_MoM2, 1.0, Visu.mesh(Γ_c)))

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, solution.u_J, meshdata.X)
display("Stromdichte bei Platten: Ebene z=0.49")
J_ana3 = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM3))
@show norm(J_MoM3-J_ana3)/norm(J_ana3)
#display(Visu.fieldplot(points2, J_MoM2, 1.0, Visu.mesh(Γ_c)))

display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(meshdata.Γ_c)))









##




























## this could be in a seperate script
using JLD2
dataname = "test1"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
all_var=load(datapath)

meshdata = all_var["meshdata"]
meshdata = load(datapath, "meshdata")






