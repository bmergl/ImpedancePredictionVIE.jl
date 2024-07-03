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
print("tehrahedrons: ", length(meshdata.Ω.faces))
#Visu.mesh(meshdata.Γ) 

##



# Quadstrat
qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,3,3,3,3) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D



solution = IP.solve(;
    meshdata = meshdata, 
    material = IP.constantmaterial(2.0, nothing), 
    κ0 = 1.0,
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D 
)
dataname = "test1" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; meshdata, solution) 






## 

# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, solution.u_J, meshdata.X)#, type=Float64)
J_ana = IP.solution_J_ana(meshdata, solution, points, J_MoM; body = meshdata.body, mat = solution.material)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, solution.u_J, meshdata.X)
J_ana2 = IP.solution_J_ana(meshdata, solution, points, J_MoM2; body = meshdata.body, mat = solution.material)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, solution.u_J, meshdata.X)
J_ana3 = IP.solution_J_ana(meshdata, solution, points, J_MoM3; body = meshdata.body, mat = solution.material)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
display("")
I_top, I_bottom = getcurrent(meshdata, solution)
@show I_top
@show I_bottom
I_ana = IP.solution_I_ana(meshdata, solution; body = meshdata.body, mat = solution.material)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
display("")

# Potential: Randknoten vs. Analytisch
u_Φ = solution.u_Φ
u_Φ_ana = IP.solution_Φ_ana(meshdata, solution; body = meshdata.body, mat = solution.material)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana) #nicht punktweise


##
display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(meshdata.Γ_c)))


## facecurrents Tests

# J_n auf Γ_c
fcr0, geo0 = facecurrents(solution.u_Jn, meshdata.w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(u_Φ, meshdata.y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

# Φ auf Γ_c    
fcr2, geo2 = facecurrents(ex, y_d)
Plotly.plot(patch(geo2, fcr2))






## this could be in a seperate script later ....
using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly
dataname = "solution_test1"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
all_var=load(datapath)

#meshdata = all_var["meshdata"]
solution = load(datapath, "solution")
meshdata = load(datapath, "meshdata")




