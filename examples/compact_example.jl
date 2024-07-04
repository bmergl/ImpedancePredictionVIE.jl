using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.09)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Γ) 

##


# Quadstrat
qs3D = BEAST.SingleNumQStrat(6)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


sol = IP.solve(;
    md = md, 
    material = IP.constant_zsplit(1000, nothing, 0.07, 100, nothing), #IP.constantmaterial(100.0, nothing), 
    κ0 = 550.0,
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D 
)

# save
dataname = "qs56_lowcontr_split_fine" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, sol) 


##



using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly

# load 
dataname = "qs56_lowcontr_split_coarse"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
sol = load(datapath, "sol")
md = load(datapath, "md")
@assert sol.material == IP.constant_zsplit(1000, nothing, 0.07, 100, nothing)



##


# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, sol.u_J, md.X)#, type=Float64)
J_ana = IP.solultion_J_ana(md, sol, points, J_MoM; body = md.body, mat = sol.material)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, sol.u_J, md.X)
J_ana2 = IP.solultion_J_ana(md, sol, points, J_MoM2; body = md.body, mat = sol.material)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solultion_J_ana(md, sol, points, J_MoM3; body = md.body, mat = sol.material)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
display("")
I_top, I_bottom = getcurrent(md, sol)
#@show I_top, I_bottom
I_ana = IP.solultion_I_ana(md, sol; body = md.body, mat = sol.material)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
display("")

# Potential: Randknoten vs. Analytisch
u_Φ = sol.u_Φ
u_Φ_ana = IP.solultion_Φ_ana(md, sol; body = md.body, mat = sol.material)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)


##
display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(md.Γ_c)))


## facecurrents Tests

# J_n auf Γ_c
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(u_Φ, md.y)
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
dataname = "sol_test1"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
all_var=load(datapath)

#md = all_var["md"]
sol = load(datapath, "sol")
md = load(datapath, "md")




