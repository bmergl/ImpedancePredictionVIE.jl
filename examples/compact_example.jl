using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.1)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 

##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

# qs3D = BEAST.SingleNumQStrat(6)
# qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
# qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


sol = IP.solve(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing), #IP.constantmaterial(1000.0, nothing),#, # #, #
    κ0 = 0.000001,
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D 
)

# save
# dataname = "test_solve1" # for JLD2 save
# jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, sol) 




##

##
using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes
using JLD2
using Plots
using Plotly
#load 
dataname = "test_solve1_fine"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
sol = load(datapath, "sol")
md = load(datapath, "md")

# dataname1 = "test_solve1"
# datapath1 = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname1.jld2"
# sol1 = load(datapath1, "sol")
# md1 = load(datapath1, "md")

norm(sol.S-sol1.S)

norm(sol1.S)
norm(sol.S)

norm(sol.u-sol1.u)/norm(sol1.u)

##

DIFF = log.(abs.(sol.S-sol1.S) .+eps(1.0))
p = Plots.heatmap(DIFF, title = "2D Rasterplot der Differenzmatrix")
yflip!(p; size=(1200,1000))



##

##


# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, sol.u_J, md.X)#, type=Float64)
J_ana = IP.solution_J_ana(md.body, sol.material, md, sol, points, J_MoM)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, sol.u_J, md.X)
J_ana2 = IP.solution_J_ana(md.body, sol.material, md, sol, points2, J_MoM2)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solution_J_ana(md.body, sol.material, md, sol, points3, J_MoM3)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
display("")
I_top, I_bottom = getcurrent(md, sol) # <-------------------------------getcurrent prüfen!!! für ALLE solve   !! Betrag, Richtung, n orientierung der dreiecke sign(volume)...
#@show I_top, I_bottom
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
display("")

# Potential: Randknoten vs. Analytisch
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)


##
display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(md.Γ_c)))


## facecurrents Tests

# J_n auf Γ_c mittels u_Jn
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(u_Φ, md.y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

# J_n auf Γ_c mittels u_J d.h. mittels ntrace der Volumenlösung
fcr0, geo0 = facecurrents(sol.u_J, md.ntrcX)
Plotly.plot(patch(geo0, fcr0))

for el in fcr0
    el != 0.0 && @show el
end