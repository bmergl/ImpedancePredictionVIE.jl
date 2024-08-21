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
#Visu.mesh(md.Ω) 


##


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

# STANDARD-TESTMATERIAL: IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2])

sol, S, R = IP.solve0(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.constantmaterial(10.0,nothing), #IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), #IP.constantmaterial(10.0, nothing), # IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), #IP.constant_zsplit(100.0, nothing, 0.0001, 10.0, nothing), #IP.pwlinx([[1.0, 20.0],[40.0, 100.0],[200.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), # # IP.pwlinx([[1.0, 2.0],[4.0, 10.0],[20.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]),   #IP.general_material(κ, nothing),  #  IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing),, #  ,#, # #, #
    κ0 = 1.0, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
    #matalloc = :avg#:center,
)

# save
#dataname = "test" # for JLD2 save
#jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, sol) 





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
dataname = "fine_h0.0004"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"
sol = load(datapath, "sol")
md = load(datapath, "md")
# dataname1 = "test_solve1"
# datapath1 = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname1.jld2"
# sol1 = load(datapath1, "sol")
# md1 = load(datapath1, "md")

#norm(sol.S-sol1.S)
# norm(sol1.S)
# norm(sol.S)
# norm(sol.u-sol1.u)/norm(sol1.u)

# DIFF = log.(abs.(sol.S-sol1.S) .+eps(1.0))
# p = Plots.heatmap(DIFF, title = "2D Rasterplot der Differenzmatrix")
# yflip!(p; size=(1200,1000))

##


gr()
MATR = log.(abs.(S) .+1.0e-10)
p = Plots.heatmap(MATR, title = "2D Rasterplot des Betrags der Systemmatrix", aspect_ratio=:equal)
yflip!(p; size=(1200,1000))


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

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.0049) for x in range_xy for y in range_xy]
J_MoM3 = IP.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solution_J_ana(md.body, sol.material, md, sol, points3, J_MoM3)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3) 
println()

# Strom durch Platten
println()
I_top, I_bottom = getcurrent(md, sol)
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
I_top2, I_bottom2 = IP.getcurrent2(md, sol)
println()
@show norm(I_top2-I_ana)/norm(I_ana)
@show norm(I_bottom2-I_ana)/norm(I_ana)
println()

# Potential: Randknoten vs. Analytisch
#@warn "check FEM settings"
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol, FEM = false)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)


##
display(Visu.fieldplot(points, J_MoM, 0.0007, Visu.mesh(md.Γ_c)))


## facecurrents Tests

# J_n auf Γ_c mittels u_Jn
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(sol.u_Φ, md.y)
Plotly.plot(patch(geo1, fcr1))

# Φ auf Γ 
u_Φ_full = vcat(sol.u_Φ,sol.v)
y_full = BEAST.LagrangeBasis{1,0,3}(md.y_d.geo, vcat(md.y.fns,md.y_d.fns), vcat(md.y.pos,md.y_d.pos))
fcr1, geo1 = facecurrents(u_Φ_full, y_full)
Plotly.plot(patch(geo1, fcr1))



# J_n auf Γ_c mittels u_J d.h. mittels ntrace der Volumenlösung
fcr3, geo3 = facecurrents(sol.u_J, md.ntrcX)
Plotly.plot(patch(geo3, fcr3))

## x-line at y0, z0 - J_z only inside the sphere mesh valid!
y0 = 0.0
z0 = 0.0
x_range = range(-md.body.L_x/2, stop = md.body.L_x/2, length = 200)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

J_MoM_x = BEAST.grideval(points_x, sol.u_J, md.X)
~, ~, J_z   = pointlist2xyzlist(J_MoM_x)
J_ana_x = IP.solution_J_ana(md.body, sol.material, md, sol, points_x, J_MoM_x)
~, ~, J_z_ana = pointlist2xyzlist(J_ana_x)

## Plot

Plots.plot(x, -J_z_ana, label = "J_z_ana", size=(700,600))
plot!(x, -J_z, label = "J_z")
#plot!(x, J_z*100, label = "J_z_modified")
#xlims!(0.0, 1.0)
#ylims!(1600, 2000)
title!("J_z(x, y0, z0)")
xlabel!("x")

