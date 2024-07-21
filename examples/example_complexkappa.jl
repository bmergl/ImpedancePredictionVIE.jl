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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0009)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)
#qs3D = BEAST.SingleNumQStrat(6)
#qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
#qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)
BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


sol, S, R = IP.solve(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.constant_xsplit(0.5, 3000.0*IP.ε0, 0.0008, 0.2, 10000.0*IP.ε0), # IP.constant_zsplit(0.5, 3000.0*IP.ε0, 0.0008, 0.2, 10000.0*IP.ε0), #  #, #IP.pwlinx([[1.0, 20.0],[40.0, 100.0],[200.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), # , # ##IP.constant_zsplit(0.2, 0.0, 0.0, 0.0, 10000.0*IP.ε0), #     #IP.general_material(κ, nothing),  #  I  ,#, # #, #
    κ0 = 0.1, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
    ϵ0 = 1.0*IP.ε0,
    ω = 2*pi*1000000.0, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
    #matalloc = :center
)

# save
dataname = "test" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md, sol) 




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
dataname = "solve1center"
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

##


# Stomdichte
range_ = range(-0.0049,stop=0.0049,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, sol.u_J, md.X)#, type=Float64)
J_ana = IP.solution_J_ana(md.body, sol.material, md, sol, points, J_MoM)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.005,stop=0.005,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, sol.u_J, md.X)
J_ana2 = IP.solution_J_ana(md.body, sol.material, md, sol, points2, J_MoM2)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.0049) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solution_J_ana(md.body, sol.material, md, sol, points3, J_MoM3)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
println()
I_top, I_bottom = getcurrent(md, sol)
#@show I_top
#@show I_bottom
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
I_top2, I_bottom2 = IP.getcurrent2(md, sol)
#@show I_top2
#@show I_bottom2
println()
@show norm(I_top2-I_ana)/norm(I_ana)
@show norm(I_bottom2-I_ana)/norm(I_ana)
println()

# Potential: Randknoten vs. Analytisch
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)

J_MoM_part = []
for J in J_MoM
    push!(J_MoM_part, -abs.(J))
end

##
display(Visu.fieldplot(points, J_MoM, 0.001, Visu.mesh(md.Γ_c)))
display(Visu.fieldplot(points, J_MoM_part, 0.001, Visu.mesh(md.Γ_c)))


## facecurrents Tests

# J_n auf Γ_c mittels u_Jn
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, abs.(fcr0)))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(abs.(sol.u_Φ), md.y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

# J_n auf Γ_c mittels u_J d.h. mittels ntrace der Volumenlösung

fcr3, geo3 = facecurrents(sol.u_J, md.ntrcX)
Plotly.plot(patch(geo3, fcr3))
Plotly.plot(patch(geo3, abs.(fcr3)))
# fcr3_corrected = [] # damit + und minus....
# for fcr_ in fcr3
#     fcr_new = abs.(fcr_)
#     abs(angle(fcr_) + pi)/pi < 0.01 && (fcr_new = -abs.(fcr_))
#     #@show angle(fcr_)
#     #@show abs(fcr_)
#     #println()
#     push!(fcr3_corrected, fcr_new)
# end
# Plotly.plot(patch(geo3, fcr3_corrected))



## x-line at y0, z0 - J_z only inside the sphere mesh valid!

y0 = 0.0
z0 = 0.0
x_range = range(-md.body.L_x/2, stop = md.body.L_x/2, length = 300)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

J_MoM_x = BEAST.grideval(points_x, sol.u_J, md.X)
~, ~, J_z   = pointlist2xyzlist(J_MoM_x)
J_ana_x = IP.solution_J_ana(md.body, sol.material, md, sol, points_x, J_MoM_x)
~, ~, J_z_ana = pointlist2xyzlist(J_ana_x)

Plots.plot(x, -abs.(J_z), label = "|J_z|", size=(700,600))
plot!(x, -abs.(J_z_ana), label = "|J_z_ana|")
#Plots.plot(x, -J_z, label = "J_z", size=(700,600))
#plot!(x, -J_z_ana, label = "J_z_ana")
#xlims!(0.0, 1.0)
ylims!(-40, 0)
#title!("|J_z(x, y0, z0)|")
xlabel!("x")

