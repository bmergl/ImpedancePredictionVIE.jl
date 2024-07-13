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
# qs3D = BEAST.SingleNumQStrat(1)
# qs4D = BEAST.DoubleNumWiltonSauterQStrat(1,1,1,1,1,1,1,1) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
# qs5D6D = BEAST.SauterSchwab3DQStrat(1,1,1,1,1,1)

qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

# qs3D = BEAST.SingleNumQStrat(6)
# qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
# qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


sol, S, R = IP.solve(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.constantmaterial(0.1, 10000*8.8541878188*1e-12), #IP.pwlinx([[1.0, 2.0],[4.0, 10.0],[20.0, 5.0]], nothing, [-md.body.L_x/2, -1/6, 1/6, md.body.L_x/2]),   #IP.general_material(κ, nothing),  #  IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing), ,#, # #, #
    κ0 = 0.05, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
    ϵ0 = 1.0*8.8541878188*1e-12, # Größenordnung??? was ist sinnvoll?
    ω = 2*pi*1000.0, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D 
)

# save
# dataname = "test" # for JLD2 save
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
dataname = "layers3_stdquad_finest"
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
display("")
I_top, I_bottom = getcurrent(md, sol) # <-------------------------------getcurrent prüfen!!! für ALLE solve   !! Betrag, Richtung, n orientierung der dreiecke sign(volume)...
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
I_top2, I_bottom2 = IP.getcurrent2(md, sol)
display("")
@show norm(I_top2-I_ana)/norm(I_ana)
@show norm(I_bottom2-I_ana)/norm(I_ana)
display("")

# Potential: Randknoten vs. Analytisch
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)


##
display(Visu.fieldplot(points, J_MoM_abscomp, 0.0007, Visu.mesh(md.Γ_c)))

# Z mittels w auf Γ_c
I_top
I_bottom
Z = (sol.potential_top-sol.potential_bottom)/I_top

# Z mittels SWG auf Γ_c
I_top2
I_bottom2
Z = (sol.potential_top-sol.potential_bottom)/I_top2

# Nenner ≈ 1 ;) wow ...
1 + sol.ω^2*((1/0.1)*0.01/(0.01^2))^2  *  (10000*8.8541878188*1e-12*(0.01^2)/0.01)^2

# C 


# Analyse Complexe Stromdichte
J_MoM_abscomp = Vector{SVector{3, Float64}}()
for J_vec in J_MoM
    push!(J_MoM_abscomp, -abs.(J_vec))

    # norm(J_vec) > 1e-12 && @show J_vec
    # @show J_vec[1]
    # @show J_vec[2]
    # @show J_vec[3]
    # @show angle(J_vec[1])*360/(3*pi)
    # @show angle(J_vec[2])*360/(3*pi)
    # @show angle(J_vec[3])*360/(3*pi)
    # println()
end


## facecurrents Tests

# J_n auf Γ_c mittels u_Jn
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(sol.u_Φ, md.y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

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
#xlims!(0.0, 1.0)
#ylims!(-0.3, 0.0)
title!("J_z(x, y0, z0)")
xlabel!("x")



## It solver test: STOP SINNLOS OHNE PRECOND... #############
# n = 1000
# A= 100*LinearAlgebra.I(n) +  randn(n,n)
# cond(A)
# b = randn(n)
#b = vcat([100],zeros(99))
#x = IterativeSolvers.gmres(A, b; maxiter=1000)
#@show norm(A*x-b)

u_it, ch = IterativeSolvers.gmres(S, sol.b; maxiter=10000, log=true)
@show norm(sol.u-u_it)
@show ch


## H-Matrix ... ja..... #######################################

B33_ΩΩ = IPVIE.B33_ΩΩ(alpha = 1.0, gammatype = Float64, chi = x->1.0)

@time M = assemble(B33_ΩΩ, md.X, md.X)
#

Mh = @views hassemble(B33_ΩΩ, md.X, md.X;
    treeoptions=BoxTreeOptions(nmin=100),
    compressor = FastBEAST.ACAOptions(tol=1e-4),
    #quadstratcbk=BEAST.defaultquadstrat(operator, test_functions, trial_functions),
    #quadstratfbk=BEAST.defaultquadstrat(operator, test_functions, trial_functions),
    multithreading=true,
    verbose=true
)
println("done!")
@show norm(M-Matrix(Mh)) # Matrix(Mh) geht gar nicht...

##