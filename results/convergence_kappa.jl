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


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumSauterQstrat(3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


mat = IP.pwlinx([[100.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-0.01/2, -0.01/6, 0.01/6, 0.01/2])

h_vec = [0.0018, 0.00135, 0.0009, 0.00045]#, 0.0003]
md_vec = []
solmom_vec = []
solfem_vec = []
err_J_mom_vec = []
err_J_fem_vec = []
err_I_mom_vec = []
err_I_fem_vec = []

for h in h_vec

    md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = h)

    solmom, S, R = IP.solve0(;
        md = md, 
        material = mat,
        κ0 = 1.0,
        ϵ0 = nothing,
        ω = nothing, 
        potential_top = 0.5, 
        potential_bottom = -0.5,
        qs3D = qs3D, 
        qs4D = qs4D, 
        qs5D6D = qs5D6D,
    )

    solfem, S, R = IP.solvefem(;
    md = md, 
    material = mat,
    κ0 = 1.0,
    ϵ0 = nothing,
    ω = nothing,
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    )

    #solmom = solfem#!!!!!!!! deaktivieren !!!!!!!!!

    push!(md_vec, md)
    push!(solmom_vec, solmom)
    push!(solfem_vec, solfem)


    # Stomdichte
    range_ = range(-0.0049,stop=0.0049,length=25)
    points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
    J_MoM = IP.grideval(points, solmom.u_J, md.X)
    J_ana = IP.solution_J_ana(md.body, solmom.material, md, solmom, points, J_MoM)
    J_FEM = IP.grideval(points, solfem.u_J, md.X)
    
    push!(err_J_mom_vec, norm(J_MoM-J_ana)/norm(J_ana))
    push!(err_J_fem_vec, norm(J_FEM-J_ana)/norm(J_ana))

    # Strom 
    I_ana = IP.solution_I_ana(md.body, solmom.material, md, solmom)
    I_top2_MoM, I_bottom2_MoM = IP.getcurrent2(md, solmom)
    I_top2_FEM, I_bottom2_FEM = IP.getcurrent2(md, solfem)

    push!(err_I_mom_vec, norm(I_top2_MoM-I_ana)/norm(I_ana))
    push!(err_I_fem_vec, norm(I_top2_FEM-I_ana)/norm(I_ana))
end

dataname = "kappaconv1" # for JLD2 save
jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md_vec, solmom_vec, solfem_vec, h_vec, err_J_mom_vec, err_J_fem_vec) 

##



## Konvergenz: h,errJ #########################################

plt = Plots.plot(h_vec, err_J_mom_vec, line=:scatter, marker = :utriangle, label = "Error J MoM", legend=:bottomright, xscale=:log10, yscale=:log10, minorgrid=true, color = 1)
plot!(plt,h_vec, err_J_mom_vec, label = "", color = 1)
plot!(plt,h_vec, err_J_fem_vec, label = "Error J FEM", color = 2, line=:scatter, marker = :dtriangle)
plot!(plt,h_vec, err_J_fem_vec, label = "", color = 2)
plot!(plt,size=(500,400))
#plot!(plt,h_vec,10*h_vec.^0.5)
#plot!(plt,h_vec,100*h_vec.^1.0)
#plot!(plt,h_vec,100000*h_vec.^2.0)
#plot!(plt,h_vec, err_J_fem_vec, line=:scatter, marker = :dtriangle, label = "Error J FEM", xscale=:log10, yscale=:log10)

xlims!(2e-4, 3e-3)
ylims!(1e-2, 1e+0)

#xlabel!("h in m")
#ylabel!("rel. Error")



## Konvergenz: h,errI (mittels SWG) ###########################

plt = Plots.plot(h_vec, err_I_mom_vec, line=:scatter, marker = :utriangle, label = "Error I MoM", legend=:bottomright, xscale=:log10, yscale=:log10, minorgrid=true, color = 1)
plot!(plt,h_vec, err_I_mom_vec, label = "", color = 1)
plot!(plt,h_vec, err_I_fem_vec, label = "Error I FEM", color = 2, line=:scatter, marker = :dtriangle)
plot!(plt,h_vec, err_I_fem_vec, label = "", color = 2)
plot!(plt,size=(500,400))

#xlims!(3e-4, 3e-3)
#ylims!(1e-2, 1e+0)

#xlabel!("h in m")
#ylabel!("rel. Error")



## Φ auf Γ ######################################################
nr = 1
sol = solmom_vec[nr]
#sol = solfem_vec[nr] 
md = md_vec[nr]

## MoM
u_Φ_full = vcat(sol.u_Φ,sol.v)
y_full = BEAST.LagrangeBasis{1,0,3}(md.y_d.geo, vcat(md.y.fns,md.y_d.fns), vcat(md.y.pos,md.y_d.pos))
fcr1, geo1 = facecurrents(u_Φ_full, y_full)
Plotly.plot(patch(geo1, fcr1))

## FEM
u_Φ_full = vcat(sol.u_Φ,sol.v)
Y_full = BEAST.LagrangeBasis{1,0,4}(md.Ω, vcat(md.Y.fns,md.Y_d.fns), vcat(md.Y.pos,md.Y_d.pos))
fcr1, geo1 = facecurrents(u_Φ_full, strace(Y_full,md.Γ))
Plotly.plot(patch(geo1, fcr1))



## J_n auf Γ_c mittels u_Jn ##########################################
nr = 1
sol = solmom_vec[nr]
#sol = solfem_vec[nr] 
md = md_vec[nr]

## MoM (only)
fcr0, geo0 = facecurrents(sol.u_Jn, md.w)
Plotly.plot(patch(geo0, fcr0))

## MoM/FEM mittels ntrace der Volumenlösung u_J  
fcr3, geo3 = facecurrents(sol.u_J, md.ntrcX)
Plotly.plot(patch(geo3, fcr3))



## J-Vektorfeld  #######################################################
nr = 1
sol = solmom_vec[nr]
#sol = solfem_vec[nr] 
md = md_vec[nr]

# Stomdichte MoM/FEM 
range_ = range(-0.0049,stop=0.0049,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = IP.grideval(points, sol.u_J, md.X)#, type=Float64)
J_ana = IP.solution_J_ana(md.body, sol.material, md, sol, points, J_MoM)
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Ω, p, nodesize = 0.02, nodecolor = "blue", linewidth = 1.0, linecolor = "green")
p = Visu.mesh(md.Γ_c, p, nodesize = 0.02, nodecolor = "blue", linewidth = 1.0, linecolor = "blue")
display(Visu.fieldplot(points, J_MoM, 0.0007, p)) #Visu.mesh(md.Γ_c)



## J_z on an x-line at y0, z0  ##################################################
nr = 1
sol2 = solmom_vec[nr]
sol1 = solfem_vec[nr] 
md = md_vec[nr]

y0 = 0.0
z0 = 0.0
x_range = range(-md.body.L_x/2, stop = md.body.L_x/2, length = 200)
points_x = [point(x, y0, z0) for x in x_range]
x = collect(x_range)

J_FEM_x = IP.grideval(points_x, sol1.u_J, md.X)
~, ~, J_z_FEM   = pointlist2xyzlist(J_FEM_x)
J_MoM_x = IP.grideval(points_x, sol2.u_J, md.X)
~, ~, J_z_MoM   = pointlist2xyzlist(J_MoM_x)
J_ana_x = IP.solution_J_ana(md.body, sol2.material, md, sol2, points_x, J_MoM_x)
~, ~, J_z_ana = pointlist2xyzlist(J_ana_x)

## Plot
Plots.plot(x, -J_z_ana, label = "-J_z_ana")#, size=(700,600))
plot!(x, -J_z_FEM, label = "-J_z_FEM")
plot!(x, -J_z_MoM, label = "-J_z_MoM")
plot!(size=(700,600))
#xlims!(0.0, 1.0)
#ylims!(1600, 2000)
title!("J_z(x, y0, z0)")
xlabel!("x")

##






