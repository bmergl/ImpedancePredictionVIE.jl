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
#mat = IP.pwlinx([[10.0, 20.0],[40.0, 100.0],[200.0, 5.0]], nothing, [-0.01/2, -0.01/6, 0.01/6, 0.01/2])

h_vec = [0.0018, 0.00135, 0.0009, 0.00045, 0.00045/2]
md_vec = []
solmom_vec = []
solfem_vec = []
err_J_mom_vec = []
err_J_fem_vec = []
err_I_mom_vec = []
err_I_fem_vec = []

for h in h_vec

    md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = h)

    # solmom, S, R = IP.solve0(;
    #     md = md, 
    #     material = mat,
    #     κ0 = 1.0,
    #     ϵ0 = nothing,
    #     ω = nothing, 
    #     potential_top = 0.5, 
    #     potential_bottom = -0.5,
    #     qs3D = qs3D, 
    #     qs4D = qs4D, 
    #     qs5D6D = qs5D6D,
    # )

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

    solmom = solfem#!!!!!!!!

    push!(md_vec, md)
    push!(solmom_vec, solmom)
    push!(solfem_vec, solfem)

    

    # Stomdichte
    range_ = range(-0.0049,stop=0.0049,length=20)
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

#dataname = "conv1" # for JLD2 save
#jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; md_vec, solmom_vec, solfem_vec, h_vec, err_J_mom_vec, err_J_fem_vec) 

##



## h-Plot für J (im Volumen)

plt = Plots.plot(h_vec, err_J_mom_vec, line=:scatter, marker = :utriangle, label = "Error J MoM", size=(600,500), legend=:bottomright)#,xaxis=:log2)
plot!(plt,h_vec, err_J_fem_vec, line=:scatter, marker = :dtriangle, label = "Error J FEM")
xlabel!("h in m")
ylabel!("rel. Error")
xlims!(0.0, 0.002)
ylims!(0.0, 0.35)

## h-Plot für I (mittels SWG)

plt = Plots.plot(h_vec, err_I_mom_vec, line=:scatter, marker = :utriangle, label = "Error I MoM", size=(600,500), legend=:topleft)#,xaxis=:log2)
plot!(plt,h_vec, err_I_fem_vec, line=:scatter, marker = :dtriangle, label = "Error I FEM")
xlabel!("h in m")
ylabel!("rel. Error")
xlims!(0.0, 0.002)
ylims!(0.0, 0.3)

## s



##









