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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
println("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 

#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)
#qs3D = BEAST.SingleNumQStrat(6)
#qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
#qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)
BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


f_list = [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0] # in Hz
Z_mom_list = []
Z_fem_list = []
Z_ana_list = []
mat = IP.constantmaterial(0.2, 10000.0*IP.ε0) #<------------------------------------------------------- KLÄREN!!!!
lastsol = nothing

for f in f_list

    solmom, S, R = IP.solve0(;
        md = md, 
        material = mat,
        κ0 = 0.1,
        ϵ0 = 1.0*IP.ε0,
        ω = 2*pi*f,  
        potential_top = 0.5, 
        potential_bottom = -0.5,
        qs3D = qs3D, 
        qs4D = qs4D, 
        qs5D6D = qs5D6D,
    )

    error("UNTEN QUADSTRAT ERHÖHEN!!!!")
    solmom_hq, S, R = IP.solve0(;
        md = md, 
        material = mat,
        κ0 = 0.1,
        ϵ0 = 1.0*IP.ε0,
        ω = 2*pi*f,  
        potential_top = 0.5, 
        potential_bottom = -0.5,
        qs3D = qs3D, 
        qs4D = qs4D, 
        qs5D6D = qs5D6D,
    )

    solfem, S, R = IP.solvefem(;
    md = md, 
    material = mat,
    κ0 = 0.1,
    ϵ0 = 1.0*IP.ε0,
    ω = 2*pi*f,
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    )


    # Strom durch Platten
    I_top2_mom, I_bottom2_mom = IP.getcurrent2(md, solmom)
    @show norm(I_top2_mom-I_bottom2_mom)# < 0.05
    I_top2_fem, I_bottom2_fem = IP.getcurrent2(md, solfem)
    @show norm(I_top2_fem-I_bottom2_fem)# < 0.05


    U = solmom.potential_top - solmom.potential_bottom

    Zmom = U/I_top2_mom
    Zfem = U/I_top2_fem
    
    Z_ana = U/IP.solution_I_ana(md.body, solmom.material, md, solmom)

    push!(Z_mom_list, Zmom)
    push!(Z_fem_list, Zfem)
    push!(Z_ana_list, Z_ana)


    lastsol = solmom
end

##

##




Z_ana_list_ext = []

L_z = md.body.L_z
L_x = md.body.L_x
L_y = md.body.L_y   
ϵ = lastsol.material.ϵ
κ = lastsol.material.κ
U = lastsol.potential_top-lastsol.potential_bottom
R = (1/κ)*L_z/(L_x*L_y)
C = ϵ*(L_x*L_y)/L_z

f_list_ext = collect(range(f_list[1], stop = f_list[end], step = 100.0))

for f in f_list_ext
    ω = 2*pi*f
    Z_ = (R - im*ω*R^2*C)/(1 + ω^2*R^2*C^2)
    push!(Z_ana_list_ext,Z_)
end




## 4 Kurven !!!!!!!!!!!!!!!!!!!!!!!!!!!


## Plot: Re

plt = Plots.plot(f_list_ext, real.(Z_ana_list_ext), line=:line, label = "Re{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, real.(Z_mom_list), line=:scatter, label = "Re{Z_mom}")
plot!(plt, f_list, real.(Z_fem_list), line=:scatter, label = "Re{Z_fem}")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")
#ylims!(0.0,510.0)

## Plot: Im

plt = Plots.plot(f_list_ext, imag.(Z_ana_list_ext), line=:line , label = "Im{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, imag.(Z_mom_list), line=:scatter, label = "Im{Z_mom}")
plot!(plt, f_list, imag.(Z_fem_list), line=:scatter, label = "Im{Z_fem}")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")


## Plot: |Z|

plt = Plots.plot(f_list_ext, abs.(Z_ana_list_ext), line=:line , label = "|Z_ana|", size=(600,500), xaxis=:log10)
plot!(plt, f_list, abs.(Z_mom_list), line=:scatter, label = "|Z_mom|")
plot!(plt, f_list, abs.(Z_fem_list), line=:scatter, label = "|Z_fem|")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")
#ylims!(0.0,510.0)


## Plot: arg(Z)

plt = Plots.plot(f_list_ext, angle.(Z_ana_list_ext)/(2*pi)*360, line=:line , label = "arg{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, angle.(Z_mom_list)/(2*pi)*360, line=:scatter, marker = :utriangle, label = "arg{Z_mom}")
plot!(plt, f_list, angle.(Z_fem_list)/(2*pi)*360, line=:scatter, marker = :dtriangle, label = "arg{Z_fem}")
xlabel!("f in Hz")
ylabel!("Angle in °")
ylims!(-90.0, 5.0)









