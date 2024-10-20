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

#md = IP.setup(; meshname = "coarsecube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0033)
md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
println("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

qs3D_high = BEAST.SingleNumQStrat(6)
qs4D_high = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D_high = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)


f_list = [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0] # in Hz
Z_mom_list = []
Z_momhq_list = []
Z_fem_list = []
#Z_ana_list = []
mat = IP.constantmaterial(0.2, 10000.0*IP.ε0) #<------------------------------------------------------- KLÄREN!!!!
lastsol = nothing
Φtop = 0.5
Φbottom = -0.5

for f in f_list

    BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
    BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
    BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D
    solmom, S, R = IP.solve0(;
        md = md, 
        material = mat,
        κ0 = 0.1,#0.1,
        ϵ0 = 1.0*IP.ε0,
        ω = 2*pi*f,  
        potential_top = Φtop, 
        potential_bottom = Φbottom,
        qs3D = qs3D, 
        qs4D = qs4D, 
        qs5D6D = qs5D6D,
    )
    #println("cond mom:")
    #println(cond(S))

    BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D_high
    BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D_high
    BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D_high
    solmom_hq, S, R = IP.solve0(;
        md = md, 
        material = mat,
        κ0 = 0.1,#0.1,
        ϵ0 = 1.0*IP.ε0,
        ω = 2*pi*f,  
        potential_top = Φtop, 
        potential_bottom = Φbottom,
        qs3D = qs3D_high, # das alleine reicht nicht!!! 
        qs4D = qs4D_high, 
        qs5D6D = qs5D6D_high,
    )

    BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
    solfem, S, R = IP.solvefem(;
    md = md, 
    material = mat,
    κ0 = 0.1,
    ϵ0 = 1.0*IP.ε0,
    ω = 2*pi*f,
    potential_top = Φtop, 
    potential_bottom = Φbottom,
    qs3D = qs3D, 
    )
    #println("cond fem:")
    #println(cond(Array(S),2))

    # # Strom durch Platten
    I_top2_mom, I_bottom2_mom = IP.getcurrent2(md, solmom)
    @show norm(I_top2_mom-I_bottom2_mom)# < 0.05

    I_top2_momhq, I_bottom2_momhq = IP.getcurrent2(md, solmom_hq)
    @show norm(I_top2_momhq-I_bottom2_momhq)# < 0.05

    I_top2_fem, I_bottom2_fem = IP.getcurrent2(md, solfem)
    @show norm(I_top2_fem-I_bottom2_fem)# < 0.05

    U =  Φtop - Φbottom

    Zmom = U/I_top2_mom
    Zmomhq = U/I_top2_momhq
    Zfem = U/I_top2_fem
    
    # #Z_ana = U/IP.solution_I_ana(md.body, solmom.material, md, solmom)
    # #push!(Z_ana_list, Z_ana)
    
    push!(Z_mom_list, Zmom)
    push!(Z_momhq_list, Zmomhq)
    push!(Z_fem_list, Zfem)
    
    #lastsol = solmom
end

##


## save

#dataname = "Z_constmat" # for JLD2 save
#jldsave("$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"; f_list, Z_mom_list, Z_momhq_list, Z_fem_list) 

## load 

dataname = "Z_constmat"
datapath = "$(pkgdir(ImpedancePredictionVIE))/data/$dataname.jld2"

f_list = load(datapath, "f_list")
Z_mom_list = load(datapath, "Z_mom_list")
Z_momhq_list = load(datapath, "Z_momhq_list") 
Z_fem_list = load(datapath, "Z_fem_list") 


##

Z_ana_list_ext = []

L_z = md.body.L_z
L_x = md.body.L_x
L_y = md.body.L_y   
ϵ = mat.ϵ
κ = mat.κ
U = Φtop - Φbottom
R = (1/κ)*L_z/(L_x*L_y)
C = ϵ*(L_x*L_y)/L_z

#f_list_ext = collect(range(f_list[1], stop = f_list[end], step = 80.0))
#f_list_ext = collect(range(1.0, stop = 10^7, step = 100.0))
f_list_ext = 10 .^ range(0, 7, length=1000)


for f in f_list_ext
    ω = 2*pi*f
    Z_ = (R - im*ω*R^2*C)/(1 + ω^2*R^2*C^2)
    push!(Z_ana_list_ext,Z_)
end



## Plot: |Z|

plt = Plots.plot(f_list_ext, abs.(Z_ana_list_ext), color = 1, line=:line , label = "                                                ", xaxis=:log10)
plot!(plt, f_list, abs.(Z_mom_list), color = 2, line=:scatter, marker = :utriangle, label = "                                    ")
plot!(plt, f_list, abs.(Z_momhq_list), color = 10,  line=:scatter, marker = :dtriangle, label = "                                    ")
plot!(plt, f_list, abs.(Z_fem_list), color = 3, line=:scatter, marker = :circle, markersize = 2.7, label = "                                    ")
plot!(plt, size=(400,300))
plot!(plt, legend=:bottomleft)
#xlabel!("f in Hz")
#ylabel!("Impedance in Ω")
ylims!(450.0,510.0)


## Plot: arg(Z)

plt = Plots.plot(f_list_ext, angle.(Z_ana_list_ext)/(2*pi)*360, color = 1, line=:line , label = "                                                ", xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, angle.(Z_mom_list)/(2*pi)*360, color = 2, line=:scatter, marker = :utriangle, label = " ")
plot!(plt, f_list, angle.(Z_momhq_list)/(2*pi)*360, color = 10, line=:scatter, marker = :dtriangle, label = " ")
plot!(plt, f_list, angle.(Z_fem_list)/(2*pi)*360, color = 3, line=:scatter, marker = :circle, markersize = 2.7, label = " ")
plot!(plt, size=(400,300))
#xlabel!("f in Hz")
#ylabel!("Angle in °")
ylims!(-90.0, 5.0)
ylims!(-1.0, 0.1)


##

## Plot: Re

plt = Plots.plot(f_list_ext, real.(Z_ana_list_ext), line=:line, label = "Re{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, real.(Z_mom_list), line=:scatter, label = "Re{Z_mom}")
plot!(plt, f_list, real.(Z_momhq_list), line=:scatter, label = "Re{Z_momhq}")
plot!(plt, f_list, real.(Z_fem_list), line=:scatter, label = "Re{Z_fem}")

xlabel!("f in Hz")
ylabel!("Impedance in Ω")
#ylims!(0.0,510.0)

## Plot: Im

plt = Plots.plot(f_list_ext, imag.(Z_ana_list_ext), line=:line , label = "Im{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, imag.(Z_mom_list), line=:scatter, label = "Im{Z_mom}")
plot!(plt, f_list, imag.(Z_momhq_list), line=:scatter, label = "Im{Z_momhq}")
plot!(plt, f_list, imag.(Z_fem_list), line=:scatter, label = "Im{Z_fem}")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")









## speichern

x_values = f_list_ext
y_values = angle.(Z_ana_list_ext)/(2*pi)*360

# Öffne eine Datei zum Schreiben
open("results/argZ_ana.txt", "w") do file
    for i in eachindex(x_values)
        # Schreibe die Werte von x und y nebeneinander
        # \t fügt einen Tabulator zwischen die Werte
        write(file, "$(x_values[i])\t$(y_values[i])\n")
    end
end

