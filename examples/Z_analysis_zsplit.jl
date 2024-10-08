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
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 

#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)
# qs3D = BEAST.SingleNumQStrat(6)
# qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
# qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)
BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


f_list = [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0] # in Hz
Z_list = []
Z2_list = []
Z_ana_list = []
for f in f_list
    sol, S, R = IP.solve0(;
        md = md, 
        material = IP.constant_zsplit(0.5, 3.0*IP.ε0, 0.0008, 0.2, 10.0*IP.ε0),
        κ0 = 0.1, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
        ϵ0 = 1.0*IP.ε0, # Größenordnung??? was ist sinnvoll?
        ω = 2*pi*f, 
        potential_top = 0.5, 
        potential_bottom = -0.5,
        qs3D = qs3D, 
        qs4D = qs4D, 
        qs5D6D = qs5D6D 
    )
    # Strom durch Platten
    I_top, I_bottom = IP.getcurrent(md, sol)
    @assert norm(I_top-I_bottom) < 0.02
    I_top2, I_bottom2 = IP.getcurrent2(md, sol)
    @assert norm(I_top2-I_bottom2) < 0.02
    U = sol.potential_top - sol.potential_bottom
    Z = U/I_top
    Z2 = U/I_top2
    Z_ana = U/IP.solution_I_ana(md.body, sol.material, md, sol)
    push!(Z_list, Z)
    push!(Z2_list, Z2) 
    push!(Z_ana_list, Z_ana)
end
##

##
## Plot: Re

plt = Plots.plot(f_list, real.(Z_ana_list), line=:line, label = "Re{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, real.(Z_list), line=:scatter, label = "Re{Z}")
plot!(plt, f_list, real.(Z2_list), line=:scatter, label = "Re{Z2}")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")
ylims!(0.0,510.0)

## Plot: Im

plt = Plots.plot(f_list_ext, imag.(Z_ana), line=:line , label = "Im{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, imag.(Z_list), line=:scatter, label = "Im{Z}")
plot!(plt, f_list, imag.(Z2_list), line=:scatter, label = "Im{Z2}")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")
#ylims!(-50.0,0.0)

## Plot: |Z|

plt = Plots.plot(f_list, abs.(Z_ana_list), line=:line , label = "|Z_ana|", size=(600,500), xaxis=:log10)
plot!(plt, f_list, abs.(Z_list), line=:scatter, label = "|Z|")
plot!(plt, f_list, abs.(Z2_list), line=:scatter, label = "|Z2|")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")
#ylims!(0.0,510.0)


## Plot: arg(Z)

plt = Plots.plot(f_list_ext, angle.(Z_ana)/(2*pi)*360, line=:line , label = "arg{Z_ana}", size=(600,500), xaxis=:log10, legend=:bottomleft)
plot!(plt, f_list, angle.(Z_list)/(2*pi)*360, line=:scatter, marker = :utriangle, label = "arg{Z}")
plot!(plt, f_list, angle.(Z2_list)/(2*pi)*360, line=:scatter, marker = :dtriangle, label = "arg{Z2}")
xlabel!("f in Hz")
ylabel!("Angle in °")
ylims!(-90.0, 5.0)










## Plot: Real - Imag

plt = Plots.plot(f_list, real.(Z_list), line=:line , label = "Re{Z}", size=(700,600), xaxis=:log10)
plot!(plt, f_list, imag.(Z_list), label = "Im{Z}")
plot!(plt, f_list, real.(Z2_list), label = "Re{Z2}")
plot!(plt, f_list, imag.(Z2_list), label = "Im{Z2}")
plot!(plt, f_list, real.(Z_ana_list), label = "Re{Z_ana}")
plot!(plt, f_list, imag.(Z_ana_list), label = "Im{Z_ana}")
#xlims!(0.0, 1.0)
#ylims!(-0.01, 0.0)
#title!("J_z(x, y0, z0)")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")


## Plot: Ab

plt = Plots.plot(f_list, abs.(Z_list), line=:line , label = "|Z|", size=(700,600), xaxis=:log10)
plot!(plt, f_list, abs.(Z2_list), label = "|Z2|")
plot!(plt, f_list, abs.(Z_ana_list), label = "|Z_ana|")

#xlims!(0.0, 1.0)
ylims!(0.0, 600.0)
#title!("J_z(x, y0, z0)")
xlabel!("f in Hz")
ylabel!("Impedance in Ω")

##

plt = Plots.plot(f_list, angle.(Z_list)/(2*pi)*360, label = "arg{Z}", line=:line , size=(700,600), xaxis=:log10)
plot!(plt, f_list, angle.(Z2_list)/(2*pi)*360, label = "arg{Z2}")
plot!(plt, f_list, angle.(Z_ana_list)/(2*pi)*360, label = "arg{Z_ana}")
xlabel!("f in Hz")
ylabel!("Angle in °")

