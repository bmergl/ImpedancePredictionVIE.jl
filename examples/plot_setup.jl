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


md = IP.setup(; meshname = "coarsecube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0033)
print("tehrahedrons: ", length(md.Ω.faces))


## Volumenmesh #######################################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Ω, p, nodesize = 0.02, nodecolor = "blue", linewidth = 1.0, linecolor = "green")
p = Visu.mesh(md.Γ_c, p, nodesize = 0.02, nodecolor = "blue", linewidth = 1.0, linecolor = "blue")



## SWG ###############################################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Ω, p, nodesize = 0.02, nodecolor = "blue", linewidth = 0.5, linecolor = "green")
p = Visu.mesh(md.Γ_c, p, nodesize = 0.02, nodecolor = "blue", linewidth = 0.5, linecolor = "blue")


els, _, _ = assemblydata(md.X)


tetnr1 = 120 # 120 geht!
s1 = els[tetnr1]
c1 = cartesian(CompScienceMeshes.center(s1))
distances = Vector{Float64}()
inds = Vector{Int64}()
for (j,el) in enumerate(els)
    j == tetnr1 && continue
    cj = cartesian(CompScienceMeshes.center(el))
    push!(distances, norm(cj-c1))  
    push!(inds, j)
end
j_min = findmin(distances)[2]
tetnr2 = inds[j_min] 
s2 = els[tetnr2]

r = rand(3)

p = Visu.simplex(p, s1; linewidth = 3.0, color = "black")
p = Visu.simplex(p, s2; linewidth = 3.0, color = "black")
p = Visu.add1(p, s1, refspace(md.X), 3e-9, 4, 1.2, r)
p = Visu.add1(p, s2, refspace(md.X), 3e-9, 3, -1.2, r)


stetnr1 = 157 #130 #110
#stetnr2 = 30
t1 = els[stetnr1]
#t2 = els[stetnr2]
p = Visu.simplex(p, t1; linewidth = 3.0, color = "black")
#p = Visu.simplex(p, t2; linewidth = 2.0, color = "black")
p = Visu.add1(p, t1, refspace(md.X), 3e-9, 3, 0.7, r)
#p = Visu.add1(p, t2, refspace(md.X), 3e-9, 3, -1.0)





## LinearLagrange #######################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Γ_nc, p, nodesize = 1.0, nodecolor = "green", linewidth = 0.8, linecolor = "green")
p = Visu.mesh(md.Γ_c, p, nodesize = 1.0, nodecolor = "blue", linewidth = 0.8, linecolor = "blue")
#p = Visu.fnspos(md.y, p, color = 1, markersize = 1.0)

y = lagrangec0d1(md.Γ)
els, ad, nr = assemblydata(y)

node = 17
point_ = y.pos[node]
shs_list = y.fns[node]
cell_numbers = Vector{Int64}()
for shs in shs_list
    push!(cell_numbers, shs.cellid)
end
s_list = [els[i] for i in cell_numbers]
for s in s_list
    p = Visu.simplex(p, s; linewidth = 6.0, color = "blue")
end


node = 61
point_ = y.pos[node]
shs_list = y.fns[node]
cell_numbers = Vector{Int64}()
for shs in shs_list
    push!(cell_numbers, shs.cellid)
end
s_list = [els[i] for i in cell_numbers]
for s in s_list
    p = Visu.simplex(p, s; linewidth = 6.0, color = "green")
end

display(p)



## PWC ####################################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Γ_nc, p, nodesize = 0.01, nodecolor = "green", linewidth = 0.8, linecolor = "green")
p = Visu.mesh(md.Γ_c, p, nodesize = 0.01, nodecolor = "blue", linewidth = 0.8, linecolor = "blue")

y = lagrangec0d1(md.Γ) #ja...
els, _, _ = assemblydata(y)
s1 = els[20]
p = Visu.simplex(p, s1; linewidth = 6.0, color = "red")
s2 = els[104]
p = Visu.simplex(p, s2; linewidth = 6.0, color = "black")



## LinLag excitation ###########################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))

els, _, _ = assemblydata(md.y)
s1 = els[20]

p = Visu.mesh(md.Γ, p)
p = Visu.simplex(p, s1)


## LinLag no excitation ############################################

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))

els, _, _ = assemblydata(md.y)

p = Visu.mesh(md.Γ_nc, p)

node = 30
point_ = md.y.pos[node]
shs_list = md.y.fns[node]
cell_numbers = Vector{Int64}()
for shs in shs_list
    push!(cell_numbers, shs.cellid)
end
s_list = els[cell_numbers]

for s in s_list
    p = Visu.simplex(p, s; linewidth = 3.0, color = "black")
end

display(p)







##



els_tree = BEAST.octree(els)
j = CompScienceMeshes.findchart(els, els_tree, SVector(0,0,0))


## PCW (Achtung, Zahlen später ergänzen: 1 in der Mitte)

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))

els, _, _ = assemblydata(md.w)
s1 = els[20]

p = Visu.mesh(md.Γ_c, p)
p = Visu.simplex(p, s1)


# Darstellung für PWC?

# Darstellung für LinLag? -> höhenlinien! über ecken!


##


# Sample data
x = 1:10
y = rand(10)

# Create a plot with Plotly backend
p = Plots.plot(x, y)

# Modify the axes properties to remove numbers
plot!(xticks = ([],[]), yticks = ([],[]))



# Sample data
x = 1:10
y = rand(10)

# Create a plot with Plotly backend
p = Plots.plot(x, y, backend = plotly())

# Customize the axes to remove numerical labels but keep ticks
p.xticks = (p.xticks[1], fill("", length(p.xticks[1])))
p.yticks = (p.yticks[1], fill("", length(p.yticks[1])))

# Display the plot
display(p)

##

