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


##

p = Visu.iplot2(; size=(800,600), xticks = ([],[]), yticks = ([],[]), zticks = ([],[]))
p = Visu.mesh(md.Ω ,p)


els, _, _ = assemblydata(md.X)

s1 = els[2]

p = Visu.simplex(p, s1)
plt = Visu.add1(p, s1, refspace(md.X),5e-9)


# Darstellung für PWC?

# Darstellung für LinLag? -> höhenlinien! über ecken!


##






display(p)


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
