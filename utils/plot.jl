using Plots

Plots.savefig("x.svg")
Plots.savefig(p1,"x")

plot!(x, -J_z, label = "J_z")


