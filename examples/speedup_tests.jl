## It solver test: STOP SINNLOS OHNE PRECOND... #############
# n = 1000
# A= 100*LinearAlgebra.I(n) +  randn(n,n)
# cond(A)
# b = randn(n)
#b = vcat([100],zeros(99))
#x = IterativeSolvers.gmres(A, b; maxiter=1000)
#@show norm(A*x-b)

# u_it, ch = IterativeSolvers.gmres(S, sol.b; maxiter=10000, log=true)
# @show norm(sol.u-u_it)
# @show ch


## H-Matrix #######################################

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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0007)
print("tehrahedrons: ", length(md.Ω.faces))

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



B33_ΩΩ = IPVIE.B33_ΩΩ(alpha = 1.0, gammatype = Float64, chi = x->1.0)


 
Mh = hassemble(B33_ΩΩ, md.X, md.X;
    treeoptions=BoxTreeOptions(nmin=100),
    compressor = FastBEAST.ACAOptions(tol=1e-4),
    #quadstratcbk=BEAST.defaultquadstrat(operator, test_functions, trial_functions),
    #quadstratfbk=BEAST.defaultquadstrat(operator, test_functions, trial_functions),
    multithreading=true,
    verbose=true
)


@time M = assemble(B33_ΩΩ, md.X, md.X)

println("done!")
#@show norm(M-Matrix(Mh)) # Matrix(Mh) geht gar nicht...

##