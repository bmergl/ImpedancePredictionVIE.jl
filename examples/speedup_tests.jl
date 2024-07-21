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

Mh = hassemble(B33_ΩΩ, md.X, md.X;
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