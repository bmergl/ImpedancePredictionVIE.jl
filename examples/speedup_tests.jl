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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Ω) 


##


#Quadstrat
qs3D = BEAST.SingleNumQStrat(4)
qs4D = BEAST.DoubleNumSauterQstrat(3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D

# STANDARD-TESTMATERIAL: IP.pwlinx([[1.0, 2000.0],[4000.0, 10000.0],[20000.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2])

sol, S, R = IP.solve0(;   # solve -> arb. Mat. / solve1 -> high contrast formulation
    md = md, 
    material = IP.constantmaterial(2.0, nothing), # IP.pwlinx([[1.0, 1.2],[1.4, 1.7],[1.9, 1.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), #IP.constant_zsplit(100.0, nothing, 0.0001, 10.0, nothing), #IP.pwlinx([[1.0, 20.0],[40.0, 100.0],[200.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]), # # IP.pwlinx([[1.0, 2.0],[4.0, 10.0],[20.0, 5.0]], nothing, [-md.body.L_x/2, -0.01/6, 0.01/6, md.body.L_x/2]),   #IP.general_material(κ, nothing),  #  IP.constant_xsplit(0.13, nothing, 0.0, 0.00007, nothing), #IP.constant_zsplit(10.0, nothing, 0.0, 0.001, nothing),, #  ,#, # #, #
    κ0 = 1.0, # möglichst in der nähe der realen Größen wählen damit cond(S) klein?
    ϵ0 = nothing,
    ω = nothing, 
    potential_top = 0.5, 
    potential_bottom = -0.5,
    qs3D = qs3D, 
    qs4D = qs4D, 
    qs5D6D = qs5D6D,
    #matalloc = :avg#:center,
)

## It solver test: STOP SINNLOS OHNE PRECOND... #############

# Zeilenskalierungsfaktoren berechnen
row_scale_factors = maximum(abs, S, dims=2)
D_r = Diagonal(Vector(1.0 ./ row_scale_factors[:,1]))

# Spaltenskalierungsfaktoren berechnen
col_scale_factors = maximum(abs, S, dims=1)
D_c = Diagonal(1.0 ./ col_scale_factors[1,:])

# Beidseitige Skalierung der Matrix
S_scaled = D_r * S * D_c

# Skaliere die RHS
b_scaled = D_r * sol.b

# Löse das skalierte System
y = S_scaled \ b_scaled

# Transformiere die Lösung zurück
x = D_c * y

sol.u
x

norm(x-sol.u)

cond(S)
cond(S_scaled)


b_test = ones(length(b_scaled))

u_it, ch = IterativeSolvers.gmres!(x, S_scaled, b_test; maxiter=4000, log=true, reltol = 1.0e-3)
rank(S_scaled)
S_scaled

@show norm(S*u_it-b)
@show ch.isconverged
@show norm(sol.u-u_it)

ilu = factorize(S_scaled)
result_ilu = IterativeSolvers.gmres(S_scaled, b_scaled, Pl=ilu.L, Pr=ilu.U, maxiter=5000, log=true)


eigenvalues = eigvals(S)




## It solver test: STOP SINNLOS OHNE PRECOND... #############
n = 1000
A= 100*LinearAlgebra.I(n) +  randn(n,n)
cond(A)
b = randn(n)
b = vcat([100],zeros(99))
x = IterativeSolvers.gmres(A, b; maxiter=1000)
@show norm(A*x-b)

u_it, ch = IterativeSolvers.gmres(S, sol.b; maxiter=10000, log=true)
@show norm(sol.u-u_it)
@show ch


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



## SS3D speedup 
using SauterSchwab3D

function k6p_cf_A_faster(f, x1::T, x2::T, x3::T, y1::T, y2::T, z1::T) where T
    #sum = 0;
    
    #Domain 1
    #=1=#   ξ1 = x1
    #=2=#   ξ2 = x2
    #=3=#   ξ3 = x3
    #=4=#   η1 = x3 * y1
    #=5=#   η2 = x3 * y2
    #=6=#   η3 = x3 * z1

    det = x3^3

    sum = f((1-ξ1,ξ1-ξ2,ξ3-η3),(1-(ξ1-ξ3+η1),ξ1-ξ2+η2,η1-η2))*det

    return sum
end 

    
n = 300000000
function f(u,v) 
    return u[1]-u[2]+u[3]-v[1]+v[2]-v[3]#  x1^2 + x2 - 10*x3
end

@time for i in 1:n
    r6 = randn(6)
    SauterSchwab3D.k6p_cf_A(f, r6...)
end

@time for i in 1:n
    r6 = randn(6)
    k6p_cf_A_faster(f, r6...)
end

