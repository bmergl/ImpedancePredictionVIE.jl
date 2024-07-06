using MKL

using LinearAlgebra
using StaticArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.18)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Γ) 

y = md.y
w = md.w
X = md.X

ntrc = md.ntrc

##

# Operators row 1
B11_Γ = IPVIE1.B11_Γ()
assemble(B11_Γ, w, y)
B11_ΓΓ = IPVIE1.B11_ΓΓ(gammatype = Float64)
assemble(B11_ΓΓ, w, y)

UB12_ΓΓ = IPVIE1.UB12_ΓΓ(gammatype = Float64)
assemble(UB12_ΓΓ, w, w)

UB13_ΓΓn = IPVIE1.UB13_ΓΓn(gammatype = Float64)
#assemble(UB13_ΓΓn, w, NEW SPACE NEEDED!!!!)
UB13_ΓΩ = IPVIE1.UB13_ΓΩ(gammatype = Float64)
assemble(UB13_ΓΩ, w, X)


# Operators row 2
B21_ΓΓ = IPVIE1.B21_ΓΓ(gammatype = Float64)
assemble(B21_ΓΓ, y, y)

UB22_Γ = IPVIE1.UB22_Γ()
assemble(UB22_Γ, y, w)
UB22_ΓΓ = IPVIE1.UB22_ΓΓ(gammatype = Float64)
assemble(UB22_ΓΓ, y, w)

UB23_ΓΓn = IPVIE1.UB23_ΓΓn(gammatype = Float64)
#assemble(B23_ΓΓ, y, NEW SPACE NEEDED!!!!)
UB23_ΓΩ = IPVIE1.UB23_ΓΩ(gammatype = Float64)
assemble(UB23_ΓΩ, y, X)


# Operators row 3
B31_ΓΓ = IPVIE1.B31_ΓΓ(gammatype = Float64)
assemble(B31_ΓΓ, ntrc(X), y)
B31_ΩΓ = IPVIE1.B31_ΩΓ(gammatype = Float64)
assemble(B31_ΩΓ, X, y)

UB32_ΓΓ = IPVIE1.UB32_ΓΓ(gammatype = Float64)
assemble(UB32_ΓΓ, ntrc(X), w)
UB32_ΩΓ = IPVIE1.UB32_ΩΓ(gammatype = Float64)
assemble(UB32_ΩΓ, X, w)

UB33_Ω = IPVIE1.UB33_Ω()
assemble(UB33_Ω, X, X)
UB33_ΓΓn = IPVIE1.UB33_ΓΓn(gammatype = Float64)
#assemble(UB33_ΓΓ, ntrc(X), NEW SPACE NEEDED!!!!)
UB33_ΓΩ = IPVIE1.UB33_ΓΩ(gammatype = Float64)
assemble(UB33_ΓΩ, ntrc(X), X)
UB33_ΩΓn = IPVIE1.UB33_ΩΓn(gammatype = Float64)
#assemble(B33_ΩΓ, X, NEW SPACE NEEDED!!!!)
UB33_ΩΩ = IPVIE1.UB33_ΩΩ(gammatype = Float64)
assemble(UB33_ΩΩ, X, X)



# Schritt 1: ntrace artige function mit input











##

elements, ad, cells = assemblydata(md.X)

ad[4]

cells[end]
ad[1][1]

for (n,b) in ad[1000][1]
    @show n, b
end
xy
md.X.fns
# erstelle X_mat basis für trial in den speziellen Fällen!



mutable struct testop
    mf::Vector{Float64}
    v::Float64
end

a = Vector([0.1, 0.9, 1.3])
v_ini = 0.0
op1 = testop(a,v_ini)

function change!(op)
    #new_op = testop(op.mf,50.8)
    #op = new_op
    op.v = 50.8
    return nothing
end

change!(op1)
op1







# Quadstrat
qs3D = BEAST.SingleNumQStrat(6)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D