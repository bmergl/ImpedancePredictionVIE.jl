module ImpedancePredictionVIE

using LinearAlgebra
using StaticArrays
using BEAST
using CompScienceMeshes
using SauterSchwab3D

using JLD2

using Plots



include("ip.jl")

include("visu.jl")

include("integrals.jl")

include("vieops1.jl")
include("ipvie1.jl")

include("vieops2.jl")
include("ipvie2.jl")


include("geometry.jl")



export Visu
export IP
export IPVIE1
export IPVIE2
export geo2mesh
#export draw_arrow!
export hello1
export showxvec
export MaterialIdentity

export realvertices
export realnodes 
export pointlist2xyzlist
export swgfaces
export gen_tau_invtau
#export dirichlet_n2f
export BoundaryOperatorΩΓ
export BoundaryOperatorΓΩ
export MaterialIdentity
export gen_tau_chi
export getcurrent

export kernelvalsdyad


# struct sample
#     Ω
#     \G

# end

end




# function __init__()

#     abstract type MaterialIdentity <: BEAST.LocalOperator end

#     scalartype(localop::MaterialIdentity) = typeof(localop.α)

# end




