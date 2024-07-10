module ImpedancePredictionVIE

using LinearAlgebra
using StaticArrays
using SparseArrays
using BEAST
using CompScienceMeshes
using SauterSchwab3D

using JLD2

using Plots


IP = ImpedancePredictionVIE
export IP

include("ip.jl")

include("visu.jl")

include("vieops.jl")
include("ipvie.jl")
include("ipvie1.jl")
include("ipvie2.jl")
include("ipvie3.jl")

include("geometry.jl")



export Visu

export IPVIE
export IPVIE1
export IPVIE2
export IPVIE3
export geo2mesh
#export draw_arrow!
export showxvec
export MaterialIdentity

export realvertices
export realnodes 
export pointlist2xyzlist
export swgfaces
#export dirichlet_n2f
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




