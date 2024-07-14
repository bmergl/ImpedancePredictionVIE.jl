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

include("analytic_solution.jl")

include("geometry.jl")


const ε0 = 8.8541878188*1e-12

export Visu

export IPVIE
export IPVIE1
export IPVIE2
export IPVIE3
export geo2mesh
export showxvec
export MaterialIdentity

export realvertices
export realnodes 
export pointlist2xyzlist
export swgfaces
export MaterialIdentity
export gen_tau_chi
export getcurrent

export kernelvalsdyad



end





