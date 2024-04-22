module ImpedancePredictionVIE

using BEAST
using CompScienceMeshes
using LinearAlgebra
using StaticArrays
using Plots


include("plotdata.jl")

export visu
export hello1
export showxvec

function hello1(a::F) where{F}
    return a + F(1.0)
end

function showxvec()
    display(x̂)
    return
end




# Write your package code here.

end
