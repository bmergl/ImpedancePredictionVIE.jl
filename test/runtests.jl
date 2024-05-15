using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays

using Test

include("test_SS3D.jl")



# Hier muss ein .geo oder ein .msh file im test Ordner liegen
# für den Unit test
# ! Der Nutzer muss später sein eigenes .geo mitbringen
# wobei bestimmte Physical-Regeln eingehalten werden müssen.




#test mesh fehlt noch und auch der MatId_Γ Operator mit den LinLag Elements...
# X = nedelecd3d(Ω)
# I = Identity()
# I2 = IPVIE1.MatId_Ω(alphatype = Float64)
# M_I = assemble(I, X, X)
# M_I2 = assemble(I2, X, X)
# @assert norm(M_I-M_I2) < 1e-15