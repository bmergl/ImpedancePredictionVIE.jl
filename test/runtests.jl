using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays

using Test

@testset "ImpedancePredictionVIE.jl" begin
    # Write your tests here.

    geoname = "cube.geo"
    geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"

    meshname = "cube.msh"
    meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"

    h = 2.0 # kleiner 0.2 sonst std
    Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h)


    function make_T4(;gamma,alpha = 1.0,tau = nothing)# umbenennen!
        tau === nothing && error()
        return MaterialDL(gamma,alpha,tau)
    end


end



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