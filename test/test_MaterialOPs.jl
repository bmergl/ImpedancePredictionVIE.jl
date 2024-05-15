using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays

using Test

@testset "Material Operators" begin

    geoname = "cube.geo"
    geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"

    meshname = "cube.msh"
    meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"

    h = 2.0
    Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h)

    y1 = lagrangec0d1(Γ) # nur allgemeines testen...

    gamma = 0.0
    alpha = 1.0
    tau = x -> 1.0

    MSL = ImpedancePredictionVIE.MaterialSL(gamma, alpha, tau)
    SL = Helmholtz3D.singlelayer(gamma=gamma, alpha = alpha)
    @test norm(assemble(MSL,y1,y1) - assemble(SL,y1,y1)) < 1e-15

    MDL = ImpedancePredictionVIE.MaterialDL(gamma, alpha, tau)
    DL = Helmholtz3D.doublelayer(gamma=gamma, alpha = alpha)
    @test norm(assemble(MDL,y1,y1) - assemble(DL,y1,y1)) < 1e-15

    MADL = ImpedancePredictionVIE.MaterialADL(gamma, alpha, tau)
    ADL = Helmholtz3D.doublelayer_transposed(gamma=gamma, alpha = alpha)
    @test norm(assemble(MADL,y1,y1) - assemble(ADL,y1,y1)) < 1e-15

    MId = ImpedancePredictionVIE.MatId(alpha, tau)
    Id = Identity()
    @test norm(assemble(MId,y1,y1) - assemble(Id,y1,y1)) < 1e-15

end


