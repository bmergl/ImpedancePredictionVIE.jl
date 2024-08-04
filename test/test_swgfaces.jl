using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays

using Test

@testset "swg faces" begin

    md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(0.01, 0.01, 0.01), h = 0.0018)
    print("tehrahedrons: ", length(md.Ω.faces))

    swg_faces = swgfaces(md.Ω, md.Γ_nc, fast = true)
    swg_faces2 = swgfaces(md.Ω, md.Γ_nc, fast = false)
    swg_faces3 = ImpedancePredictionVIE.swgfaces2(md.Ω, md.Γ_nc)

    # Andere Reihenfolge beim MultiTheading möglich mache daher normtest
    @show sum(norm.(swg_faces))
    @show sum(norm.(swg_faces2))
    @test sum(norm.(swg_faces)) ≈ sum(norm.(swg_faces2))
    @test sum(norm.(swg_faces2)) ≈ sum(norm.(swg_faces3))
    
    error("STOP")
end