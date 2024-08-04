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

    swg_faces_fast = swgfaces(md.Ω, md.Γ_nc, fast = true)
    swg_faces_slow = swgfaces(md.Ω, md.Γ_nc, fast = false)
    swg_faces_sa = ImpedancePredictionVIE.swgfaces_set_approach(md.Ω, md.Γ_nc)

    # Andere Reihenfolge beim MultiTheading möglich mache daher normtest
    @test sum(norm.(swg_faces_slow)) ≈ sum(norm.(swg_faces_fast))
    @test sum(norm.(swg_faces_slow)) ≈ sum(norm.(swg_faces_sa))
    
    error("STOP")
end