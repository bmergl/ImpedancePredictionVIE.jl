using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using SphericalScattering
using Test
using ImpedancePredictionVIE

#@testset "SauterSchwab 5D/6D integrals" begin


    # Neudefinition ... unklar wie es anders gehen würde...
    struct VIEhhVolume2{T,U,P} <: ImpedancePredictionVIE.VolumeOperatorΩΩ#BEAST.VolumeOperator
        gamma::T
        α::U
        tau::P
    end
    struct VIEhhBoundary2{T,U,P} <: ImpedancePredictionVIE.BoundaryOperatorΓΩ
        gamma::T
        α::U
        tau::P
    end
    struct VIEhhVolumegradG2{T,U,P} <: ImpedancePredictionVIE.VolumeOperatorΩΩ #BEAST.VolumeOperator
        gamma::T
        α::U
        tau::P
    end

    # Konstruktoren
    function hhvolume2(; tau=nothing)

        return VIEhhVolume2(0.0, -1.0, tau)
    end

    function hhboundary2(; tau=nothing)

        return VIEhhBoundary2(0.0, 1.0, tau)
    end

    function hhvolumegradG2(; tau=nothing)

        return VIEhhVolumegradG2(0.0, -1.0, tau)
    end

    
    function BEAST.integrand(viop::VIEhhVolume2, kerneldata, tvals, tgeo, bvals, bgeo)

        gx = @SVector[tvals[i].value for i in 1:4]
        fy = @SVector[bvals[i].value for i in 1:4]

        dgx = @SVector[tvals[i].gradient for i in 1:4]
        dfy = @SVector[bvals[i].gradient for i in 1:4]

        G = kerneldata.green
        gradG = kerneldata.gradgreen

        Ty = kerneldata.tau

        α = viop.α

        return @SMatrix[α * dot(dgx[i], G*Ty*dfy[j]) for i in 1:4, j in 1:4]
    end

    function BEAST.integrand(viop::VIEhhBoundary2, kerneldata, tvals, tgeo, bvals, bgeo)

        gx = @SVector[tvals[i].value for i in 1:3]
        dfy = @SVector[bvals[i].gradient for i in 1:4]

        G = kerneldata.green
        gradG = kerneldata.gradgreen

        Ty = kerneldata.tau

        α = viop.α

        return @SMatrix[α * dot( tgeo.patch.normals[1]*gx[i],G*Ty*dfy[j]) for i in 1:3, j in 1:4]
    end

    function BEAST.integrand(viop::VIEhhVolumegradG2, kerneldata, tvals, tgeo, bvals, bgeo)

        gx = @SVector[tvals[i].value for i in 1:4]
        fy = @SVector[bvals[i].value for i in 1:4]

        dgx = @SVector[tvals[i].gradient for i in 1:4]
        dfy = @SVector[bvals[i].gradient for i in 1:4]

        G = kerneldata.green
        gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

        Ty = kerneldata.tau

        α = viop.α

        return @SMatrix[α * gx[i] * dot(gradG, Ty*dfy[j]) for i in 1:4, j in 1:4]
    end








    # Homogeneous Dielectic Sphere Unit Test from lsvie

    # Environment
    ε1 = 1.0*ε0
    μ1 = μ0

    # Dielectic Sphere
    ε2 = 5.0*ε0
    μ2 = μ0
    r = 1.0
    sp = DielectricSphere(; radius = r, filling = Medium(ε2, μ2))


    # Mesh, Basis
    h = 0.35
    mesh = CompScienceMeshes.tetmeshsphere(r,h)
    bnd = boundary(mesh)
    X = lagrangec0d1(mesh; dirichlet = false)
    #@show numfunctions(X)
    strc = X -> strace(X,bnd)

    # VIE operators
    function generate_tau(ε_ins, ε_env)
        contr = 1.0 - ε_ins/ε_env
        function tau(x::SVector{U,T}) where {U,T}
            return T(contr)
        end
        return tau
    end
    τ = generate_tau(ε2, ε1)

    I = Identity()
    V = hhvolume2(tau = τ)
    B = hhboundary2(tau = τ)
    Y = hhvolumegradG2(tau = τ)

    # I, V, B =  Identity(), VIE.hhvolume(tau = τ, wavenumber = 0.0), VIE.hhboundary(tau = τ, wavenumber = 0.0)
    # Y = VIE.hhvolumegradG(tau = τ, wavenumber = 0.0)


    #Overwrite:
    #V = VIE.hhvolume(tau = τ, wavenumber = 0.0)

    # Exitation
    dirE = SVector(1.0, 0.0, 0.0)
    dirgradΦ = -dirE
    amp = 1.0

    # SphericalScattering exitation
    ex = UniformField(direction = dirE, amplitude = amp, embedding = Medium(ε1, μ1))

    # VIE exitation
    Φ_inc = VIE.linearpotential(direction=dirgradΦ, amplitude=amp)


    # Assembly
    b = assemble(Φ_inc, X)

    Z_I = assemble(I, X, X)

    Z_V = assemble(V, X, X)
    Z_B = assemble(B, strc(X), X)

    Z_Y = assemble(Y, X, X)

    Z_version1 = Z_I + Z_V + Z_B
    Z_version2 = Z_I + Z_Y

    # MoM solution
    u_version1 = Z_version1 \ Vector(b)
    u_version2 = Z_version2 \ Vector(b)

    # Observation points
    range_ = range(-1.0*r,stop=1.0*r,length=20)
    points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
    points_sp=[]
    for p in points
        norm(p)<0.97*r && push!(points_sp,p)
    end

    # SphericalScattering solution inside the dielectric sphere
    Φ = field(sp, ex, ScalarPotential(points_sp))

    # MoM solution inside the dielectric sphere
    Φ_MoM_version1 = BEAST.grideval(points_sp, u_version1, X, type=Float64)
    Φ_MoM_version2 = BEAST.grideval(points_sp, u_version2, X, type=Float64)



    err_Φ_version1 = norm(Φ - Φ_MoM_version1) / norm(Φ)
    err_Φ_version2 = norm(Φ - Φ_MoM_version2) / norm(Φ)

    @show err_Φ_version1
    @show err_Φ_version2

    @test err_Φ_version1 < 0.02
    @test err_Φ_version2 < 0.01


    #Symmetrie test BoundaryOperatorΓΩ vs. BoundaryOperatorΩΓ #########################
    struct TestOpBoundaryΓΩ{T,U,P} <: ImpedancePredictionVIE.BoundaryOperatorΓΩ
        gamma::T
        α::U
        tau::P
    end    
    struct TestOpBoundaryΩΓ{T,U,P} <: ImpedancePredictionVIE.BoundaryOperatorΩΓ
        gamma::T
        α::U
        tau::P
    end

    function testOpBoundaryΓΩ(; tau=nothing)

        return TestOpBoundaryΓΩ(0.0, 1.0, tau)
    end

    function testOpBoundaryΩΓ(; tau=nothing)

        return TestOpBoundaryΩΓ(0.0, 1.0, tau)
    end

    function BEAST.integrand(viop::TestOpBoundaryΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

        gx = @SVector[tvals[i].value for i in 1:3]
        fy = @SVector[bvals[i].value for i in 1:4]

        G = kerneldata.green

        Ty = kerneldata.tau

        α = viop.α

        return @SMatrix[α * gx[i] * G * fy[j] for i in 1:3, j in 1:4]
    end

    function BEAST.integrand(viop::TestOpBoundaryΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)

        gx = @SVector[tvals[i].value for i in 1:4]
        fy = @SVector[bvals[i].value for i in 1:3]

        G = kerneldata.green

        Ty = kerneldata.tau

        α = viop.α

        return @SMatrix[α * gx[i] * G * fy[j] for i in 1:4, j in 1:3]
    end

    T_ΓΩ = testOpBoundaryΓΩ(tau = τ)
    T_ΩΓ = testOpBoundaryΩΓ(tau = τ)

    #BEAST.defaultquadstrat(op::TestOpBoundaryΓΩ, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)
    #BEAST.defaultquadstrat(op::TestOpBoundaryΩΓ, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)

    M_ΓΩ = assemble(T_ΓΩ, strc(X), X)
    M_ΩΓ = assemble(T_ΩΓ, X, strc(X))
    #@show M_ΓΩ-transpose(M_ΩΓ)

    @test norm(M_ΓΩ - transpose(M_ΩΓ)) < 1e-14




    ####### BoundaryOperatorΓΩ vs. BoundaryOperator (oben minimal genauer - warum?) ###################################

    B_std = VIE.hhboundary(wavenumber = 0.0, tau = τ)
    B_ΓΩ = hhboundary2(tau = τ)

    #BEAST.defaultquadstrat(op::BEAST.VIEhhBoundary, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,6,6,6,6)
    #BEAST.defaultquadstrat(op::VIEhhBoundary2, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,6,6,6,6)

    Z_B_std = assemble(B_std, strc(X), X)
    Z_B_ΓΩ = assemble(B_ΓΩ, strc(X), X)

    @test norm(Z_B_std - Z_B_ΓΩ) < 9.0e-4           # Warum ist das eigentlich nicht exakt?
    norm(Z_B_std)
    norm(Z_B_ΓΩ)
    @test abs(maximum(Z_B_std - Z_B_ΓΩ)) < 3.0e-4
    @test abs(minimum(Z_B_std - Z_B_ΓΩ)) < 3.0e-4
    #displax(Z_B_std - Z_B_ΓΩ)



    ####### VolumeOperatorΩΩ vs. VolumeOperator (oben minimal genauer - warum?) ###################################

    V_std = VIE.hhvolume(wavenumber = 0.0, tau = τ)
    V_ΩΩ = hhvolume2(tau = τ)

    Z_V_std = assemble(V_std, X, X)
    Z_V_ΩΩ = assemble(V_ΩΩ, X, X)
    
    @test norm(Z_V_std - Z_V_ΩΩ) < 6.0e-3
    norm(Z_V_std)
    norm(Z_V_ΩΩ)
    @test abs(maximum(Z_V_std - Z_V_ΩΩ)) < 1.0e-3
    @test abs(minimum(Z_V_std - Z_V_ΩΩ)) < 1.0e-3
    #display(Z_V_std - Z_V_ΩΩ)


#end