#module IP

#using ..ImpedancePredictionVIE #damit alles exportierte hier funktioniert

using ..ImpedancePredictionVIE
using CompScienceMeshes
using BEAST
using LinearAlgebra
using StaticArrays

using Plots


abstract type geometric_body end
struct cuboid <: geometric_body
    L_x
    L_y
    L_z
end

struct cylinder <: geometric_body
    H
    R
end

struct meshdata
    body::geometric_body
    h::Float64
    Ω::Mesh{3, 4, Float64}
    Γ::Mesh{3, 3, Float64}
    Γ_c::Mesh{3, 3, Float64}
    Γ_c_t::Mesh{3, 3, Float64}
    Γ_c_b::Mesh{3, 3, Float64}
    Γ_nc::Mesh{3, 3, Float64}
    topnodes::Vector{Int64}
    bottomnodes::Vector{Int64}
    dirichletnodes::Vector{Int64}
    y_d::BEAST.LagrangeBasis{1, 0, Mesh{3, 3, Float64}, Float64, 3, SVector{3, Float64}}
    y::BEAST.LagrangeBasis{1, 0, Mesh{3, 3, Float64}, Float64, 3, SVector{3, Float64}}
    w::BEAST.LagrangeBasis{0, -1, Mesh{3, 3, Float64}, Float64, 1, SVector{3, Float64}}
    swg_faces::Vector{SVector{3, Int64}}
    X::BEAST.NDLCDBasis{Float64, Mesh{3, 4, Float64}, SVector{3, Float64}}
    ntrcX::BEAST.LagrangeBasis{0, -1, CompScienceMeshes.SubMesh{3, 3, Float64}, Float64, 1, SVector{3, Float64}}
    ntrc::Function
end

abstract type material end
struct constantmaterial <: material
    κ::Union{Float64, Nothing}
    ϵ::Union{Float64, Nothing}
end



struct solution
    meshdata::meshdata
    material::material
    #κ::Function
    κ0::Union{Float64, Nothing}
    #ϵ::Function
    ϵ0::Union{Float64, Nothing}
    ω::Union{Float64, Nothing}
    τ0::Union{Float64, ComplexF64}
    # τ
    # inv_τ
    # χ
    potential_top
    potential_bottom
    qs3D
    qs4D
    qs5D6D
    R
    v
    b
    S
    u
    u_Φ
    u_Jn
    u_J
end


function setup(; geoname::String = "cube.geo", meshname::String = "cube.msh",
    body::geometric_body = cuboid(1.0, 1.0, 1.0), h = 0.18)

    #geoname = "name.geo"
    geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"

    #meshname = "name.msh"
    meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"

    Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h, body = body)

    
    # LinearLag  on Γ_c (dirichlet)
    topnodes = realnodes(Γ_c_t) # wichtig denn später z.B. 10V
    bottomnodes = realnodes(Γ_c_b) # ""  z.B. 0V
    dirichletnodes = vcat(topnodes, bottomnodes)
    y_d = lagrangec0d1(Γ, dirichletnodes, Val{3})
    # dirichletnodes[i] <=> y_d.pos[i] ====> [10V 10V ..... 10V 0V 0V ..... 0V 0V] (example)

    # LinearLag on Γ_nc
    y = lagrangec0d1(Γ_nc, dirichlet = true) 

    # SWG on Ω (without faces on Γ_nc)
    swg_faces = swgfaces(Ω, Γ_nc, fast = true)
    X = nedelecd3d(Ω, Mesh(Ω.vertices, swg_faces))
    @assert length(X.pos) == length(swg_faces)
    ntrc = X -> BEAST.ntrace(X, Γ)
    ntrcX = ntrc(X)

    # PWC on Γ_c (vanish on Γ_nc)
    w = lagrangecxd0(Γ_c)

    @show numfunctions(y)
    @show numfunctions(y_d)
    @show numfunctions(w)
    @show numfunctions(X)

    return meshdata(body, h, Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc, topnodes,
    bottomnodes, dirichletnodes, y_d, y, w, swg_faces, X, ntrcX, ntrc)
end


function solve(;
    meshdata::meshdata, 
    κ::Union{Function, Nothing} = nothing, 
    κ0::Union{Float64, Nothing} = nothing, 
    ϵ::Union{Function, Nothing} = nothing, 
    ϵ0::Union{Float64, Nothing} = nothing,
    ω::Union{Float64, Nothing} = nothing, 
    potential_top::Float64, 
    potential_bottom::Float64,
    qs3D = BEAST.SingleNumQStrat(3), 
    qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4), #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
    qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4))


    # Material
    τ, inv_τ, τ0, χ = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    τ(p) < 1e-12 && error("Disable the following lines...")
    @assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) <1e-10
    @assert abs(1/inv_τ(p) -  τ(p)) <1e-10

    # Excitation
    v_top = ones(length(meshdata.topnodes)) * potential_top
    v_bottom = ones(length(meshdata.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1

    B11_Γ = IPVIE.B11_Γ(alpha = 1.0) #Verschwindet!!! <---- klären....
    #assemble(B11_Γ, w, y)
    B11_ΓΓ = IPVIE.B11_ΓΓ(alpha = 1.0, gammatype = Float64)
    #assemble(B11_ΓΓ, w, y)

    B12_ΓΓ = IPVIE.B12_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
    #assemble(B12_ΓΓ, w, w)

    B13_ΓΓ = IPVIE.B13_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
    #assemble(B13_ΓΓ, w, ntrc(X))
    B13_ΓΩ = IPVIE.B13_ΓΩ(alpha = -1.0, gammatype = Float64, chi = χ)
    #assemble(B13_ΓΩ, w, X)


    # Operators row 2

    B21_ΓΓ = IPVIE.B21_ΓΓ(beta = -1.0, gammatype = Float64) # -1.0 siehe BEAST & Steinbach 6.5
    #assemble(B21_ΓΓ, y, y)

    B22_Γ = IPVIE.B22_Γ(alpha = -1.0, invtau = inv_τ)
    #assemble(B22_Γ, y, w)
    B22_ΓΓ = IPVIE.B22_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
    #assemble(B22_ΓΓ, y, w)

    B23_ΓΓ = IPVIE.B23_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ) #VZ? sollte passen
    #assemble(B23_ΓΓ, y, ntrc(X))
    B23_ΓΩ = IPVIE.B23_ΓΩ(alpha = 1.0, gammatype = Float64, chi = χ)
    #assemble(B23_ΓΩ, y, X)


    # Operators row 3

    B31_ΓΓ = IPVIE.B31_ΓΓ(alpha = 1.0, gammatype = Float64)
    #assemble(B31_ΓΓ, ntrc(X), y)
    B31_ΩΓ = IPVIE.B31_ΩΓ(alpha = -1.0, gammatype = Float64)
    #assemble(B31_ΩΓ, X, y)

    B32_ΓΓ = IPVIE.B32_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
    #assemble(B32_ΓΓ, ntrc(X), w)
    B32_ΩΓ = IPVIE.B32_ΩΓ(alpha = -1.0, gammatype = Float64, invtau = inv_τ)
    #assemble(B32_ΩΓ, X, w)

    B33_Ω = IPVIE.B33_Ω(alpha = -1.0, invtau = inv_τ)
    #assemble(B33_Ω, X, X)
    B33_ΓΓ = IPVIE.B33_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
    #assemble(B33_ΓΓ, ntrc(X), ntrc(X))
    B33_ΓΩ = IPVIE.B33_ΓΩ(alpha = -1.0, gammatype = Float64, chi = χ)
    #assemble(B33_ΓΩ, ntrc(X), X)
    B33_ΩΓ = IPVIE.B33_ΩΓ(alpha = -1.0, gammatype = Float64, chi = χ)
    #assemble(B33_ΩΓ, X, ntrc(X))
    B33_ΩΩ = IPVIE.B33_ΩΩ(alpha = 1.0, gammatype = Float64, chi = χ)
    #assemble(B33_ΩΩ, X, X)
    

    # fns spaces
    y_d = meshdata.y_d
    y = meshdata.y
    w = meshdata.w
    X = meshdata.X
    ntrc = meshdata.ntrc 

    # LHS
    @hilbertspace i j k # row    -> test
    @hilbertspace l m n # col    -> trial
    lhs = @varform(
        B11_Γ[i,l] + B11_ΓΓ[i,l] +
        B12_ΓΓ[i,m] +
        B13_ΓΓ[i,ntrc(n)] + B13_ΓΩ[i,n] + 

        B21_ΓΓ[j,l] + 
        B22_Γ[j,m] + B22_ΓΓ[j,m] +
        B23_ΓΓ[j,ntrc(n)] + 
        B23_ΓΩ[j,n] +

        B31_ΓΓ[ntrc(k),l] + B31_ΩΓ[k,l] +
        B32_ΓΓ[ntrc(k),m] + B32_ΩΓ[k,m] +
        B33_Ω[k,n] + B33_ΓΓ[ntrc(k),ntrc(n)] + B33_ΓΩ[ntrc(k),n] + B33_ΩΓ[k,ntrc(n)] + B33_ΩΩ[k,n]
    )
    lhsd = @discretise lhs i∈w j∈y k∈X l∈y m∈w n∈X # w und y swapped for test!
    lhsd_test = lhsd.test_space_dict
    lhsd_trial = lhsd.trial_space_dict
    testSpace_lhs = BEAST._spacedict_to_directproductspace(lhsd_test)
    trialSpace_lhs = BEAST._spacedict_to_directproductspace(lhsd_trial)
    M = assemble(lhs, testSpace_lhs, trialSpace_lhs)
    S = Matrix(M)

    # RHS
    @hilbertspace o # col, only one blockcol
    rhs = @varform(
        -B11_Γ[i,o] -B11_ΓΓ[i,o] +

        -B21_ΓΓ[j,o] + 

        -B31_ΓΓ[ntrc(k),o] -B31_ΩΓ[k,o]
    )
    rhsd = @discretise rhs i∈w j∈y k∈X o∈y_d # w und y swapped for test!
    rhsd_test = rhsd.test_space_dict
    rhsd_trial = rhsd.trial_space_dict
    testSpace_rhs = BEAST._spacedict_to_directproductspace(rhsd_test)
    trialSpace_rhs = BEAST._spacedict_to_directproductspace(rhsd_trial)
    R = Matrix(assemble(rhs, testSpace_rhs, trialSpace_rhs))

    # S*u = R*v, solve for u
    b = R*v
    u = S \ b # it solver...
    #@assert norm(S*u - b) < 1e-8
    u_Φ = u[1:length(y)]
    u_Jn = u[length(y)+1:length(y)+length(w)]
    u_J = u[length(y)+length(w)+1:end]
    @assert length(u_Φ) == length(y.fns)
    @assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)


    @warn "Add control feature for the material functions!"
    # kappa epsilon tau invtau chi alternative hier....
    return solution(meshdata, material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, R, v, b, S, u, u_Φ, u_Jn, u_J)
end





#end















# struct fnsdata


# end


# mutable struct 