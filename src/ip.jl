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
function (mat::constantmaterial)()
    κ = x -> mat.κ
    ϵ = x -> mat.ϵ
    mat.κ === nothing && (κ = nothing)
    mat.ϵ === nothing && (ϵ = nothing)
    return κ, ϵ
end

struct constant_zsplit <: material
    κ_p::Union{Float64, Nothing}
    ϵ_p::Union{Float64, Nothing}
    z0::Float64
    κ_m::Union{Float64, Nothing}
    ϵ_m::Union{Float64, Nothing}
end
function (mat::constant_zsplit)()

    function gen()
        function kappa(x) 
            x[3] >= mat.z0 && (return mat.κ_p)
            return mat.κ_m
        end
        function epsilon(x)
            x[3] >= mat.z0 && (return mat.ϵ_p)
            return mat.ϵ_m
        end
        (mat.κ_p === nothing || mat.κ_m === nothing) && (kappa = nothing)
        (mat.ϵ_p === nothing || mat.ϵ_m === nothing) && (epsilon = nothing)
        return kappa, epsilon
    end
    κ, ϵ = gen()

    return κ, ϵ
end

struct constant_xsplit <: material
    κ_p::Union{Float64, Nothing}
    ϵ_p::Union{Float64, Nothing}
    x0::Float64
    κ_m::Union{Float64, Nothing}
    ϵ_m::Union{Float64, Nothing}
end
function (mat::constant_xsplit)()

    function gen()
        function kappa(x) 
            x[1] >= mat.x0 && (return mat.κ_p)
            return mat.κ_m
        end
        function epsilon(x)
            x[1] >= mat.x0 && (return mat.ϵ_p)
            return mat.ϵ_m
        end
        (mat.κ_p === nothing || mat.κ_m === nothing) && (kappa = nothing)
        (mat.ϵ_p === nothing || mat.ϵ_m === nothing) && (epsilon = nothing)
        return kappa, epsilon
    end
    κ, ϵ = gen()

    return κ, ϵ
end

struct general_material <: material # => reconstruction is a problem!!!
    κ::Union{Function, Nothing}
    ϵ::Union{Function, Nothing}
    # savesomepointswithmat::Union{....Vektoren... , Nothing} => init mit nothing
end
function (mat::general_material)
    return mat.κ, mat.ϵ
end




struct solution
    material::material # κ::Function, ϵ::Function -> τ, inv_τ, χ über funktion bestimmbar...
    κ0::Union{Float64, Nothing}
    ϵ0::Union{Float64, Nothing}
    ω::Union{Float64, Nothing}
    τ0::Union{Float64, ComplexF64}
    
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
    md::meshdata,
    material::material,
    κ0::Union{Float64, Nothing} = nothing, 
    ϵ0::Union{Float64, Nothing} = nothing, # VORSICHT zsh zum echten epsilon0 nicht unbed. gegeben!
    ω::Union{Float64, Nothing} = nothing, 
    potential_top::Float64, 
    potential_bottom::Float64,
    qs3D = BEAST.SingleNumQStrat(3), 
    qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4), #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
    qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4))

    # if typeof(material) == general_material
        # save some points, and ...
        # material = general_material(...,...,This time not nothing here=>some test points with mat values)
    # end

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    @show τ(p)

    τ(p) < 1e-12 && error("Disable the following lines...")
    @assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) < 1e-10
    @assert abs(1/inv_τ(p) -  τ(p)) < 1e-10

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
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
    y_d = md.y_d
    y = md.y
    w = md.w
    X = md.X
    ntrc = md.ntrc 

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


    # qs3Dhigh = qs3D
    # qs4Dhigh = qs4D
    # qs5D6Dhigh = qs5D6D

    # qs3Dhigh = BEAST.SingleNumQStrat(6)
    # qs4Dhigh = BEAST.DoubleNumWiltonSauterQStrat(6,6,6,6,6,6,6,6)
    # qs5D6Dhigh = BEAST.SauterSchwab3DQStrat(6,6,6,6,6,6)
    
    # O11 = -B11_Γ
    # R11 = O11.coeffs[1]*assemble(O11.ops[1], w, y_d; quadstrat = qs3Dhigh)

    # O11_ = -B11_ΓΓ
    # R11_= O11_.coeffs[2]*assemble(O11_.ops[2], w, y_d; quadstrat = qs4Dhigh) +
    #         O11_.coeffs[1]*assemble(O11_.ops[1], w, y_d; quadstrat = qs3Dhigh)

    # O21 = -B21_ΓΓ
    # R21 = O21.coeffs[1]*assemble(O21.ops[1], y, y_d; quadstrat = qs4Dhigh)

    # O31 = -B31_ΓΓ
    # R31 = O31.coeffs[2]*assemble(O31.ops[2], ntrc(X), y_d; quadstrat = qs4Dhigh) + 
    #         O31.coeffs[1]*assemble(O31.ops[1], ntrc(X), y_d; quadstrat = qs3Dhigh)

    # O31_ = -B31_ΩΓ     
    # R31_= O31_.coeffs[1]*assemble(O31_.ops[1], X, y_d; quadstrat = qs5D6Dhigh)

    # R = vcat(R11+R11_, R21, R31+R31_)


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


    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, R, v, b, S, u, u_Φ, u_Jn, u_J)
end








###### constantmaterial - analytic solution ####################

# function solution_J_ana(m::IP.meshdata, s::IP.solution, points, J_MoM; body::cuboid, mat::constantmaterial)
    
#     Jz_ana = -mat.κ*(s.potential_top-s.potential_bottom)/body.L_z

#     J_ana = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM))

#     return J_ana
# end

# function solution_I_ana(m::IP.meshdata, s::IP.solution; body::cuboid, mat::constantmaterial)

#     I_ana = mat.κ*body.L_x*body.L_y*(1/body.L_z)*(s.potential_top-s.potential_bottom)

#     return I_ana
# end

# function solution_Φ_ana(m::IP.meshdata, s::IP.solution; body::cuboid, mat::constantmaterial)
#     u_Φ = s.u_Φ

#     function CapacitorPot(x, v_top, v_bottom) # hom. medium
#         L_z = body.L_z
#         z = x[3]

#         return ((v_top-v_bottom)/L_z) * (z + L_z/2) + v_bottom
#     end
#     u_Φ_ana = Vector{Float64}(undef,length(u_Φ))
#     y = m.y
#     for (i, pos) in enumerate(y.pos)
#         u_Φ_ana[i] = CapacitorPot(pos, s.potential_top, s.potential_bottom)
#     end

#     return u_Φ_ana
# end


###### constant_zsplit Material - analytic solution ##################

function solution_J_ana(m::IP.meshdata, s::IP.solution, points, J_MoM; body::IP.cuboid, mat::IP.constant_zsplit)
    
    A = body.L_x * body.L_y
    R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    U = s.potential_top - s.potential_bottom
    I = U/R
    Jz_ana = -I/A 

    J_ana = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM))

    return J_ana
end

function solution_I_ana(m::IP.meshdata, s::IP.solution; body::IP.cuboid, mat::IP.constant_zsplit)

    A = body.L_x * body.L_y
    R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    U = s.potential_top - s.potential_bottom
    I = U/R

    I_ana = I

    return I_ana
end

function solution_Φ_ana(m::IP.meshdata, s::IP.solution; body::IP.cuboid, mat::IP.constant_zsplit)
    u_Φ = s.u_Φ

    function linPot_z(x, z1, z2, Φ1, Φ2) # linear in section [z1,z2] -> [Φ1, Φ2]
        z = x[3] # obervation plane
        L_z = body.L_z

        return ((Φ2 - Φ1)/(z2 - z1))*(z - z1) + Φ1
    end

    u_Φ_ana = Vector{Float64}(undef,length(u_Φ))

    A = body.L_x * body.L_y
    R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    U = s.potential_top - s.potential_bottom
    I = U/R
    R_p = (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    R_m = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A
    U_p = U * R_p/R
    U_m = U * R_m/R
    
    
    z0 = mat.z0
    y = m.y
    for (i, pos) in enumerate(y.pos)
        if pos[3] >= z0 # R_p region

            z1 = z0
            Φ1 = s.potential_top - U_p
            z2 = body.L_z/2
            Φ2 = s.potential_top

        elseif pos[3] < z0 # R_m region

            z1 = -body.L_z/2
            Φ1 = s.potential_bottom
            z2 = z0
            Φ2 = s.potential_bottom + U_m

        else
            error("")
        end

        u_Φ_ana[i] = linPot_z(pos, z1, z2, Φ1, Φ2)
    end

    return u_Φ_ana
end


###### constant_xsplit Material - analytic solution ##################

# function solution_J_ana(m::IP.meshdata, s::IP.solution, points, J_MoM; body::IP.cuboid, mat::IP.constant_xsplit)
    
#     x0 = mat.x0
#     A_p = body.L_y * (body.L_x/2 - x0)
#     A_m = body.L_y * (body.L_x/2 + x0)

#     R_p = (1/mat.κ_p)*body.L_z/A_p
#     R_m = (1/mat.κ_m)*body.L_z/A_m

#     U = s.potential_top - s.potential_bottom
#     I_p = U/R_p
#     I_m = U/R_m

#     Jz_p_ana = -I_p/A_p
#     Jz_m_ana = -I_m/A_m 

#     @assert length(points) == length(J_MoM)

#     J_ana = Vector{SVector{3, Float64}}(undef, length(J_MoM)) #(0.0, 0.0, Jz_ana), length(J_MoM))
#     for (i,pos) in enumerate(points)
#         x = pos[1]
#         if x >= x0
#             J_ana[i] = SVector{3, Float64}(0.0, 0.0, Jz_p_ana)
#         elseif x < x0
#             J_ana[i] = SVector{3, Float64}(0.0, 0.0, Jz_m_ana)
#         else
#             error("")
#         end
#     end

#     return J_ana
# end

# function solution_I_ana(m::IP.meshdata, s::IP.solution; body::IP.cuboid, mat::IP.constant_xsplit)

#     x0 = mat.x0
#     A_p = body.L_y * (body.L_x/2 - x0)
#     A_m = body.L_y * (body.L_x/2 + x0)

#     R_p = (1/mat.κ_p)*body.L_z/A_p
#     R_m = (1/mat.κ_m)*body.L_z/A_m

#     U = s.potential_top - s.potential_bottom
#     I_p = U/R_p
#     I_m = U/R_m
#     I = I_p + I_p

#     I_ana = I

#     return I_ana
# end

# function solution_Φ_ana(m::IP.meshdata, s::IP.solution; body::IP.cuboid, mat::IP.constant_xsplit)
#     u_Φ = s.u_Φ

#     function linPot_z(x, z1, z2, Φ1, Φ2) # linear in section [z1,z2] -> [Φ1, Φ2]
#         z = x[3] # obervation plane
#         L_z = body.L_z

#         return ((Φ2 - Φ1)/(z2 - z1))*(z - z1) + Φ1
#     end

#     u_Φ_ana = Vector{Float64}(undef,length(u_Φ))
#     y = m.y

#     z1 = -body.L_z/2
#     Φ1 = s.potential_bottom
#     z2 = body.L_z/2
#     Φ2 = s.potential_top

#     for (i, pos) in enumerate(y.pos)
#         u_Φ_ana[i] = linPot_z(pos, z1, z2, Φ1, Φ2)
#     end

#     return u_Φ_ana
# end






# analogen wäre Dn =? Flächenladunsdichte
function getcurrent(m::IP.meshdata, s::IP.solution) # ging das jetzt noch einfacher oder nicht? assemblydata???
    
    w = m.w
    Γ_c_t = m.Γ_c_t
    Γ_c_b = m.Γ_c_b
    u_Jn = s.u_Jn
    
    @assert length(u_Jn) == length(w.fns)

    numtopcharts = length(Γ_c_t.faces)
    top_ccs = Vector{SVector{3, Float64}}(undef, numtopcharts)
    top_charts = Vector{CompScienceMeshes.Simplex{3, 2, 1, 3, Float64}}(undef, numtopcharts)
    for i in 1:length(top_ccs)
        chart = CompScienceMeshes.chart(Γ_c_t, i)
        center = CompScienceMeshes.center(chart)
        top_ccs[i] = cartesian(center)
        top_charts[i] = chart
    end

    numbottomcharts = length(Γ_c_b.faces)
    bottom_ccs = Vector{SVector{3, Float64}}(undef, numbottomcharts)
    bottom_charts = Vector{CompScienceMeshes.Simplex{3, 2, 1, 3, Float64}}(undef, numtopcharts)
    for i in 1:length(bottom_ccs)
        chart = CompScienceMeshes.chart(Γ_c_b, i)
        center = CompScienceMeshes.center(chart)
        bottom_ccs[i] = cartesian(center)
        bottom_charts[i] = chart
    end

    chart_tree_top = BEAST.octree(top_charts)
    chart_tree_bottom = BEAST.octree(bottom_charts)

    I_top = 0.0
    I_bottom = 0.0

    cnt_top = 0
    cnt_bottom = 0

    for (j, pos) in enumerate(w.pos)

        i = CompScienceMeshes.findchart(top_charts, chart_tree_top, pos)      # ja...hier wäre das nicht nötig assemblydata?
        k = CompScienceMeshes.findchart(bottom_charts, chart_tree_bottom, pos)


        if i !== nothing
            A = top_charts[i].volume 
            I_top += -A * u_Jn[j] * w.fns[j][1].coeff # "-" because dÂ of ∫∫J_vec*dÂ in opposite dir
            cnt_top += 1
        elseif k !== nothing
            A = bottom_charts[k].volume 
            I_bottom += A * u_Jn[j] * w.fns[j][1].coeff 
            cnt_bottom += 1
        else
            error("Neither top nor bottom")
        end

    end


    @assert cnt_top == numtopcharts
    @assert cnt_bottom == numbottomcharts

    return I_top, I_bottom
end







# function get_Dz_ref(zbottom, ztop, Φbottom, Φtop, ϵ_vec, d_vec) # ein skalarer Wert zurück, ϵ_vec and d_vec bottom to top!
#     n = length(ϵ_vec)
#     @assert n == length(d_vec)
#     @assert ztop > zbottom

#     d = abs(ztop - zbottom)

#     @assert abs(sum(d_vec)-d) < 1e-12 

#     z_vec = Vector{n, Float64}()

#     for (i,d_i) in enumerate(d_vec)

#     end

#     z_vec[1] = zbottom + d[1]
#     for i = 2:n-1
#         d_i = d[i]
        

#     end

# end




#end















# struct fnsdata


# end


# mutable struct 