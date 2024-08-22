using ..ImpedancePredictionVIE
using CompScienceMeshes
using BEAST
using LinearAlgebra
using StaticArrays
using SparseArrays

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
    Y_d::BEAST.LagrangeBasis{1, 0, Mesh{3, 4, Float64}, Float64, 4, SVector{3, Float64}}
    Y::BEAST.LagrangeBasis{1, 0, Mesh{3, 4, Float64}, Float64, 4, SVector{3, Float64}}
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
        (mat.κ_p === nothing && mat.κ_m === nothing) && (kappa = nothing)  # doch mal weglassen
        (mat.ϵ_p === nothing && mat.ϵ_m === nothing) && (epsilon = nothing)

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

struct pwlinx <: material
        κ::Union{Vector{Vector{Float64}}, Nothing}
        ϵ::Union{Vector{Vector{Float64}}, Nothing}
        xpos::Vector{Float64} # ! mit anfang und ende d.h -L.x/2 und L.x/2
end
function (mat::pwlinx)()
    κ_vec = mat.κ
    ϵ_vec = mat.ϵ
    xpos = mat.xpos
    l = length(xpos)
    κ_vec !== nothing && @assert l == length(κ_vec) + 1
    ϵ_vec !== nothing && @assert l == length(ϵ_vec) + 1

    function gen(l_)
        if κ_vec !== nothing
            function kappa(x)
                x1 = x[1]
                for i in 1:l_ # faster
                    x1 < xpos[i] && continue
                    x1 > xpos[i+1] && continue
                    return ((κ_vec[i][2] - κ_vec[i][1])/(xpos[i+1] - xpos[i]))*(x1 - xpos[i]) + κ_vec[i][1]
                end
                error("")
            end
        else
            kappa = nothing 
        end
        if ϵ_vec !== nothing
            function epsilon(x)
                x1 = x[1]
                for i in 1:l_ # faster
                    x1 < xpos[i] && continue
                    x1 > xpos[i+1] && continue
                    return ((ϵ_vec[i][2] - ϵ_vec[i][1])/(xpos[i+1] - xpos[i]))*(x1 - xpos[i]) + ϵ_vec[i][1]
                end
                error("")
            end
        else
            epsilon = nothing 
        end
        return kappa, epsilon
    end
    κ, ϵ = gen(l)

    return κ, ϵ
end


struct general_material <: material # => reconstruction is a problem!!!
    κ::Union{Function, Nothing}
    ϵ::Union{Function, Nothing}
    # savesomepointswithmat::Union{....Vektoren... , Nothing} => init mit nothing
end
function (mat::general_material)()
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
    v
    b
    u
    u_Φ
    u_Jn
    u_J
end





function setup(; geoname::Union{String, Nothing} = nothing, meshname::Union{String, Nothing} = nothing,
    body::geometric_body = cuboid(1.0, 1.0, 1.0), h = 0.18)

    if geoname === nothing && meshname === nothing
        error("Missing geoname/meshname")
    elseif geoname === nothing && meshname !== nothing
        meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"
        Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = msh2mesh(meshpath)
    elseif meshname === nothing && geoname !== nothing
        geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"
        meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/($geoname)_"
        Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h, body = body)
    else  
        geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"
        meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"
        Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h, body = body)
    end

    # LinearLag  on Γ_c (dirichlet)
    topnodes = realnodes(Γ_c_t) # wichtig denn später z.B. 10V
    bottomnodes = realnodes(Γ_c_b) # ""  z.B. 0V
    dirichletnodes = vcat(topnodes, bottomnodes)
    y_d = lagrangec0d1(Γ, dirichletnodes, Val{3})
    # dirichletnodes[i] <=> y_d.pos[i] ====> [10V 10V ..... 10V 0V 0V ..... 0V 0V] (example)

    # LinearLag on Γ_nc
    #y = lagrangec0d1(Γ_nc, dirichlet = true) <----problem beim plotten...
    bndnodes = realnodes(Γ)
    set1 = Set(bndnodes) # Konvertiere Arrays in Sets für schnelles Nachschlagen
    set2 = Set(dirichletnodes)
    result_set = setdiff(set1, set2) # Berechne die Differenzmenge
    dofnodes = collect(result_set) # Ausgabe des Ergebnisses
    y = lagrangec0d1(Γ, dofnodes, Val{3})


    # SWG on Ω (without faces on Γ_nc)
    #swg_faces = swgfaces(Ω, Γ_nc, fast = true)
    println("Time for calculation of swg_faces:")
    @time swg_faces = swgfaces_set_approach(Ω, Γ_nc)
    println()
    X = nedelecd3d(Ω, Mesh(Ω.vertices, swg_faces))
    @assert length(X.pos) == length(swg_faces)
    ntrc = X -> BEAST.ntrace(X, Γ)
    ntrcX = ntrc(X)

    # PWC on Γ_c (vanish on Γ_nc)
    #w = lagrangecxd0(Γ_c)
    #@warn "If w is involved problems can occur now - see setup()"
    newfns = Vector{Vector{BEAST.Shape{Float64}}}()
    newpos = Vector{SVector{3, Float64}}()
    for (i,shs) in enumerate(ntrcX.fns)
        newshs = Vector{BEAST.Shape{Float64}}()
        shs == [] && continue
        @assert length(shs) == 1
        push!(newfns, shs)
        push!(newpos, ntrcX.pos[i])
    end
    w = BEAST.LagrangeBasis{0,-1,1}(ntrcX.geo.supermesh, newfns, newpos)

    # FEM only
    # nondirichletnodes = Vector{Int64}()
    # println("Time for FEM node preparation")
    # @time for k in 1:length(Ω.vertices)
    #     push_ = true
    #     for n in dirichletnodes 
    #         if k == n 
    #             push_ = false
    #             break
    #         end
    #     end
    #     push_ && push!(nondirichletnodes, k)
    # end
    # @assert length(nondirichletnodes) + length(dirichletnodes) == length(Ω.vertices)
    # Y_d = lagrangec0d1(Ω, dirichletnodes, Val{4})
    # Y = lagrangec0d1(Ω, nondirichletnodes, Val{4})
    allnodes = collect(1:length(Ω.vertices))
    set1 = Set(allnodes) # Konvertiere Arrays in Sets für schnelles Nachschlagen
    set2 # = Set(dirichletnodes)
    result_set = setdiff(set1, set2) # Berechne die Differenzmenge
    nondirichletnodes = collect(result_set) # Ausgabe des Ergebnisses
    @assert length(nondirichletnodes) + length(dirichletnodes) == length(Ω.vertices)
    Y_d = lagrangec0d1(Ω, dirichletnodes, Val{4})
    Y = lagrangec0d1(Ω, nondirichletnodes, Val{4})


    @show numfunctions(y)
    @show numfunctions(y_d)
    @show numfunctions(w)
    @show numfunctions(X)
    @show numfunctions(Y)
    @show numfunctions(Y_d)

    return meshdata(body, h, Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc, topnodes,
    bottomnodes, dirichletnodes, y_d, y, w, swg_faces, X, ntrcX, ntrc, Y_d, Y)
end





function solve(; # low contrast formulation
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

    println("3×3 block matrix - classical formulation")

    # if typeof(material) == general_material
        # save some points, and ...
        # material = general_material(...,...,This time not nothing here=>some test points with mat values)
    # end

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    @show τ(p)
    @show inv_τ(p)
    @show τ0
    @show χ(p)
    @show T

    # τ(p) < 1e-12 && error("Disable the following lines...")
    # @assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) < 1e-10
    # @assert abs(1/inv_τ(p) -  τ(p)) < 1e-10

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1

    B11_Γ = IPVIE.B11_Γ(alpha = 1.0) #Verschwindet!!! <---- klären....
    #assemble(B11_Γ, w, y)
    B11_ΓΓ = IPVIE.B11_ΓΓ(alpha = 1.0, gammatype = Float64) # Float64 wenn kein Material
    #assemble(B11_ΓΓ, w, y)

    B12_ΓΓ = IPVIE.B12_ΓΓ(alpha = 1.0, gammatype = T, invtau = inv_τ)
    #assemble(B12_ΓΓ, w, w)

    B13_ΓΓ = IPVIE.B13_ΓΓ(alpha = 1.0, gammatype = T, chi = χ)
    #assemble(B13_ΓΓ, w, ntrc(X))
    B13_ΓΩ = IPVIE.B13_ΓΩ(alpha = -1.0, gammatype = T, chi = χ)
    #assemble(B13_ΓΩ, w, X)


    # Operators row 2

    B21_ΓΓ = IPVIE.B21_ΓΓ(beta = -1.0, gammatype = Float64) # -1.0 siehe BEAST & Steinbach 6.5, Ja, Float64 statt T
    #assemble(B21_ΓΓ, y, y)

    B22_Γ = IPVIE.B22_Γ(alpha = -1.0, invtau = inv_τ) # !!! ggf hier T(-1.0) statt -1.0 ...fraglich wieso nur unten Fehler
    #assemble(B22_Γ, y, w)
    B22_ΓΓ = IPVIE.B22_ΓΓ(alpha = 1.0, gammatype = T, invtau = inv_τ)
    #assemble(B22_ΓΓ, y, w)

    B23_ΓΓ = IPVIE.B23_ΓΓ(alpha = 1.0, gammatype = T, chi = χ) #VZ? sollte passen
    #assemble(B23_ΓΓ, y, ntrc(X))
    B23_ΓΩ = IPVIE.B23_ΓΩ(alpha = 1.0, gammatype = T, chi = χ)
    #assemble(B23_ΓΩ, y, X)


    # Operators row 3

    B31_ΓΓ = IPVIE.B31_ΓΓ(alpha = 1.0, gammatype = Float64)
    #assemble(B31_ΓΓ, ntrc(X), y)
    B31_ΩΓ = IPVIE.B31_ΩΓ(alpha = -1.0, gammatype = Float64)
    #assemble(B31_ΩΓ, X, y)

    B32_ΓΓ = IPVIE.B32_ΓΓ(alpha = 1.0, gammatype = T, invtau = inv_τ)
    #assemble(B32_ΓΓ, ntrc(X), w)
    B32_ΩΓ = IPVIE.B32_ΩΓ(alpha = -1.0, gammatype = T, invtau = inv_τ)
    #assemble(B32_ΩΓ, X, w)

    B33_Ω = IPVIE.B33_Ω(alpha = T(-1.0), invtau = inv_τ) # def scalartype(MaterialIdentity) = alpha ...
    #assemble(B33_Ω, X, X)
    B33_ΓΓ = IPVIE.B33_ΓΓ(alpha = 1.0, gammatype = T, chi = χ)
    #assemble(B33_ΓΓ, ntrc(X), ntrc(X))
    B33_ΓΩ = IPVIE.B33_ΓΩ(alpha = -1.0, gammatype = T, chi = χ)
    #assemble(B33_ΓΩ, ntrc(X), X)
    B33_ΩΓ = IPVIE.B33_ΩΓ(alpha = -1.0, gammatype = T, chi = χ)
    #assemble(B33_ΩΓ, X, ntrc(X))
    B33_ΩΩ = IPVIE.B33_ΩΩ(alpha = 1.0, gammatype = T, chi = χ)
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


    # qs3Dhigh = BEAST.SingleNumQStrat(10)
    # qs4Dhigh = BEAST.DoubleNumWiltonSauterQStrat(10,10,10,10,10,10,10,10)
    # qs5D6Dhigh = BEAST.SauterSchwab3DQStrat(10,10,10,1,1,1)
    
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

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, u_Jn, u_J), S, R
end


function solve0(; # low contrast formulation, same as solve but operators are easier accessible
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

    println("3×3 block matrix - classical formulation")

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)

    # Operators row 1
    B11_Γ = Identity() #Verschwindet!!! <---- klären....
    B11_ΓΓ = Helmholtz3D.doublelayer(gamma = T(0.0), alpha = 1.0)# + (-1/2)*Identity()

    B12_ΓΓ =  IP.MaterialSL(T(0.0), 1.0, inv_τ)

    B13_ΓΓ = IP.MaterialSL(T(0.0), 1.0, χ)
    B13_ΓΩ = IP.gradG_ΓΩ(T(0.0), -1.0, χ)


    # Operators row 2
    B21_ΓΓ = Helmholtz3D.hypersingular(gamma = T(0.0), beta = -1.0) # das versteckt alpha std. Null für gamma =0.0!!!  # -1.0 siehe BEAST & Steinbach 6.5, Ja, Float64 statt T

    B22_Γ = IP.MatId(T(-1.0), inv_τ)
    B22_ΓΓ = IP.MaterialADL(T(0.0), 1.0, inv_τ)

    B23_ΓΓ = IP.MaterialADL(T(0.0), 1.0, χ)
    B23_ΓΩ = IP.ncgrad_gradGc_ΓΩ(T(0.0), 1.0, χ)


    # Operators row 3
    B31_ΓΓ = Helmholtz3D.doublelayer(gamma = T(0.0), alpha = 1.0)# + (-1/2)*Identity()
    B31_ΩΓ = IP.div_ngradG_ΩΓ(T(0.0), -1.0, x->1.0)

    B32_ΓΓ = IP.MaterialSL(T(0.0), 1.0, inv_τ)
    B32_ΩΓ = IP.div_G_ΩΓ(T(0.0), -1.0, inv_τ)

    B33_Ω =  IP.MatId(T(-1.0), inv_τ) # different_tau scalartype(MaterialIdentity) = alpha ...
    B33_ΓΓ = IP.MaterialSL(T(0.0), 1.0, χ)
    B33_ΓΩ = IP.n_gradG_ΓΩ(T(0.0), -1.0, χ)
    B33_ΩΓ = IP.div_G_ΩΓ(T(0.0), -1.0, χ)
    B33_ΩΩ = IP.div_gradG_ΩΩ(T(0.0), 1.0, χ)

    # fns spaces
    y_d = md.y_d
    y = md.y
    w = md.w
    X = md.X
    ntrcX = md.ntrcX 
    #wdual = duallagrangec0d1(md.Γ, md.Γ_c) # DoubleNumSauterQstrat ist ohne wilton!
    #ydual = duallagrangecxd0(md.Γ, md.Γ_nc)
    t_eq1 = w
    t_eq2 = y

    # assemble
    B11 = assemble(B11_Γ, t_eq1, y) + assemble(B11_ΓΓ, t_eq1, y) -(1/2)*assemble(Identity(), t_eq1, y)
    B12 = assemble(B12_ΓΓ, t_eq1, w)
    B13 = assemble(B13_ΓΓ, t_eq1, ntrcX) + assemble(B13_ΓΩ, t_eq1, X) 

    B21 = assemble(B21_ΓΓ, t_eq2, y)
    B22 = assemble(B22_Γ, t_eq2, w) + assemble(B22_ΓΓ, t_eq2, w)
    B23 = assemble(B23_ΓΓ, t_eq2, ntrcX) + assemble(B23_ΓΩ, t_eq2, X)
    
    B31 = assemble(B31_ΓΓ, ntrcX, y) -(1/2)*assemble(Identity(), ntrcX, y) + assemble(B31_ΩΓ, X, y) 
    B32 = assemble(B32_ΓΓ, ntrcX, w) + assemble(B32_ΩΓ, X, w)
    B33 = assemble(B33_Ω, X, X) +
        assemble(B33_ΓΓ, ntrcX, ntrcX) +
        assemble(B33_ΓΩ, ntrcX, X) +
        assemble(B33_ΩΓ, X, ntrcX) + 
        assemble(B33_ΩΩ, X, X)

    R11 = -assemble(B11_Γ, t_eq1, y_d) -assemble(B11_ΓΓ, t_eq1, y_d) +(1/2)*assemble(Identity(), t_eq1, y_d)
    R21 = -assemble(B21_ΓΓ, t_eq2, y_d)
    R31 = -assemble(B31_ΓΓ, ntrcX, y_d) +(1/2)*assemble(Identity(), ntrcX, y_d) -assemble(B31_ΩΓ, X, y_d)

    #error("STOP HERE!")

    #ROW1 = hcat(B21,B22,B23) #swap 1&2
    #ROW2 = hcat(B11,B12,B13) # "-"
    #ROW3 = hcat(B31,B32,B33)
    #S = vcat(ROW1,ROW2,ROW3)
    R = vcat(R21,R11,R31) # "-"
    b = R*v


    rownrmB21 = norm(B21) #[norm(B21[i, :]) for i in 1:size(B21, 1)] #
    rownrmB12 = norm(B12) #[norm(B12[i, :]) for i in 1:size(B12, 1)] #
    rownrmB33 = norm(B33) #[norm(B33[i, :]) for i in 1:size(B33, 1)] #

    α_T = 1.0 #./norm(b[1:length(w.fns)]) 
    α_B = 1.0 #./(α_T.*rownrmB21)

    β_T = 1.0 #./norm(b[length(w.fns)+1:length(w.fns)+length(y.fns)]) 
    β_B = 1.0 #./(β_T.*rownrmB12)

    γ_T = 1.0 #./norm(b[length(w.fns)+length(y.fns)+1:end]) 
    γ_B = 1.0 #./(γ_T.*rownrmB33)

    # @show α_T
    # @show α_B
    # @show β_T
    # @show β_B
    # @show γ_T
    # @show γ_B

    ROW1 = α_T.*hcat(α_B.*B21, β_B.*B22, γ_B.*B23) #swap 1&2
    ROW2 = β_T.*hcat(α_B.*B11, β_B.*B12, γ_B.*B13) # "-"
    ROW3 = γ_T.*hcat(α_B.*B31, β_B.*B32, γ_B.*B33)
    S = vcat(ROW1,ROW2,ROW3)
    R = vcat(α_T.*R21, β_T.*R11, γ_T.*R31) # new
    b = R*v

    # @show opnorm(B11)
    # @show opnorm(B12)
    # @show opnorm(B13)

    # @show opnorm(B21)
    # @show opnorm(B22)
    # @show opnorm(B23)

    # @show opnorm(B31)
    # @show opnorm(B32)
    # @show opnorm(B33)



    # @show cond(B11)
    # @show cond(B12)
    # @show cond(B13)

    # @show cond(B21)
    # @show cond(B22)
    # @show cond(B23)

    # @show cond(B31)
    # @show cond(B32)
    # @show cond(B33)


    # S*u = b, solve for u
    u = S \ b # it solver...
    @assert norm(S*u - b) < 1e-8
    u_Φ = u[1:length(y)].*α_B
    u_Jn = u[length(y)+1:length(y)+length(w)].*β_B
    u_J = u[length(y)+length(w)+1:end].*γ_B
    #u = [u_Φ, u_Jn, u_J]
    @assert length(u_Φ) == length(y.fns)
    @assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, u_Jn, u_J), S, R
end


function solvefembem(;
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

    println("2×2 block matrix - fembem formulation")

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    @show τ(p)
    @show inv_τ(p)
    @show τ0
    @show χ(p)
    @show T

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1
    B11_Ω = IP.grad_grad_Ω(T(1.0), τ)

    B12_Ω = IP.grad__Ω(T(1.0), x->1.0)

    # Operators row 2
    B21_ΓΓ = Helmholtz3D.doublelayer(gamma = T(0.0), alpha = 1.0) # + (-1/2)*Identity()
    B21_ΩΓ = IP.div_ngradG_ΩΓ(T(0.0), -1.0, x->1.0)

    B22_ΓΓ_ = IP.MaterialSL(T(0.0), 1.0, inv_τ)
    B22_ΩΓ_ = IP.div_G_ΩΓ(T(0.0), -1.0, inv_τ)
    B22_Ω =  IP.MatId(T(-1.0), inv_τ) # different_tau scalartype(MaterialIdentity) = alpha ...
    B22_ΓΓ = IP.MaterialSL(T(0.0), 1.0, χ)
    B22_ΓΩ = IP.n_gradG_ΓΩ(T(0.0), -1.0, χ)
    B22_ΩΓ = IP.div_G_ΩΓ(T(0.0), -1.0, χ)
    B22_ΩΩ = IP.div_gradG_ΩΩ(T(0.0), 1.0, χ)

    # fns spaces
    y_d = md.y_d
    y = md.y
    w = md.w
    X = md.X
    ntrcX = md.ntrcX 
    X = md.X
    Y = md.Y
    Y_d = md.Y_d
    strcY = strace(Y, md.Γ_nc)
    strcY_d = strace(Y_d, md.Γ_c)

    # assemble
    B11 = assemble(B11_Ω, Y, Y)
    B12 = assemble(B12_Ω, Y, X)
    
    B21 = assemble(B21_ΓΓ, ntrcX, strcY) -(1/2)*assemble(Identity(), ntrcX, strcY) + assemble(B21_ΩΓ, X, strcY)
    B22 = assemble(B22_ΓΓ_, ntrcX, ntrcX) +
        assemble(B22_ΩΓ_, X, ntrcX) +
        assemble(B22_Ω, X, X) +
        assemble(B22_ΓΓ, ntrcX, ntrcX) +
        assemble(B22_ΓΩ, ntrcX, X) +
        assemble(B22_ΩΓ, X, ntrcX) + 
        assemble(B22_ΩΩ, X, X)
        

    R11 = -assemble(B11_Ω, Y, Y_d)
    R21 = -assemble(B21_ΓΓ, ntrcX, strcY_d) +(1/2)*assemble(Identity(), ntrcX, strcY_d) - assemble(B21_ΩΓ, X, strcY_d)

    #error("STOP HERE!")
    
    ROW1 = hcat(B11,B12)
    ROW2 = hcat(B21,B22)
    S = vcat(ROW1,ROW2)
    R = vcat(R11,R21)

    # display(B11)
    # display(B22)
    # display(B33)

    # S*u = R*v, solve for u
    b = R*v
    u = S \ b # it solver...
    @assert norm(S*u - b) < 1e-8
    u_Φ = u[1:length(Y)]
    u_J = u[length(Y)+1:end]
    @assert length(u_Φ) == length(Y.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, nothing, u_J), S, R
end


function solvefem(; # FEM formulation
    md::meshdata,
    material::material,
    κ0::Union{Float64, Nothing} = nothing, 
    ϵ0::Union{Float64, Nothing} = nothing, # VORSICHT zsh zum echten epsilon0 nicht unbed. gegeben!
    ω::Union{Float64, Nothing} = nothing, 
    potential_top::Float64, 
    potential_bottom::Float64,
    qs3D = BEAST.SingleNumQStrat(3))

    println("2×2 block matrix - FEM formulation")

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)

    # Basis
    X = md.X
    Y = md.Y
    Y_d = md.Y_d

    # Operators

    #O = zeros(T,length(Y.fns),length(Y.fns))
    #Õ = zeros(T,length(Y.fns),length(Y_d.fns))
    O = spzeros(T, length(Y.fns),length(Y.fns))
    Õ = spzeros(T, length(Y.fns),length(Y_d.fns))
    I = assemble(BEAST.Identity(), X, X)


    Op_A = IP._grad_Ω(T(1.0), τ)
    A = assemble(Op_A, X, Y)
    Ã = -assemble(Op_A, X, Y_d)

    Op_B = IP._div_Ω(T(1.0), τ)
    B = assemble(Op_B, Y, X)

    row1 = hcat(O,B)
    row2 = hcat(A,I)
    S = vcat(row1,row2)
    R = vcat(Õ,Ã)
    
    b = R*v
    u = S \ b
    #@assert norm(S*u-b) < 1.0e-4


    u_Φ = u[1:length(Y)]
    u_Jn = nothing #u[length(y)+1:length(y)+length(w)]
    u_J = u[length(Y)+1:end]
    @assert length(u_Φ) == length(Y.fns)
    #@assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, nothing, nothing, v, b, u, u_Φ, u_Jn, u_J), S, R
end





function solve1(;   # high contrast formulation
    md::meshdata,
    material::material,
    κ0::Union{Float64, Nothing} = nothing, 
    ϵ0::Union{Float64, Nothing} = nothing, # VORSICHT zsh zum echten epsilon0 nicht unbed. gegeben!
    ω::Union{Float64, Nothing} = nothing, 
    potential_top::Float64, 
    potential_bottom::Float64,
    qs3D = BEAST.SingleNumQStrat(3), 
    qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4), #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
    qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4),
    matalloc = :center
    )

    println("3×3 block matrix - high contrast formulation")

    # if typeof(material) == general_material
        # save some points, and ...
        # material = general_material(...,...,This time not nothing here=>some test points with mat values)
    # end

    # fns spaces
    y_d = md.y_d
    y = md.y
    w = md.w
    X = md.X
    ntrc = md.ntrc

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    @show τ(p)
    @show inv_τ(p)
    @show τ0
    @show χ(p)
    @show T

    # τ(p) < 1e-12 && error("Disable the following lines...")
    # @assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) < 1e-10
    # @assert abs(1/inv_τ(p) -  τ(p)) < 1e-10

    # Material allocation to cells 
    if matalloc == :center
        cell2mat_inv_τ, cell2mat_χ = IP.gen_cell2mat(τ, inv_τ, τ0, χ, T, X) # <--- element order comes from assemblydata(X)
    elseif matalloc == :avg
        elementsX, _, _ = assemblydata(X)
        cellvolumes = CompScienceMeshes.volume.(elementsX)

        help = lagrangecxd0(md.Ω)
        helpOP1 = IP.MatId(1.0, inv_τ)
        helpOP2 = IP.MatId(1.0, χ)

        cell2mat_inv_τ = Vector(diag(assemble(helpOP1, help, help))) ./ cellvolumes
        cell2mat_χ = Vector(diag(assemble(helpOP2, help, help))) ./ cellvolumes
    else
        error("matalloc type not defined!")
    end


    X_mat = IP.gen_X_mat(X, cell2mat_χ)
    X_mat_I = IP.gen_X_mat(X, cell2mat_inv_τ)
    w_mat = IP.gen_w_mat(w, X, cell2mat_inv_τ)
    y_ctype = IP.gen_y_complextype(y, T)
    if T == Float64 
        #@assert y_ctype == y
        @assert  y_ctype.geo == y.geo
        @assert  y_ctype.fns == y.fns
        @assert  y_ctype.pos == y.pos
    end

    swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)
    intrcX_mat = IP.inner_mat_ntrace(X, swg_faces_mesh, cell2mat_χ)


    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1
    B11_Γ = IPVIE1.B11_Γ(alpha = 1.0)
    B11_ΓΓ = IPVIE1.B11_ΓΓ(alpha = 1.0, gammatype = Float64)
    B11 = assemble(B11_Γ, w, y) + assemble(B11_ΓΓ, w, y)

    UB12_ΓΓ = IPVIE1.UB12_ΓΓ(alpha = 1.0, gammatype = Float64)
    B12 = assemble(UB12_ΓΓ, w, w_mat)

    UB13_ΓΓn = IPVIE1.UB13_ΓΓn(alpha = 1.0, gammatype = Float64)
    UB13_ΓΩ = IPVIE1.UB13_ΓΩ(alpha = 1.0, gammatype = Float64)
    B13 =  assemble(UB13_ΓΓn, w, intrcX_mat) + assemble(UB13_ΓΩ, w, X_mat)


    # Operators row 2
    B21_ΓΓ = IPVIE1.B21_ΓΓ(beta = -1.0, gammatype = Float64)
    B21 = assemble(B21_ΓΓ, y, y)

    UB22_Γ = IPVIE1.UB22_Γ(alpha = -1.0)
    UB22_ΓΓ = IPVIE1.UB22_ΓΓ(alpha = 1.0, gammatype = Float64)
    B22 = assemble(UB22_Γ, y_ctype, w_mat) + assemble(UB22_ΓΓ, y, w_mat)  #  <--- T and T needed....

    UB23_ΓΓn = IPVIE1.UB23_ΓΓn(alpha = 1.0, gammatype = Float64)
    UB23_ΓΩ = IPVIE1.UB23_ΓΩ(alpha = 1.0, gammatype = Float64)
    B23 = assemble(UB23_ΓΓn, y, intrcX_mat) + assemble(UB23_ΓΩ, y, X_mat)


    # Operators row 3
    B31_ΓΓ = IPVIE1.B31_ΓΓ(alpha = 1.0, gammatype = Float64)
    B31_ΩΓ = IPVIE1.B31_ΩΓ(alpha = -1.0, gammatype = Float64)
    B31 = assemble(B31_ΓΓ, ntrc(X), y) + assemble(B31_ΩΓ, X, y)

    UB32_ΓΓ = IPVIE1.UB32_ΓΓ(alpha = 1.0, gammatype = Float64)
    UB32_ΩΓ = IPVIE1.UB32_ΩΓ(alpha = -1.0, gammatype = Float64)
    B32 = assemble(UB32_ΓΓ, ntrc(X), w_mat) + assemble(UB32_ΩΓ, X, w_mat)

    UB33_Ω = IPVIE1.UB33_Ω(alpha = -1.0)
    UB33_ΓΩ = IPVIE1.UB33_ΓΩ(alpha = 1.0, gammatype = Float64)
    UB33_ΩΩ = IPVIE1.UB33_ΩΩ(alpha = -1.0, gammatype = Float64)
    UB33_ΓΓn = IPVIE1.UB33_ΓΓn(alpha = 1.0, gammatype = Float64)
    UB33_ΩΓn = IPVIE1.UB33_ΩΓn(alpha = -1.0, gammatype = Float64)

    B33 = assemble(UB33_Ω, X, X_mat_I) +
            assemble(UB33_ΓΩ, ntrc(X), X_mat) +
            assemble(UB33_ΩΩ, X, X_mat) +
            assemble(UB33_ΓΓn, ntrc(X), intrcX_mat) + 
            assemble(UB33_ΩΓn, X, intrcX_mat)
    

    R11 = assemble(-B11_Γ, w, y_d) + assemble(-B11_ΓΓ, w, y_d)
    R21 = assemble(-B21_ΓΓ, y, y_d) 
    R31 = assemble(-B31_ΓΓ, ntrc(X), y_d) + assemble(-B31_ΩΓ, X, y_d)

    ROW1 = hcat(B11,B12,B13)
    ROW2 = hcat(B21,B22,B23)
    ROW3 = hcat(B31,B32,B33)
    S = vcat(ROW1,ROW2,ROW3)
    R = vcat(R11,R21,R31)

    # S*u = R*v, solve for u
    b = R*v
    u = S \ b # it solver?

    @assert norm(S*u - b) < 1e-5
    u_Φ = u[1:length(y)]
    u_Jn = u[length(y)+1:length(y)+length(w)]
    u_J = u[length(y)+length(w)+1:end]
    @assert length(u_Φ) == length(y.fns)
    @assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, u_Jn, u_J), S, R

end




function solve2(;   # high contrast formulation - 2 × 2 Block
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

    println("2×2 block matrix - high contrast formulation")

    # if typeof(material) == general_material
        # save some points, and ...
        # material = general_material(...,...,This time not nothing here=>some test points with mat values)
    # end

    # fns spaces
    y_d = md.y_d
    y = md.y
    w = md.w # nur zum Anzeigen von u_Jn wichtig....
    X = md.X
    ntrc = md.ntrc
    ntrcX = md.ntrcX

    

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
    p = point(0.0,0.0,0.0)
    # @show τ(p)
    # @show inv_τ(p)
    # @show τ0
    # @show χ(p)

    τ(p) < 1e-12 && error("Disable the following lines...")
    @assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) < 1e-10
    @assert abs(1/inv_τ(p) -  τ(p)) < 1e-10

    # Material -> cells, Material fns
    cell2mat_inv_τ, cell2mat_χ = IP.gen_cell2mat(τ, inv_τ, τ0, χ, T, X)

    X_mat = IP.gen_X_mat(X, cell2mat_χ)
    X_mat_I = IP.gen_X_mat(X, cell2mat_inv_τ)
    #w_mat = IP.gen_w_mat(w, X, cell2mat_inv_τ)


    swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)
    intrcX_mat = IP.inner_mat_ntrace(X, swg_faces_mesh, cell2mat_χ)
    ntrcX_mat = IP.gen_ntrcX_mat(ntrcX, cell2mat_inv_τ)

    @assert ntrcX_mat.geo ==  ntrcX.geo
    @assert ntrcX_mat.pos ==  ntrcX.pos
    @assert length(ntrcX_mat.fns) == length(ntrcX.fns)

    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1
    B11_Γ = IPVIE2.B11_Γ(alpha = 1.0)
    B11_ΓΓ = IPVIE2.B11_ΓΓ(alpha = 1.0, gammatype = T)
    B11 = assemble(B11_Γ, y, y) + assemble(B11_ΓΓ, y, y)

    UB12_ΓΓ = IPVIE2.UB12_ΓΓ(alpha = 1.0, gammatype = T)
    UB12_ΓΩ = IPVIE2.UB12_ΓΩ(alpha = 1.0, gammatype = T)
    UB12_ΓΓn = IPVIE2.UB12_ΓΓn(alpha = 1.0, gammatype = T)
    B12 = assemble(UB12_ΓΓ, y, ntrcX_mat) + assemble(UB12_ΓΩ, y, X_mat) + assemble(UB12_ΓΓn, y, intrcX_mat)


    # Operators row 2
    B21_ΓΓ = IPVIE2.B21_ΓΓ(alpha = 1.0, gammatype = T)
    B21_ΩΓ = IPVIE2.B21_ΩΓ(alpha = -1.0, gammatype = T)
    B21 = assemble(B21_ΓΓ, ntrc(X), y) + assemble(B21_ΩΓ, X, y)

    UB22_Ω = IPVIE2.UB22_Ω(alpha = -1.0)
    UB22_ΓΓ = IPVIE2.UB22_ΓΓ(alpha = 1.0, gammatype = T)
    UB22_ΩΓ = IPVIE2.UB22_ΩΓ(alpha = -1.0, gammatype = T)
    UB22_ΓΩ = IPVIE2.UB22_ΓΩ(alpha = 1.0, gammatype = T)
    UB22_ΩΩ = IPVIE2.UB22_ΩΩ(alpha = -1.0, gammatype = T)
    UB22_ΓΓn = IPVIE2.UB22_ΓΓn(alpha = 1.0, gammatype = T)
    UB22_ΩΓn = IPVIE2.UB22_ΩΓn(alpha = -1.0, gammatype = T)

    B22 = assemble(UB22_Ω, X, X_mat_I) +
            assemble(UB22_ΓΓ, ntrc(X), ntrcX_mat) + 
            assemble(UB22_ΩΓ, X, ntrcX_mat) +
            assemble(UB22_ΓΩ, ntrc(X), X_mat) +
            assemble(UB22_ΩΩ, X, X_mat) +
            assemble(UB22_ΓΓn, ntrc(X), intrcX_mat) + 
            assemble(UB22_ΩΓn, X, intrcX_mat)
    

    R11 = assemble(-B11_Γ, y, y_d) + assemble(-B11_ΓΓ, y, y_d)
    R21 = assemble(-B21_ΓΓ, ntrc(X), y_d) + assemble(-B21_ΩΓ, X, y_d)

    ROW1 = hcat(B11,B12)
    ROW2 = hcat(B21,B22)
    S = vcat(ROW1,ROW2)
    R = vcat(R11,R21)

    # S*u = R*v, solve for u
    b = R*v
    u = S \ b # it solver?

    @show norm(S*u - b)

    u_Φ = u[1:length(y)]
    u_J = u[length(y)+1:end]

    @assert X.pos == ntrcX.pos # => X.pos[1] == ntrcX.pos[1], ...
    #u_Jn_not_w_ready = Vector{T}()
    index_list = Vector{Int64}()
    for (i,fns) in enumerate(ntrcX.fns)
        fns == [] && continue
        #push!(u_Jn_not_w_ready, u_J[i])
        push!(index_list, i)
    end

    u_Jn = Vector{T}()
    for pw in w.pos
        for index in index_list #ntrcX...
            pn = ntrcX.pos[index]
            if norm(pw-pn) < 1.0e-10
                push!(u_Jn, u_J[index])
            end
        end
    end

    @assert length(u_Φ) == length(y.fns)
    @assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, u_Jn, u_J), S, R




end






function solve3(;   # arb. Material formulation - 2 × 2 Block
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

    println("2×2 block matrix - arb. material formulation (not hc)")

    # if typeof(material) == general_material
        # save some points, and ...
        # material = general_material(...,...,This time not nothing here=>some test points with mat values)
    # end

    # fns spaces
    y_d = md.y_d
    y = md.y
    w = md.w # nur zum Anzeigen von u_Jn wichtig....
    X = md.X
    ntrc = md.ntrc
    ntrcX = md.ntrcX

    

    # Material
    κ, ϵ = material()
    τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)


    # Excitation
    v_top = ones(length(md.topnodes)) * potential_top
    v_bottom = ones(length(md.bottomnodes)) * potential_bottom
    v = vcat(v_top, v_bottom)


    # Operators row 1
    B11_Γ = IPVIE3.B11_Γ(alpha = 1.0)
    B11_ΓΓ = IPVIE3.B11_ΓΓ(alpha = 1.0, gammatype = T)
    B11 = assemble(B11_Γ, y, y) + assemble(B11_ΓΓ, y, y)

    B12_ΓΓ = IPVIE3.B12_ΓΓ(alpha = 1.0, gammatype = T, invtau0 = 1/τ0)
    B12_ΓΩ = IPVIE3.B12_ΓΩ(alpha = 1.0, gammatype = T, chi = χ)
    B12 = assemble(B12_ΓΓ, y, ntrcX) + assemble(B12_ΓΩ, y, X)


    # Operators row 2
    B21_ΓΓ = IPVIE3.B21_ΓΓ(alpha = 1.0, gammatype = T)
    B21_ΩΓ = IPVIE3.B21_ΩΓ(alpha = -1.0, gammatype = T)
    B21 = assemble(B21_ΓΓ, ntrcX, y) + assemble(B21_ΩΓ, X, y)

    B22_Ω = IPVIE3.B22_Ω(alpha = -1.0, invtau = inv_τ)
    B22_ΓΓ = IPVIE3.B22_ΓΓ(alpha = 1.0, gammatype = T, invtau0 = 1/τ0)
    B22_ΓΩ = IPVIE3.B22_ΓΩ(alpha = -1.0, gammatype = T, chi = χ)
    B22_ΩΓ = IPVIE3.B22_ΩΓ(alpha = -1.0, gammatype = T, invtau0 = 1/τ0)
    B22_ΩΩ = IPVIE3.B22_ΩΩ(alpha = 1.0, gammatype = T, chi = χ)

    B22 = assemble(B22_Ω, X, X) + 
            assemble(B22_ΓΓ, ntrcX, ntrcX) + 
            assemble(B22_ΓΩ, ntrcX, X) +
            assemble(B22_ΩΓ, X, ntrcX) +
            assemble(B22_ΩΩ, X, X)
    
    R11 = assemble(-B11_Γ, y, y_d) + assemble(-B11_ΓΓ, y, y_d)
    R21 = assemble(-B21_ΓΓ, ntrcX, y_d) + assemble(-B21_ΩΓ, X, y_d)

    ROW1 = hcat(B11,B12)
    ROW2 = hcat(B21,B22)
    S = vcat(ROW1,ROW2)
    R = vcat(R11,R21)

    # S*u = R*v, solve for u
    b = R*v
    u = S \ b # it solver?

    @show norm(S*u - b)

    u_Φ = u[1:length(y)]
    u_J = u[length(y)+1:end]

    @assert X.pos == ntrcX.pos # => X.pos[1] == ntrcX.pos[1], ...
    index_list = Vector{Int64}()
    for (i,fns) in enumerate(ntrcX.fns)
        fns == [] && continue
        push!(index_list, i)
    end
    u_Jn = Vector{T}()
    for pw in w.pos
        for index in index_list #ntrcX...
            pn = ntrcX.pos[index]
            if norm(pw-pn) < 1.0e-10
                push!(u_Jn, u_J[index])
            end
        end
    end

    @assert length(u_Φ) == length(y.fns)
    @assert length(u_Jn) == length(w.fns)
    @assert length(u_J) == length(X.fns)

    return solution(material, κ0, ϵ0, ω, τ0, potential_top, potential_bottom, qs3D, qs4D, qs5D6D, v, b, u, u_Φ, u_Jn, u_J), S, R
end








# getcurrent 2 => mittels ntrc!
function getcurrent2(m::IP.meshdata, s::IP.solution)

    Γ_c_t = m.Γ_c_t
    Γ_c_b = m.Γ_c_b
    u_J = s.u_J # not u_Jn!

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

    ntrcX = m.ntrcX

    for (j, pos) in enumerate(ntrcX.pos)
        ntrcX.fns[j] == [] && continue

        i = CompScienceMeshes.findchart(top_charts, chart_tree_top, pos)
        k = CompScienceMeshes.findchart(bottom_charts, chart_tree_bottom, pos)

        if i !== nothing
            A = top_charts[i].volume 
            @assert A > 0.0
            I_part = -A * u_J[j] * ntrcX.fns[j][1].coeff # "-" because dÂ of ∫∫J_vec*dÂ in opposite dir
            I_top += I_part
            @assert length(ntrcX.fns[j]) == 1
            @assert I_part >= -0.0
            #@show -A * u_J[j] * ntrcX.fns[j][1].coeff
            cnt_top += 1
        elseif k !== nothing
            A = bottom_charts[k].volume 
            @assert A > 0.0
            I_part = A * u_J[j] * ntrcX.fns[j][1].coeff
            I_bottom += I_part
            @assert length(ntrcX.fns[j]) == 1
            @assert I_part >= -0.0
            #@show A * u_J[j] * ntrcX.fns[j][1].coeff
            cnt_bottom += 1
        else
            error("Neither top nor bottom")
        end
    
    end

    @assert cnt_top == numtopcharts
    @assert cnt_bottom == numbottomcharts

    return I_top, I_bottom
end




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
            @assert top_charts[i].normals[1][3] > 0.0
            @assert A > 0.0
            I_top += -A * u_Jn[j] * w.fns[j][1].coeff # "-" because dÂ of ∫∫J_vec*dÂ in opposite dir
            cnt_top += 1
        elseif k !== nothing
            A = bottom_charts[k].volume
            @assert bottom_charts[k].normals[1][3] < 0.0
            @assert A > 0.0
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