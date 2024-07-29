using ..ImpedancePredictionVIE
using CompScienceMeshes
using BEAST
using LinearAlgebra
using StaticArrays
using SparseArrays

##### constantmaterial - analytic solution ####################

function solution_J_ana(body::cuboid, mat::constantmaterial, m::IP.meshdata, s::IP.solution, points, J_MoM)
    
    if mat.ϵ !== nothing && mat.κ !== nothing
        L_z = body.L_z
        ω = s.ω
        ϵ = mat.ϵ
        κ = mat.κ
        U = s.potential_top-s.potential_bottom
        Jz_ana = -(κ + im*ω*ϵ)*U/L_z
        J_ana = fill(SVector{3, ComplexF64}(0.0, 0.0, Jz_ana), length(J_MoM))

        return J_ana
    end


    Jz_ana = -mat.κ*(s.potential_top-s.potential_bottom)/body.L_z

    J_ana = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM))

    return J_ana
end

function solution_I_ana(body::cuboid, mat::constantmaterial, m::IP.meshdata, s::IP.solution)
    
    if mat.ϵ !== nothing && mat.κ !== nothing
        L_z = body.L_z
        L_x = body.L_x
        L_y = body.L_y
        ω = s.ω
        ϵ = mat.ϵ
        κ = mat.κ
        U = s.potential_top-s.potential_bottom
        R = (1/κ)*L_z/(L_x*L_y)
        C = ϵ*(L_x*L_y)/L_z
        #@show C
        Z = (R - im*ω*R^2*C)/(1 + ω^2*R^2*C^2)

        I_ana = U/Z

        return I_ana
    end

    I_ana = mat.κ*body.L_x*body.L_y*(1/body.L_z)*(s.potential_top-s.potential_bottom)

    return I_ana
end

function solution_Φ_ana(body::cuboid, mat::constantmaterial, m::IP.meshdata, s::IP.solution; FEM = false)
    u_Φ = s.u_Φ

    function CapacitorPot(x, v_top, v_bottom) # hom. medium
        L_z = body.L_z
        z = x[3]

        return ((v_top-v_bottom)/L_z) * (z + L_z/2) + v_bottom
    end
    u_Φ_ana = Vector{Float64}(undef,length(u_Φ))
    y = m.y
    FEM == true && (y = m.Y)
    for (i, pos) in enumerate(y.pos)
        u_Φ_ana[i] = CapacitorPot(pos, s.potential_top, s.potential_bottom)
    end

    return u_Φ_ana
end


##### constant_zsplit Material - analytic solution ##################

function solution_J_ana(body::IP.cuboid, mat::IP.constant_zsplit, m::IP.meshdata, s::IP.solution, points, J_MoM)
    
    if mat.ϵ_p !== nothing && mat.κ_p !== nothing && mat.ϵ_m !== nothing && mat.κ_m !== nothing
        ϵ_p = mat.ϵ_p
        κ_p = mat.κ_p
        ϵ_m = mat.ϵ_m
        κ_m = mat.κ_m
        A = body.L_x * body.L_y
        lz = body.L_z
        ω = s.ω
        z0 = mat.z0
        l_p = lz/2 - z0
        l_m = lz/2 + z0

        Z1 = 1/(κ_p*A/l_p + im*ω*ϵ_p*A/l_p)
        Z2 = 1/(κ_m*A/l_m + im*ω*ϵ_m*A/l_m)
        Z = Z1 + Z2

        U = s.potential_top-s.potential_bottom
        
        Jz_ana = -U/(Z*A)

        J_ana = fill(SVector{3, ComplexF64}(0.0, 0.0, Jz_ana), length(J_MoM))
        return J_ana
    end

    A = body.L_x * body.L_y
    R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    U = s.potential_top - s.potential_bottom
    I = U/R
    Jz_ana = -I/A 

    J_ana = fill(SVector{3, Float64}(0.0, 0.0, Jz_ana), length(J_MoM))

    return J_ana
end

function solution_I_ana(body::IP.cuboid, mat::IP.constant_zsplit, m::IP.meshdata, s::IP.solution)

    if mat.ϵ_p !== nothing && mat.κ_p !== nothing && mat.ϵ_m !== nothing && mat.κ_m !== nothing
        ϵ_p = mat.ϵ_p
        κ_p = mat.κ_p
        ϵ_m = mat.ϵ_m
        κ_m = mat.κ_m
        A = body.L_x * body.L_y
        lz = body.L_z
        ω = s.ω

        z0 = mat.z0
        l_p = lz/2 - z0
        l_m = lz/2 + z0

        Z1 = 1/(κ_p*A/l_p + im*ω*ϵ_p*A/l_p)
        Z2 = 1/(κ_m*A/l_m + im*ω*ϵ_m*A/l_m)
        Z = Z1 + Z2

        U = s.potential_top-s.potential_bottom
        
        I_ana = U/Z

        return I_ana
    end

    A = body.L_x * body.L_y
    R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
    U = s.potential_top - s.potential_bottom
    I = U/R

    I_ana = I

    return I_ana
end

function solution_Φ_ana(body::IP.cuboid, mat::IP.constant_zsplit, m::IP.meshdata, s::IP.solution; FEM = false)
    u_Φ = s.u_Φ

    function linPot_z(x, z1, z2, Φ1, Φ2) # linear in section [z1,z2] -> [Φ1, Φ2]
        z = x[3] # obervation plane

        return ((Φ2 - Φ1)/(z2 - z1))*(z - z1) + Φ1
    end

    
    if mat.ϵ_p !== nothing && mat.κ_p !== nothing && mat.ϵ_m !== nothing && mat.κ_m !== nothing
        ϵ_p = mat.ϵ_p
        κ_p = mat.κ_p
        ϵ_m = mat.ϵ_m
        κ_m = mat.κ_m
        A = body.L_x * body.L_y
        lz = body.L_z
        ω = s.ω

        Z1 = 1/(κ_p*A/lz + im*ω*ϵ_p*A/lz)
        Z2 = 1/(κ_m*A/lz + im*ω*ϵ_m*A/lz)
        Z = Z1 + Z2

        U = s.potential_top-s.potential_bottom
        U_p = (Z1/Z)*U
        U_m = (Z2/Z)*U
        u_Φ_ana = Vector{ComplexF64}(undef,length(u_Φ))
    else 
        A = body.L_x * body.L_y
        R = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A + (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
        U = s.potential_top - s.potential_bottom
        R_p = (1/mat.κ_p)*(body.L_z/2 - mat.z0)/A
        R_m = (1/mat.κ_m)*(body.L_z/2 + mat.z0)/A
        U_p = U * R_p/R
        U_m = U * R_m/R
        u_Φ_ana = Vector{Float64}(undef,length(u_Φ))
    end
    
    
    z0 = mat.z0
    y = m.y
    FEM == true && (y = m.Y)
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

function solution_J_ana(body::IP.cuboid, mat::IP.constant_xsplit, m::IP.meshdata, s::IP.solution, points, J_MoM)
    
    if mat.ϵ_p !== nothing && mat.κ_p !== nothing && mat.ϵ_m !== nothing && mat.κ_m !== nothing
        ϵ_p = mat.ϵ_p
        κ_p = mat.κ_p
        ϵ_m = mat.ϵ_m
        κ_m = mat.κ_m

        x0 = mat.x0
        A1 = (body.L_x/2 - x0) * body.L_y
        A2 = (body.L_x/2 + x0) * body.L_y
        lz = body.L_z

        ω = s.ω

        #invZ = κ_p*A1/lz + im*ω*ϵ_p*A1/lz + κ_m*A2/lz + im*ω*ϵ_m*A2/lz
        invZ1 = κ_p*A1/lz + im*ω*ϵ_p*A1/lz
        invZ2 = κ_m*A2/lz + im*ω*ϵ_m*A2/lz
        
        U = s.potential_top-s.potential_bottom

        Jz_p_ana = -(U/A1)*invZ1
        Jz_m_ana = -(U/A2)*invZ2

        J_ana = Vector{SVector{3, ComplexF64}}(undef, length(J_MoM)) #(0.0, 0.0, Jz_ana), length(J_MoM))
        for (i,pos) in enumerate(points)
            x = pos[1]
            if x >= x0
                J_ana[i] = SVector{3, ComplexF64}(0.0, 0.0, Jz_p_ana)
            elseif x < x0
                J_ana[i] = SVector{3, ComplexF64}(0.0, 0.0, Jz_m_ana)
            else
                error("")
            end
        end
    
        return J_ana

    else

        x0 = mat.x0
        A_p = body.L_y * (body.L_x/2 - x0)
        A_m = body.L_y * (body.L_x/2 + x0)

        R_p = (1/mat.κ_p)*body.L_z/A_p
        R_m = (1/mat.κ_m)*body.L_z/A_m

        U = s.potential_top - s.potential_bottom
        I_p = U/R_p
        I_m = U/R_m

        Jz_p_ana = -I_p/A_p
        Jz_m_ana = -I_m/A_m 
    
    end

    @assert length(points) == length(J_MoM)

    J_ana = Vector{SVector{3, Float64}}(undef, length(J_MoM)) #(0.0, 0.0, Jz_ana), length(J_MoM))
    for (i,pos) in enumerate(points)
        x = pos[1]
        if x >= x0
            J_ana[i] = SVector{3, Float64}(0.0, 0.0, Jz_p_ana)
        elseif x < x0
            J_ana[i] = SVector{3, Float64}(0.0, 0.0, Jz_m_ana)
        else
            error("")
        end
    end

    return J_ana
end

function solution_I_ana(body::IP.cuboid, mat::IP.constant_xsplit, m::IP.meshdata, s::IP.solution)

    if mat.ϵ_p !== nothing && mat.κ_p !== nothing && mat.ϵ_m !== nothing && mat.κ_m !== nothing
        ϵ_p = mat.ϵ_p
        κ_p = mat.κ_p
        ϵ_m = mat.ϵ_m
        κ_m = mat.κ_m

        x0 = mat.x0
        A1 = (body.L_x/2 - x0) * body.L_y
        A2 = (body.L_x/2 + x0) * body.L_y
        lz = body.L_z

        ω = s.ω

        invZ = κ_p*A1/lz + im*ω*ϵ_p*A1/lz + κ_m*A2/lz + im*ω*ϵ_m*A2/lz
    
        U = s.potential_top-s.potential_bottom

        I_ana = U*invZ
    
        return I_ana
    end

    x0 = mat.x0
    A_p = body.L_y * (body.L_x/2 - x0)
    A_m = body.L_y * (body.L_x/2 + x0)

    R_p = (1/mat.κ_p)*body.L_z/A_p
    R_m = (1/mat.κ_m)*body.L_z/A_m

    U = s.potential_top - s.potential_bottom
    I_p = U/R_p
    I_m = U/R_m
    I = I_p + I_m

    I_ana = I

    return I_ana
end

function solution_Φ_ana(body::IP.cuboid, mat::IP.constant_xsplit, m::IP.meshdata, s::IP.solution; FEM = false)
    u_Φ = s.u_Φ

    function linPot_z(x, z1, z2, Φ1, Φ2) # linear in section [z1,z2] -> [Φ1, Φ2]
        z = x[3] # obervation plane

        return ((Φ2 - Φ1)/(z2 - z1))*(z - z1) + Φ1
    end

    u_Φ_ana = Vector{Any}(undef,length(u_Φ))
    y = m.y
    FEM == true && (y = m.Y)

    z1 = -body.L_z/2
    Φ1 = s.potential_bottom
    z2 = body.L_z/2
    Φ2 = s.potential_top

    for (i, pos) in enumerate(y.pos)
        u_Φ_ana[i] = linPot_z(pos, z1, z2, Φ1, Φ2)
    end

    return u_Φ_ana
end


###### pwlinx - analytic solution ##################

function solution_J_ana(body::IP.cuboid, mat::IP.pwlinx, m::IP.meshdata, s::IP.solution, points, J_MoM)
    
    κ, ϵ = mat() # Funktionen

    U = s.potential_top - s.potential_bottom
    Lz = body.L_z

    @assert length(points) == length(J_MoM)

    J_ana = Vector{SVector{3, Float64}}(undef, length(J_MoM))
    for (i,pos) in enumerate(points)
        J_ana[i] = SVector{3, Float64}(0.0, 0.0, -κ(pos)*U/Lz)
    end

    return J_ana
end

function solution_I_ana(body::IP.cuboid, mat::IP.pwlinx, m::IP.meshdata, s::IP.solution)

    U = s.potential_top - s.potential_bottom
    Ly = body.L_y
    Lz = body.L_z
    xpos = mat.xpos
    κ_vec = mat.κ

    I_ana = 0.0
    for i in 1:length(xpos)-1
        x1 = xpos[i]
        x2 = xpos[i+1]
        I_part = U*(Ly/Lz)*(κ_vec[i][2]+κ_vec[i][1])*(x2-x1)/2
        @assert I_part >= 0.0
        I_ana += I_part
    end

    return I_ana
end

function solution_Φ_ana(body::IP.cuboid, mat::IP.pwlinx, m::IP.meshdata, s::IP.solution; FEM = false)
    u_Φ = s.u_Φ

    function linPot_z(x, z1, z2, Φ1, Φ2) # linear in section [z1,z2] -> [Φ1, Φ2]
        z = x[3] # obervation plane

        return ((Φ2 - Φ1)/(z2 - z1))*(z - z1) + Φ1
    end

    u_Φ_ana = Vector{Float64}(undef,length(u_Φ))
    y = m.y
    FEM == true && (y = m.Y)

    z1 = -body.L_z/2
    Φ1 = s.potential_bottom
    z2 = body.L_z/2
    Φ2 = s.potential_top

    for (i, pos) in enumerate(y.pos)
        u_Φ_ana[i] = linPot_z(pos, z1, z2, Φ1, Φ2)
    end

    return u_Φ_ana
end
