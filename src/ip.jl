module IP

using ..ImpedancePredictionVIE #damit alles exportierte hier funktioniert

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
end

struct filling


end

# was brauchen wir noch zugang zu den Materialparametern nachdem die Lösung da ist und mit 
# JLD2 gesichert wurde? Ja irgendwie schon... 
# Achtung irgendwann kommt material ja nicht mehr über die Funktion

struct solution

    κ::Function
    κ0::Union{Float64, Nothing}
    ϵ::Function
    ϵ0::Union{Float64, Nothing}
    ω::Union{Float64, Nothing}
    qs3D
    qs4D
    qs5D6D

    τ0::Union{Float64, ComplexF64}
    τ
    inv_τ
    χ

    R
    v
    b
    S

    u
    u_Φ
    u_Jn
    u_J

end



end















# struct fnsdata


# end


# mutable struct 