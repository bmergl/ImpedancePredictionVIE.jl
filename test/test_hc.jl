using LinearAlgebra
using StaticArrays
using SparseArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using Test

md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.045)
print("tehrahedrons: ", length(md.Ω.faces))

# Select Test Material
material = IP.constant_zsplit(10.0, nothing, 0.0, 2.0, nothing) #IP.constantmaterial(1.0, nothing),# IP.constant_xsplit(10.0, nothing, 0.0, 5.0, nothing), ,   #IP.constant_xsplit(1.0, nothing, 0.1, 1.0, nothing), #IP.constant_zsplit(1000, nothing, 0.07, 100, nothing) #IP.constantmaterial(100.0, nothing), 
κ0 = 1.0
ϵ0 = nothing
ω = nothing 


# fns spaces
y_d = md.y_d
y = md.y
w = md.w
X = md.X
ntrc = md.ntrc 

# Material
κ, ϵ = material()
τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)
p = point(0.0,0.0,0.3)
@show τ(p)
@show inv_τ(p)
@show τ0
@show χ(p)

τ(p) < 1e-12 && error("Disable the following lines...")
@assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) < 1e-10
@assert abs(1/inv_τ(p) -  τ(p)) < 1e-10

# Material -> cells, Material fns
cell2mat_inv_τ, cell2mat_χ = IP.gen_cell2mat(τ, inv_τ, τ0, χ, T, X)

p_up = point(0,0,0.1)
p_down = point(0,0,-0.1)
inv_τ(p_up)
χ(p_up)
inv_τ(p_down)
χ(p_down)

elementsX, adX, cellsX = assemblydata(X)
for (n,el) in enumerate(elementsX)
    center = cartesian(CompScienceMeshes.center(el))
    if center[3] >= 0.0
        @test cell2mat_inv_τ[n] == inv_τ(p_up)
        @test cell2mat_χ[n] == χ(p_up)
    elseif center[3] < 0.0
        @test cell2mat_inv_τ[n] == inv_τ(p_down)
        @test cell2mat_χ[n] == χ(p_down)
    end
end

X_mat = IP.gen_X_mat(X, cell2mat_χ)
X_mat_ = IP.gen_X_mat(X, cell2mat_inv_τ)
w_mat = IP.gen_w_mat(w, X, cell2mat_inv_τ)


for (n,fns) in enumerate(w_mat.fns)
    pos = w_mat.pos[n]
    if pos[3] >= 0.0
        @test fns[1].coeff / w.fns[n][1].coeff ≈ inv_τ(p_up)
    elseif pos[3] < 0.0
        @test fns[1].coeff / w.fns[n][1].coeff ≈ inv_τ(p_down)
    end
end

swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)
intrcX_mat = IP.inner_mat_ntrace(X, swg_faces_mesh, cell2mat_χ)

# intrcX_mat check (constant_zsplit)
for (n,fns) in enumerate(intrcX_mat.fns)

    pos = intrcX_mat.pos[n]
    
    if pos[3] > 0.05
        @test fns == []
        continue
    elseif pos[3] < -0.05
        @test fns == []
        continue
    end

end
newfns = Vector{Vector{BEAST.Shape{Float64}}}() # remove [] and corresponding pos
newpos = Vector{SVector{3, Float64}}()
newfaces = Vector{SVector{3, Int}}()
elements, ad, nr = assemblydata(intrcX_mat)

for (index,n) in enumerate(nr) # existente
    @assert cartesian(CompScienceMeshes.center(elements[index])) ≈ intrcX_mat.pos[n]
    push!(newfaces, intrcX_mat.geo.supermesh.faces[n])
    fnsn = intrcX_mat.fns[n][1]
    
    newshs = [BEAST.Shape(index ,fnsn.refid, fnsn.coeff)]
    push!(newfns, newshs)
    push!(newpos, intrcX_mat.pos[n])
end
newmesh = Mesh(intrcX_mat.geo.supermesh.vertices, newfaces)
intrcX_mat_reduced = BEAST.LagrangeBasis{0,-1,1}(newmesh, newfns, newpos)

cf = ones(length(intrcX_mat_reduced))
fcr0, geo0 = facecurrents(cf, intrcX_mat_reduced)
Plotly.plot(patch(geo0, fcr0))



##


##
t
##








## assemble comparison IPVIE vs. IPVIE1 - constant_zsplit















## assemble test with constant material

# Select Test Material
material = IP.constantmaterial(1.0, nothing) # IP.constant_xsplit(10.0, nothing, 0.0, 5.0, nothing), ,   #IP.constant_xsplit(1.0, nothing, 0.1, 1.0, nothing), #IP.constant_zsplit(1000, nothing, 0.07, 100, nothing) #IP.constantmaterial(100.0, nothing), 
κ0 = 1.0
ϵ0 = nothing
ω = nothing 


# Material
κ, ϵ = material()
τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω)

# Material -> cells, Material fns
cell2mat_inv_τ, cell2mat_χ = IP.gen_cell2mat(τ, inv_τ, τ0, χ, T, X)





# overwrite cell2mat
cell2mat_inv_τ = 2*ones(length(cell2mat_inv_τ))
cell2mat_χ = 2*ones(length(cell2mat_χ))
w_mat = IP.gen_w_mat(w, X, cell2mat_inv_τ)
X_mat = IP.gen_X_mat(X, cell2mat_χ)
X_mat_ = IP.gen_X_mat(X, cell2mat_inv_τ)

UB13_ΓΩ = IPVIE1.UB13_ΓΩ(gammatype=Float64)
@test norm(2*assemble(UB13_ΓΩ, w, X) - assemble(UB13_ΓΩ, w, X_mat)) < 10e-14

UB22_ΓΓ = IPVIE1.UB22_ΓΓ(gammatype=Float64) 
@test norm(2*assemble(UB22_ΓΓ, y, w) - assemble(UB22_ΓΓ, y, w_mat)) < 10e-14

UB33_ΓΩ = IPVIE1.UB33_ΓΩ(gammatype=Float64)
@test norm(2*assemble(UB33_ΓΩ, ntrc(X), X) - assemble(UB33_ΓΩ, ntrc(X), X_mat)) < 10e-14




norm(assemble(UB13_ΓΩ, y_d, w))


elementsX, adX, cellsX = assemblydata(X)
for (n,el) in enumerate(elementsX)
    center = cartesian(CompScienceMeshes.center(el))
    if center[3] >= 0.0
        @test cell2mat_inv_τ[n] == inv_τ(p_up)
        @test cell2mat_χ[n] == χ(p_up)
    elseif center[3] < 0.0
        @test cell2mat_inv_τ[n] == inv_τ(p_down)
        @test cell2mat_χ[n] == χ(p_down)
    end
end






for (n,fns) in enumerate(w_mat.fns)
    pos = w_mat.pos[n]
    if pos[3] >= 0.0
        @test fns[1].coeff == inv_τ(p_up)
    elseif pos[3] < 0.0
        @test fns[1].coeff == inv_τ(p_down)
    end
end



swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)

intrcX_mat = IP.inner_mat_ntrace(X, swg_faces_mesh, cell2mat_χ)




# for (n,fns) in enumerate(intrcX_mat.fns)
#     fns == [] && continue
#     @show fns[1] 
# end

# (χ(p_up)-χ(p_down))*49.0

# ntrace(X,md.Γ)


# ntrc(X).fns
# for (n,fns) in enumerate(ntrc(X).fns)
#     fns == [] && continue
#     @show fns[1] 
# end

# intrcX_mat_test = IP.inner_mat_ntrace(X, swg_faces_mesh, ones(length(cell2mat_χ)))
# intrcX_mat_test.fns

# for (n,fns) in enumerate(intrcX_mat_test.fns)
#     fns == [] && continue
#     @show fns[1] 
# end

