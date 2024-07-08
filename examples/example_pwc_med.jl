using MKL

using LinearAlgebra
using StaticArrays
using SparseArrays
using ImpedancePredictionVIE
using BEAST
using CompScienceMeshes

using JLD2

using Plots
using Plotly


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.12)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Γ) 

y_d = md.y_d
y = md.y
w = md.w
X = md.X

ntrc = md.ntrc


testmat = IP.constant_xsplit(100000.0, nothing, 0.0, 0.001, nothing)
κ, ϵ = testmat() # hier noch Funktionen
κ0 = 100.0
ϵ0 = nothing
ω = nothing

τ, inv_τ, τ0, χ, T = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω) # hier noch Funktionen

cell2mat_inv_τ, cell2mat_χ = IP.gen_cell2mat(τ, inv_τ, τ0, χ, T, X)

X_mat = IP.gen_X_mat(X, cell2mat_χ)
X_mat_ = IP.gen_X_mat(X, cell2mat_inv_τ)
w_mat = IP.gen_w_mat(w, X, cell2mat_inv_τ)

swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)
intrcX_mat = IP.inner_mat_ntrace(X, swg_faces_mesh, cell2mat_χ)


# Quadstrat
qs3D = BEAST.SingleNumQStrat(3)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,4,4,4,4) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(3,3,4,4,4,4)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D


# Operators row 1
B11_Γ = IPVIE1.B11_Γ()
B11_ΓΓ = IPVIE1.B11_ΓΓ(gammatype = Float64)
B11 = assemble(B11_Γ, w, y) + assemble(B11_ΓΓ, w, y)

UB12_ΓΓ = IPVIE1.UB12_ΓΓ(gammatype = Float64)
B12 = assemble(UB12_ΓΓ, w, w_mat)

UB13_ΓΓn = IPVIE1.UB13_ΓΓn(gammatype = Float64)
UB13_ΓΩ = IPVIE1.UB13_ΓΩ(gammatype = Float64)
B13 = assemble(UB13_ΓΓn, w, intrcX_mat) + assemble(UB13_ΓΩ, w, X_mat)


# Operators row 2
B21_ΓΓ = IPVIE1.B21_ΓΓ(gammatype = Float64)
B21 = assemble(B21_ΓΓ, y, y)

UB22_Γ = IPVIE1.UB22_Γ()
UB22_ΓΓ = IPVIE1.UB22_ΓΓ(gammatype = Float64)
B22 = assemble(UB22_Γ, y, w_mat) + assemble(UB22_ΓΓ, y, w_mat)

UB23_ΓΓn = IPVIE1.UB23_ΓΓn(gammatype = Float64)
UB23_ΓΩ = IPVIE1.UB23_ΓΩ(gammatype = Float64)
B23 = assemble(UB23_ΓΓn, y, intrcX_mat) + assemble(UB23_ΓΩ, y, X_mat)


# Operators row 3
B31_ΓΓ = IPVIE1.B31_ΓΓ(gammatype = Float64)
B31_ΩΓ = IPVIE1.B31_ΩΓ(gammatype = Float64)
B31 = assemble(B31_ΓΓ, ntrc(X), y) + assemble(B31_ΩΓ, X, y)

UB32_ΓΓ = IPVIE1.UB32_ΓΓ(gammatype = Float64)
UB32_ΩΓ = IPVIE1.UB32_ΩΓ(gammatype = Float64)
B32 = assemble(UB32_ΓΓ, ntrc(X), w_mat) + assemble(UB32_ΩΓ, X, w_mat)

UB33_Ω = IPVIE1.UB33_Ω()
UB33_ΓΓn = IPVIE1.UB33_ΓΓn(gammatype = Float64)
UB33_ΓΩ = IPVIE1.UB33_ΓΩ(gammatype = Float64)
UB33_ΩΓn = IPVIE1.UB33_ΩΓn(gammatype = Float64)
UB33_ΩΩ = IPVIE1.UB33_ΩΩ(gammatype = Float64)
B33 = assemble(UB33_Ω, X, X_mat_) +
        assemble(UB33_ΓΓn, ntrc(X), intrcX_mat) + 
        assemble(UB33_ΓΩ, ntrc(X), X_mat) +
        assemble(UB33_ΩΓn, X, intrcX_mat) + 
        assemble(UB33_ΩΩ, X, X_mat)

R11 = assemble(-B11_Γ, w, y_d) + assemble(-B11_ΓΓ, w, y_d)
R21 = assemble(-B21_ΓΓ, y, y_d) 
R31 = assemble(-B31_ΓΓ, ntrc(X), y_d) + assemble(-B31_ΩΓ, X, y_d)


##
potential_top = 0.5
potential_bottom = -0.5

# Excitation
v_top = ones(length(md.topnodes)) * potential_top
v_bottom = ones(length(md.bottomnodes)) * potential_bottom
v = vcat(v_top, v_bottom)

ROW1 = hcat(B11,B12,B13)
ROW2 = hcat(B21,B22,B23)
ROW3 = hcat(B31,B32,B33)
S = vcat(ROW1,ROW2,ROW3)
R = vcat(R11,R21,R31)

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

##

# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, u_J, X)#, type=Float64)
J_ana = IP.solution_J_ana(md.body, sol.material, md, sol, points, J_MoM)
display("Stomdichte Gesamtvolumen")
@show norm(J_MoM-J_ana)/norm(J_ana)# = norm(norm.(J_MoM-J_ana))/norm(J_ana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, sol.u_J, md.X)
J_ana2 = IP.solution_J_ana(md.body, sol.material, md, sol, points2, J_MoM2)
display("Stromdichte Mitte: Ebene z=0.0")
@show norm(J_MoM2-J_ana2)/norm(J_ana2)

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, sol.u_J, md.X)
J_ana3 = IP.solution_J_ana(md.body, sol.material, md, sol, points3, J_MoM3)
display("Stromdichte bei Platten: Ebene z=0.49")
@show norm(J_MoM3-J_ana3)/norm(J_ana3)

# Strom durch Platten
display("")
I_top, I_bottom = getcurrent(md, sol)
#@show I_top, I_bottom
I_ana = IP.solution_I_ana(md.body, sol.material, md, sol)
@show norm(I_top-I_ana)/norm(I_ana)
@show norm(I_bottom-I_ana)/norm(I_ana)
display("")

# Potential: Randknoten vs. Analytisch
u_Φ = sol.u_Φ
u_Φ_ana = IP.solution_Φ_ana(md.body, sol.material, md, sol)
@show norm(u_Φ-u_Φ_ana)/norm(u_Φ_ana)


##
display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(md.Γ_c)))














##  OLD: ######################################

dasda
    # for (p,el) in enumerate(E)

    #     for (q,fc) in enumerate(faces(el))
    #         on_target(fc) || continue

    #         r = 0
    #         for k in nzrange(D,P[p]) # P ist Liste mit tet-Nummern (hier chronologisch aber nicht immer deshalb so wie hier...)
    #             vals[k] == q && (r = rows[k]; break)
    #         end
    #         @assert r != 0
            
    #         fc1 = chart(ogeo, r)

    #         # Betrachte alle nur von dem tet aus, dessen außernormale der betroffen Fläche in die selbe Richtung zeigt
    #         dotp = dot(fc.normals[1], fc1.normals[1])
    #         dotp ≈ -1.0 && continue  # entweder ≈1 oder ≈-1 (siehe q in nzrange...)
    #         @assert dotp ≈ 1.0
            
    #         Q = ntrace(x, el, q, fc1)

    #         for i in 1:size(Q,1)
    #             for j in 1:size(Q,2)
    #                 for (m,a) in ad[p,j]

    #                     v = a*Q[i,j]

    #                     # Berechne die 1-2 angrenzenden Tetraeder von fc1
    #                     tets = Vector{Int64}()
    #                     for k in nzrange(Dt,r) # r ist der index von fc1
    #                         push!(tets,rowst[k])
    #                     end
    #                     n̂_n = fc1.normals[1]

    #                     if length(tets) == 2
    #                         i_tet1 = tets[1]
    #                         i_tet2 = tets[2]
    #                         tet1 = E[i_tet1]
    #                         tet2 = E[i_tet2]
    #                         center1 = cartesian(CompScienceMeshes.center(tet1))
    #                         center2 = cartesian(CompScienceMeshes.center(tet2))
    #                         centerf = cartesian(CompScienceMeshes.center(fc1))
    #                         v_f1 = center1 - centerf
    #                         v_f2 = center2 - centerf

    #                         if dot(v_f1, n̂_n) > 0.0
    #                             tet_minus = tet1
    #                             tet_plus = tet2
    #                             i_tet_minus = i_tet1 
    #                             i_tet_plus = i_tet2
    #                             @assert dot(v_f2, n̂_n) < 0.0
    #                         elseif dot(v_f1, n̂_n) < 0.0
    #                             tet_minus = tet2
    #                             tet_plus = tet1
    #                             i_tet_minus = i_tet2
    #                             i_tet_plus = i_tet1
    #                             @assert dot(v_f2, n̂_n) > 0.0
    #                         end
    #                         χ_minus = cell2mat_χ[i_tet_minus]
    #                         χ_plus = cell2mat_χ[i_tet_plus]
    #                         δχ = χ_minus - χ_plus

    #                     elseif length(tets) == 1
    #                         i_tet1 = tets[1]
    #                         tet1 = E[i_tet1]
    #                         center1 = cartesian(CompScienceMeshes.center(tet1))
    #                         centerf = cartesian(CompScienceMeshes.center(fc1))
    #                         v_f1 = center1 - centerf
                            
    #                         if dot(v_f1, n̂_n) > 0.0
    #                             tet_minus = tet1
    #                             tet_plus = nothing
    #                             i_tet_minus = i_tet1 
    #                             i_tet_plus = nothing
    #                             χ_minus = cell2mat_χ[i_tet_minus]
    #                             χ_plus = 0.0
    #                         elseif dot(v_f1, n̂_n) < 0.0
    #                             tet_minus = nothing
    #                             tet_plus = tet1
    #                             i_tet_minus = nothing
    #                             i_tet_plus = i_tet1
    #                             χ_minus = 0.0
    #                             χ_plus = cell2mat_χ[i_tet_plus]
    #                         end
    #                         δχ = χ_minus - χ_plus

    #                     else
    #                         error("length(tets) problem!")
    #                     end
    #                     #@show δχ
    #                     v_new = v * δχ
    #                     isapprox(v_new,0,atol=sqrt(eps(T))) && continue # WICHTIG sonst unnötig viele Nullen
    #                     push!(fns[m], BEAST.Shape(r, i, v_new))
    #                 end
    #             end
    #         end

    #     end

    # end

    

cnt = 0
for el in F.fns
    el == [] && (cnt += 1)
end
@show cnt


r = 0
for k in nzrange(D,1000) # P ist Liste mit tet-Nummern (hier chronologisch aber nicht immer deshalb so wie hier...)
    @show k, rows[k], vals[k]
    #vals[k] == 3 && (r = rows[k]; break)
end
@show r
@assert r != 0


Dt = connectivity(ogeo, igeo, abs)
rowst, valst = rowvals(Dt), nonzeros(Dt)
for k in nzrange(Dt,1000) # P ist Liste mit tet-Nummern (hier chronologisch aber nicht immer deshalb so wie hier...)
    @show k, rowst[k], valst[k]
    #vals[k] == 3 && (r = rows[k]; break)
end



# X.geo.faces
# elementsX, adX, cellsX = assemblydata(md.X)
# w.geo.faces
# elementsw, adw, cellsw = assemblydata(md.w)

#face2cell[170]
#plt = Visu.iplot()
#plt = Visu.points(elementsw[170].vertices)
#Visu.simplex(plt, elementsX[1025])







