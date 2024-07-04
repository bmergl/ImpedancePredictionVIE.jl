using MKL

using ImpedancePredictionVIE
using CompScienceMeshes
using Plots
using BEAST
using LinearAlgebra
using StaticArrays
using Plotly





geoname = "cube.geo"
geopath = "$(pkgdir(ImpedancePredictionVIE))/geo/$geoname"

meshname = "cube.msh"
meshpath = "$(pkgdir(ImpedancePredictionVIE))/geo/$meshname"

h = 0.18 # kleiner 0.18 sonst std   0.18 -> 0.09 -> 0.045 für Konvergenztest
Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc = geo2mesh(geopath, meshpath, h)

# Visu.mesh(Ω)
# Visu.mesh(Γ_c)
# Visu.mesh(Γ_c_t)
# Visu.mesh(Γ_c_b)
# Visu.mesh(Γ_nc)
# Visu.mesh(Γ)


# linearlag (dirichlet) auf Γ_c
topnodes = realnodes(Γ_c_t) # wichtig denn später z.B. 10V
bottomnodes = realnodes(Γ_c_b) # ""  z.B. 0V
dirichletnodes = vcat(topnodes, bottomnodes)
y_d = lagrangec0d1(Γ, dirichletnodes, Val{3})
# INFO: dirichletnodes[i] gehört zu y_d.pos[i] also "Spannungsvektor": [10V 10V ..... 10V 0V 0V ..... 0V 0V]
#Visu.fnspos(y_d, Visu.mesh(Γ))

# linearlag auf Γ_nc
y = lagrangec0d1(Γ_nc, dirichlet = true) 
#Visu.fnspos(y, Visu.mesh(Γ))

# SWG auf Ω (ohne Γ_nc Flächen)
swg_faces = swgfaces(Ω, Γ_nc, fast = true)
#swg_faces2 = swgfaces(Ω, Γ_nc, fast = false)
#norm(norm.(swg_faces)-norm.(swg_faces2))
X = nedelecd3d(Ω, Mesh(Ω.vertices, swg_faces))
@assert length(X.pos) == length(swg_faces)
#Visu.fnspos(X, Visu.mesh(Γ))
ntrc = X -> BEAST.ntrace(X, Γ)
ntrcX = ntrc(X)
#Visu.fnspos(ntrc(X), Visu.mesh(Γ))

# PWC auf Γ
w_ = lagrangecxd0(Γ_c) # War doch alles egal... hat keinen Einfluss auf die Lösung  ;(
# w_fns = Vector{Vector{BEAST.Shape{Float64}}}(undef, length(w_.fns))
# w_pos = Vector{SVector{3, Float64}}(undef, length(w_.fns))
# w_charts = assemblydata(w_)[1]
# w_charts_tree = BEAST.octree(w_charts)
# for (i,fns) in enumerate(ntrcX.fns)
#     if fns != []
#         pos = ntrcX.pos[i]

#         j = CompScienceMeshes.findchart(w_charts, w_charts_tree, pos) # orientieren an w_
#         #s = sign(fns[1].coeff)
#         s = -sign(randn())

#         w_fns[j] = [BEAST.Shape(j, fns[1].refid, s)]    # evtl. doch nur VZ weitergeben....
#         w_pos[j] = pos
#     end
# end
# w_.pos = w_pos
# w_.fns = w_fns 
w = w_
#Visu.fnspos(w, Visu.mesh(Γ))



@show numfunctions(y)
@show numfunctions(y_d)
@show numfunctions(w)
@show numfunctions(X)
## ########################################################

κ = x -> 2.0
κ0 = 1.0
# function genkappa()
#     function kappa(x)
#         # @show x
#         # rand() <0.1  && error("STOP")
#         x[3] < -0.2 && return 100.0 
#         # return 10.0
        
#         #sqrt(x[1]^2+x[2]^2) < 0.2 && return 1000.0
#         return 10.0
#     end
#     return kappa
# end
# κ = genkappa()

τ, inv_τ, τ0, χ = gen_tau_chi(problemtype = :current, kappa = κ, kappa0 = κ0)
p = SVector(0.0,0.0,0.0)
τ(p)
inv_τ(p)
τ0
χ(p)
#@assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) <1e-10


BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = BEAST.SingleNumQStrat(3)

#BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(4,4,8,6,6,6)
#BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = BEAST.SauterSchwab3DQStrat(3,3,3,3,3,3)

#BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,5,5,5,5)
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = BEAST.DoubleNumWiltonSauterQStrat(3,3,3,3,3,3,3,3)
#BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = BEAST.DoubleNumWiltonSauterQStrat(9,9,9,9,9,9,9,9)

# Anregung
u_top = ones(length(topnodes)) * 0.5  # Volle Symmetrie!
u_bottom = ones(length(bottomnodes)) * (-0.5)
v = vcat(u_top, u_bottom)

# Definiton der Operatoren

B11_Γ = IPVIE2.B11_Γ(alpha = 1.0) #Verschwindet!!!
#assemble(B11_Γ, w, y)
B11_ΓΓ = IPVIE2.B11_ΓΓ(alpha = 1.0, gammatype = Float64)
#norm(assemble(B11_ΓΓ, w, y))
B12_ΓΓ = IPVIE2.B12_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
#norm(assemble(B12_ΓΓ, w, w))
B13_ΓΓ = IPVIE2.B13_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
#norm(assemble(B13_ΓΓ, w, ntrc(X)))
B13_ΓΩ = IPVIE2.B13_ΓΩ(alpha = -1.0, gammatype = Float64, chi = χ)
#norm(assemble(B13_ΓΩ, w, X))


B21_ΓΓ = IPVIE2.B21_ΓΓ(beta = -1.0, gammatype = Float64) # -1.0 siehe BEAST & Steinbach 6.5
#norm(assemble(B21_ΓΓ, y, y))
B22_Γ = IPVIE2.B22_Γ(alpha = -1.0, invtau = inv_τ)
#norm(assemble(B22_Γ, y, w))
B22_ΓΓ = IPVIE2.B22_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
#norm(assemble(B22_ΓΓ, y, w))
B23_ΓΓ = IPVIE2.B23_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ) #VZ? sollte passen
#norm(assemble(B23_ΓΓ, y, ntrc(X)))
B23_ΓΩ = IPVIE2.B23_ΓΩ(alpha = 1.0, gammatype = Float64, chi = χ)
#norm(assemble(B23_ΓΩ, y, X))

B23_dyad = IPVIE2.B23_dyad(alpha = 1.0, gammatype = Float64, chi = χ) # !!! Hypersingular - keine Konvergenz

B23_constmed = IPVIE2.B23_constmed(alpha = 1.0, gammatype = Float64, chi = χ)


# # KONTROLLE DER MATRIXELEMENTE ####################################################

# S1 = assemble(B23_ΓΓ, y, ntrc(X))
# S2 = assemble(B23_ΓΩ, y, X) # cross Formulierung
# S0 = S1+S2

# S3 = assemble(B23_constmed, y, X)

# @show norm(S0-S3) # bei ! Alle KErnel anpassen


# ## 

# #@show maximum(abs.(S0-S3)./abs.(S3))
# relErrM = abs.(S0-S3)./abs.(S3)

# ind2Mind = CartesianIndices(relErrM)
# pair_list = []
# cnt = 0
# for (index, element) in enumerate(relErrM) 
#     i = ind2Mind[index][1]
#     j = ind2Mind[index][2]
#     s0 = S0[i,j]
#     s3 = S3[i,j]
    
#     lim = 0.02

#     if element > lim # 50.0 && element > 0.3 #&& abs(S1[i,j]) == 0.0 #
#         if norm(y.pos[i]-X.pos[j]) > -1.0 && s3 > maximum(S3)/1000
#         println("relErr: $element, S0[$i,$j]=$s0, S3[$i,$j]=$s3")
#         push!(pair_list,[i,j])
#         cnt += 1
#         end
#     end

# end
# @show cnt/(size(S0,1)*size(S0,2))

# display(pair_list)

# ##

# pair = pair_list[1]
# i = pair[1]
# j = pair[2]

# # i = 61
# # j = 8

# S1[i,j]
# S2[i,j]

# S3[i,j]

# (abs(S2[i,j])-S3[i,j])/S3[i,j]


# c1 = X.fns[j][1].cellid
# c2 = X.fns[j][2].cellid
# vertc1 = Ω.vertices[Ω.faces[c1]]
# vertc2 = Ω.vertices[Ω.faces[c2]]

# plt = Visu.mesh(Γ_nc)
# Visu.simplex(plt,vertc1)
# Visu.simplex(plt,vertc2)
# #Visu.points([y.pos[i],],plt)
# Visu.points([y.pos[i],X.pos[j]],plt)

#

# y.fns[i]
# y.pos[i]

# c_1 = y.fns[i][1].cellid
# c_2 = y.fns[i][2].cellid
# c_3 = y.fns[i][3].cellid
# c_4 = y.fns[i][4].cellid
# c_5 = y.fns[i][5].cellid
# c_6 = y.fns[i][6].cellid

# #Visu.points(vertc_5,plt)
# Γ.faces[c_1]
# Γ.vertices[16]

# ##
# i = 76
# j = 1883

# S1[i,j]
# S2[i,j]
# S0[i,j]

# S3[i,j]
# ~,sv0,~ = svd(S0)
# p1 = Plots.plot(sv0)
# ~,sv3,~ = svd(S3)
# plot!(p1,sv3)

#


B31_ΓΓ = IPVIE2.B31_ΓΓ(alpha = 1.0, gammatype = Float64, dyad = false)
#norm(assemble(B31_ΓΓ, ntrc(X), y))
B31_ΩΓ = IPVIE2.B31_ΩΓ(alpha = -1.0, gammatype = Float64, dyad = false)
#norm(assemble(B31_ΩΓ, X, y))
B32_ΓΓ = IPVIE2.B32_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ, dyad = false)
#norm(assemble(B32_ΓΓ, ntrc(X), w))
B32_ΩΓ = IPVIE2.B32_ΩΓ(alpha = -1.0, gammatype = Float64, invtau = inv_τ, dyad = false)
#norm(assemble(B32_ΩΓ, X, w))
B33_Ω = IPVIE2.B33_Ω(alpha = -1.0, invtau = inv_τ)
#norm(assemble(B33_Ω, X, X))
B33_ΓΓ = IPVIE2.B33_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΓΓ, ntrc(X), ntrc(X))
B33_ΓΩ = IPVIE2.B33_ΓΩ(alpha = -1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΓΩ, ntrc(X), X)
B33_ΩΓ = IPVIE2.B33_ΩΓ(alpha = -1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΩΓ, X, ntrc(X))
B33_ΩΩ = IPVIE2.B33_ΩΩ(alpha = 1.0, gammatype = Float64, chi = χ)
#norm(assemble(B33_ΩΩ, X, X))


# LHS
@hilbertspace i j k # Zeilen    ->Test
@hilbertspace l m n # Spalten   ->Basis
lhs = @varform(
    B11_Γ[i,l] + B11_ΓΓ[i,l] +
    B12_ΓΓ[i,m] +
    B13_ΓΓ[i,ntrc(n)] + B13_ΓΩ[i,n] + 

    B21_ΓΓ[j,l] + 
    B22_Γ[j,m] + B22_ΓΓ[j,m] +
    B23_ΓΓ[j,ntrc(n)] + 
    B23_ΓΩ[j,n] +
    #B23_dyad[j,n] +
    #B23_constmed[j,n] +
    

    B31_ΓΓ[ntrc(k),l] + B31_ΩΓ[k,l] +
    B32_ΓΓ[ntrc(k),m] + B32_ΩΓ[k,m] +
    B33_Ω[k,n] + B33_ΓΓ[ntrc(k),ntrc(n)] + B33_ΓΩ[ntrc(k),n] + B33_ΩΓ[k,ntrc(n)] + B33_ΩΩ[k,n]
)
lhsd = @discretise lhs i∈w j∈y k∈X l∈y m∈w n∈X #! w und y Tausch!
lhsd_test = lhsd.test_space_dict
lhsd_trial = lhsd.trial_space_dict
testSpace_lhs = BEAST._spacedict_to_directproductspace(lhsd_test)
trialSpace_lhs = BEAST._spacedict_to_directproductspace(lhsd_trial)
M = assemble(lhs, testSpace_lhs, trialSpace_lhs)
S = Matrix(M)

# RHS
@hilbertspace o # Spalten, !!! Nur eine Blockspalte
rhs = @varform( # Vorlage für nicht-quadratische Matrix ...
    -B11_Γ[i,o] -B11_ΓΓ[i,o] +

    -B21_ΓΓ[j,o] + 

    -B31_ΓΓ[ntrc(k),o] -B31_ΩΓ[k,o]
)
rhsd = @discretise rhs i∈w j∈y k∈X o∈y_d
rhsd_test = rhsd.test_space_dict
rhsd_trial = rhsd.trial_space_dict
testSpace_rhs = BEAST._spacedict_to_directproductspace(rhsd_test)
trialSpace_rhs = BEAST._spacedict_to_directproductspace(rhsd_trial)
R = Matrix(assemble(rhs, testSpace_rhs, trialSpace_rhs))

# S*u = R*v
b = R*v
u = S \ b
#@assert norm(S*u - b) < 1e-8
u_Φ = u[1:length(y)]
u_Jn = u[length(y)+1:length(y)+length(w)]
u_J = u[length(y)+length(w)+1:end]
@assert length(u_Φ) == length(y.fns)
@assert length(u_Jn) == length(w.fns)
@assert length(u_J) == length(X.fns)

#@show cond(S)
## selbst höchste genauigkeit bringt nichts 7% 7% 13% wird nicht unterschritten
# ... suche also für χ ungleich null nach χ-Fehlerterm - ist es eine Dyade???


# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, u_J, X)#, type=Float64)
display("Stomdichte Gesamtvolumen")
Jx, Jy, Jz = pointlist2xyzlist(J_MoM)
Jana = -ones(length(Jz))*κ(p)*(u_top[1]-u_bottom[1])
@show norm(Jz - Jana)/norm(Jana)

# Stromdichte Mitte: Ebene z=0.0
range_xy = range(-0.5,stop=0.5,length=9)
points2 = [point(x,y,0.0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, u_J, X)
display("Stromdichte Mitte: Ebene z=0.0")
Jx2, Jy2, Jz2 = pointlist2xyzlist(J_MoM2)
Jana2 = -ones(length(Jz2))*κ(p)*(u_top[1]-u_bottom[1])
@show norm(Jz2 - Jana2)/norm(Jana2)
#display(Visu.fieldplot(points2, J_MoM2, 1.0, Visu.mesh(Γ_c)))

# Stromdichte bei Platten: Ebene z=0.49
points3 = [point(x,y,0.49) for x in range_xy for y in range_xy]
J_MoM3 = BEAST.grideval(points3, u_J, X)
display("Stromdichte bei Platten: Ebene z=0.49")
Jx3, Jy3, Jz3 = pointlist2xyzlist(J_MoM3)
Jana3 = -ones(length(Jz3))*κ(p)*(u_top[1]-u_bottom[1])
@show norm(Jz3 - Jana3)/norm(Jana3)
#display(Visu.fieldplot(points2, J_MoM2, 1.0, Visu.mesh(Γ_c)))

display(Visu.fieldplot(points, J_MoM, 0.06, Visu.mesh(Γ_c)))


# Strom durch Platten
display("")
I_top, I_bottom = getcurrent(u_Jn, w, Γ_c_t, Γ_c_b)
@show I_top
@show I_bottom
Iana = 1.0 * κ(p) *(u_top[1]-u_bottom[1])
@show norm(I_top-Iana)/norm(Iana)
@show norm(I_bottom-Iana)/norm(Iana)
display("")

# Potential: Randknoten vs. Analytisch

function CapacitorPot(x, u_top, u_bottom; scale = 1.0) # vereinfacht für hom. medium
    lz = 1.0
    m = scale * (u_top[1]-u_bottom[1]) / lz

    z = x[3]
    return m * (z + lz/2) + u_bottom[1]
end
u_Φana = Vector{Float64}(undef,length(u_Φ))
for (i,pos) in enumerate(y.pos)
    u_Φana[i] = CapacitorPot(pos, u_top, u_bottom, scale = 1)
end
@show norm(u_Φ-u_Φana)/norm(u_Φana) #nicht punktweise



## facecurrents Tests

# J_n auf Γ_c
fcr0, geo0 = facecurrents(u_Jn, w)
Plotly.plot(patch(geo0, fcr0))

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(u_Φ, y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

# Φ auf Γ_c    
fcr2, geo2 = facecurrents(ex, y_d)
Plotly.plot(patch(geo2, fcr2))          


