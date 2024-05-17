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

h = 2.0 # kleiner 0.2 sonst std
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
dirichletnodes = append!(deepcopy(topnodes), bottomnodes)
y_d = lagrangec0d1(Γ, dirichletnodes, Val{3})
# INFO: dirichletnodes[i] gehört zu y_d.pos[i] also "Spannungsvektor": [10V 10V ..... 10V 0V 0V ..... 0V 0V]
#Visu.fnspos(y_d, Visu.mesh(Γ))

# linearlag auf Γ_nc
y = lagrangec0d1(Γ_nc, dirichlet = true) 
#Visu.fnspos(y, Visu.mesh(Γ))

# PWC auf Γ
w = lagrangecxd0(Γ_c)# w=ntrc(X) geht nicht!
#Visu.fnspos(w, Visu.mesh(Γ))

# SWG auf Ω (ohne Γ_nc Flächen)
swgfaces = SWGfaces(Ω, Γ_nc) # Quadratische Komplexität ist extrem langsam!!!! ----> Octree???
X = nedelecd3d(Ω, Mesh(Ω.vertices, swgfaces))#X = nedelecd3d(Ω)
@assert length(X.pos) == length(swgfaces)
ntrc = X -> BEAST.ntrace(X, Γ)
#Visu.fnspos(X, Visu.mesh(Γ))


@show numfunctions(y)
@show numfunctions(y_d)
@show numfunctions(w)
@show numfunctions(X)
## #########################################################


κ = x -> 1.0          # 1.0  Katastrophe
κ0 = 1.0           # -1.0 ...

τ, inv_τ, τ0, χ = gen_tau_chi(problemtype = :current, kappa = κ, kappa0 = κ0)
p = SVector(0.0,0.0,0.0)
τ(p)
inv_τ(p)
τ0
χ(p)
@assert χ(p) - (τ(p)/τ0 - 1)*1/τ(p) <1e-10



# Anregung
u_top = ones(length(topnodes)) * 0.5  # Volle Symmetrie!
u_bottom = ones(length(bottomnodes)) * (-0.5)
ex = append!(deepcopy(u_top), u_bottom)

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


B21_ΓΓ = IPVIE2.B21_ΓΓ(beta = -1.0, gammatype = Float64) # <----- noch konkret Begründen -1...
#norm(assemble(B21_ΓΓ, y, y))
B22_Γ = IPVIE2.B22_Γ(alpha = -1.0, invtau = inv_τ)
#norm(assemble(B22_Γ, y, w))
B22_ΓΓ = IPVIE2.B22_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
#norm(assemble(B22_ΓΓ, y, w))
B23_ΓΓ = IPVIE2.B23_ΓΓ(alpha = 1.0, gammatype = Float64, chi=χ) #VZ? sollte passen
#norm(assemble(B23_ΓΓ, y, ntrc(X)))
B23_ΓΩ = IPVIE2.B23_ΓΩ(alpha = 1.0, gammatype = Float64, chi=χ) #+extra Term?
#norm(assemble(B23_ΓΩ, y, X))


B31_ΓΓ = IPVIE2.B31_ΓΓ(alpha = 1.0, gammatype = Float64)
norm(assemble(B31_ΓΓ, ntrc(X), y))
B31_ΩΓ = IPVIE2.B31_ΩΓ(alpha = -1.0, gammatype = Float64)
norm(assemble(B31_ΩΓ, X, y))
B32_ΓΓ = IPVIE2.B32_ΓΓ(alpha = 1.0, gammatype = Float64, invtau = inv_τ)
norm(assemble(B32_ΓΓ, ntrc(X), w))
B32_ΩΓ = IPVIE2.B32_ΩΓ(alpha = -1.0, gammatype = Float64, invtau = inv_τ)
norm(assemble(B32_ΩΓ, X, w))
B33_Ω = IPVIE2.B33_Ω(alpha = -1.0, invtau = inv_τ)
#assemble(B33_Ω, X, X)
B33_ΓΓ = IPVIE2.B33_ΓΓ(alpha = 1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΓΓ, ntrc(X), ntrc(X))
B33_ΓΩ = IPVIE2.B33_ΓΩ(alpha = -1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΓΩ, ntrc(X), X)
B33_ΩΓ = IPVIE2.B33_ΩΓ(alpha = -1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΩΓ, X, ntrc(X))
B33_ΩΩ = IPVIE2.B33_ΩΩ(alpha = 1.0, gammatype = Float64, chi = χ)
#assemble(B33_ΩΩ, X, X)


# RHS assemble test bzgl Überlapp der Identity
# norm(assemble(B11_Γ, w, y))  
# norm(assemble(B11_Γ, w, y_d))

# norm(assemble(B11_ΓΓ, w, y))
# norm(assemble(B11_ΓΓ, w, y_d))

# norm(assemble(B21_ΓΓ, y, y))
# norm(assemble(B21_ΓΓ, y, y_d))

# norm(assemble(B31_ΓΓ, ntrc(X), y))
# norm(assemble(B31_ΓΓ, ntrc(X), y_d))

# norm(assemble(B31_ΩΓ, X, y))
# norm(assemble(B31_ΩΓ, X, y_d))


# LHS
@hilbertspace i j k # Zeilen    ->Test
@hilbertspace l m n # Spalten   ->Basis

lhs = @varform(
    B11_Γ[i,l] + B11_ΓΓ[i,l] +
    B12_ΓΓ[i,m] +
    B13_ΓΓ[i,ntrc(n)] + B13_ΓΩ[i,n] + 

    B21_ΓΓ[j,l] + 
    B22_Γ[j,m] + B22_ΓΓ[j,m] +
    B23_ΓΓ[j,ntrc(n)] + B23_ΓΩ[j,n] +

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
#@hilbertspace i j k # Zeilen
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


# S*u = R*ex
b = R*ex
u = S \ b
#@assert norm(S*u - b) < 1e-8
u_Φ = u[1:length(y)]
u_Jn = u[length(y)+1:length(y)+length(w)]
u_J = u[length(y)+length(w)+1:end]
@assert length(u_Φ) == length(y.fns)
@assert length(u_Jn) == length(w.fns)
@assert length(u_J) == length(X.fns)


# Ränder
@show maximum(u_Φ)
@show minimum(u_Φ)

@show maximum(u_Jn)
@show minimum(u_Jn)

# Stomdichte
range_ = range(-0.49,stop=0.49,length=9)
points = [point(x,y,z) for x in range_ for y in range_ for z in range_]
J_MoM = BEAST.grideval(points, u_J, X)#, type=Float64)

display("Volume Current Density - Total Volume:")
Jallx, Jally, Jallz = pointlist2xyzlist(J_MoM)
@show sum(Jallz)/length(Jallz)
@show sum(Jallx)/length(Jallx)
@show sum(Jally)/length(Jally)

display(Visu.fieldplot(points, J_MoM, 0.25, Visu.mesh(Γ_c)))


# Stromdichte in Ebene z0
range_xy = range(-0.5,stop=0.5,length=9)
z0 = -0.4999999
points2 = [point(x,y,z0) for x in range_xy for y in range_xy]
J_MoM2 = BEAST.grideval(points2, u_J, X)

display("Volume Current Density - z0 Plane:")
Jallx, Jally, Jallz = pointlist2xyzlist(J_MoM2)
@show sum(Jallz)/length(Jallz)  
@show sum(Jallx)/length(Jallx)
@show sum(Jally)/length(Jally)
#display(Visu.fieldplot(points2, J_MoM2, 1.0, Visu.mesh(Γ_c)))
display("")

I_top, I_bottom = getcurrent(u_Jn, w, Γ_c_t, Γ_c_b)

@show I_top
@show I_bottom

@show I_top*τ0
@show I_bottom*τ0


# Potential: Randknoten vs. Analytisch

function CapacitorPot(x, u_top, u_bottom; scale = 1.0) # vereinfacht!!!
    lz = 1.0
    m = scale * (u_top[1]-u_bottom[1]) / lz

    z = x[3]
    return m * z 
end
u_Φana = Vector{Float64}(undef,length(u_Φ))
for (i,pos) in enumerate(y.pos)
    u_Φana[i] = CapacitorPot(pos, u_top, u_bottom, scale = 1.0)
end
@show norm(u_Φ-u_Φana)/norm(u_Φana) #nicht punktweise

## facecurrents Tests

# J_n auf Γ_c
fcr0, geo0 = facecurrents(u_Jn, w)
Plotly.plot(patch(geo0, fcr0))      # FALSCH ORIENTIERT!!!

# Φ auf Γ_nc -> Achtung an Plattengrenzen fehlt noch Dirichlet Beitrag!
fcr1, geo1 = facecurrents(u_Φ, y)
Plotly.plot(patch(geo1, fcr1))      #MANCHMAL FALSCH ORIENTIERT!!! je nach tau0+-   => vmtl doch irgendwie * τ0

# Φ auf Γ_c    
fcr2, geo2 = facecurrents(ex, y_d)
Plotly.plot(patch(geo2, fcr2))          






""
##
t






















# ERGEBNISSE bei τ0=1.0 : 

# J_z - Komponente ist in der Nähe der Kontaktflächen unphysikalisch! 
#      d.h. +/- Sprünge benachbarter Zellen...

# Doppelte Leitfähigkeit => noch nicht doppelter Strom!!!
# Doppelte Spannung => noch nicht doppelter Strom!!!
# Halbe Leitfähigkeit => Halber Strom!!! STIMMT
# irgendwas passt bei tau0 noch nicht...

""
##
# Richtungstest J an Übergängen der swgfaces
for face in swgfaces
    patch = simplex(Ω.vertices[face]) 
    n = patch.normals[1]
    c = CompScienceMeshes.center(patch) # mp
    p = cartesian(c)
    p1 = p + n*1e-7
    p2 = p - n*1e-7
    J_list = BEAST.grideval([p1,p2], u[length(y)+1:end], X)
    
    J1 = J_list[1]
    J2 = J_list[2]
    J1n = dot(J1,n)
    J2n = dot(J2,n)

    if (J1n != 0.0) && (J2n != 0.0)
        @assert sign(J1n) == sign(J1n)
    end
end


##

# ZUSAMMENHANG charts - DOF
charts, ad, ag = assemblydata(y)
length(charts)
r=1
chart_tree = BEAST.octree(charts)
i1 = CompScienceMeshes.findchart(charts, chart_tree, p)
i1



for (m,we) in ad[37, r]
    @show m, we
    @show w.pos[m]
    #values[j] += we * coeffs[m] * vals[r][1]
end
cartesian(CompScienceMeshes.center(charts[37]))


## 


    
       

# npos=0
# nneg=0
# for el in u_Jn
#     sign(el) < 0.0 && (nneg += 1)
#     sign(el) > 0.0 && (npos += 1) 
# end
# npos
# nneg




# S1 = BEAST._spacedict_to_directproductspace(lhs.test_space)
# S2 = BEAST._spacedict_to_directproductspace(lhs.trial_space)
# assemble(lhs,S1,S2)

# vie = @discretise lhs k∈y, l∈X, j∈y, m∈X, n∈y_d


#     ==
#     (-TL_Γ-TL_ΓΓ)[k,n] +                 #??? Wie bekommt man da Matrix*gegebenen Vektor hinein??? 
#     (-BL_ΓΓ)[ntrc(l),n] + (-BL_ΩΓ)[l,n], #??? Wie bekommt man da Matrix*gegebenen Vektor hinein??? 
#     k∈y, l∈X, j∈y, m∈X, n∈y_d
# )
# u = solve(vie)
#@assert lengths...


#vie = @discretise( )




#a1, a2, a3 = assemblydata(X)



##
# t1=[[u,v,w] for w in 0.0:0.09:1 for v in 0:0.09:1-w for u in 0:0.09:1-v-w]

# range_ = collect(range(start=0.0, stop=1.0, length=3))
# bary_list = [[u, v, w] for u in range_ for v in range_-u for w in range_-u-v]

##

s1 = simplex(Ω.vertices[Ω.faces[3]])

plt = Visu.mesh(Γ_c_t)
plt = Visu.simplex(plt, s1)
plt = Visu.add1(plt, s1, refspace(X))


##
plt=Visu.mesh(Ω)

Visu.add1(plt, s1, refspace(X))


##



@show length(skeleton(Ω,2))
@show numfunctions(X)
#y = boundary(Ω)

Visu.mesh(Ω)

r=refspace(X)
p=point(0.0,0.1,0.3)
#r(p)
##

barycoord1=point(0.0,0.1,0.3)
@assert norm(barycoord1) <= 1.0 
s1 = simplex(Ω.vertices[Ω.faces[3]])
n1 = neighborhood(s1,barycoord1)

# jetzt kann man r(CSM-MP) anwenden
r(n1)[4]






# α, α′ = 1/η, 1/η′
# vie = @discretise(
#     (η*T+η′*T′)[k,j] -      (K+K′)[k,m] +
#          (K+K′)[l,j] + (α*T+α′*T′)[l,m] == -e[k] - h[l],
#     j∈X, m∈X, k∈X, l∈X)

# u = solve(vie)
#@assert lengths...









##


# plotly()
# plt = plot(legend=false, xlims=(-1,1), ylims=(-1,1), zlims=(-1,1), size=(850,850))

# # Zeichnen Sie einen Pfeil von Punkt a nach Punkt b
# pnt = [0.1,-0.4,0.3]
# dir = [-0.3,-0.3,0.8]

# Visu.draw_arrow!(plt, pnt, dir; scale = 0.1, arrcolor = "blue", arrwidth = 3)

# # Zeigen Sie den Plot an
# range_ = range(start=0.0, stop=1.0, length=3)




##






default(fmt=:svg)
default(fmt=:png)

###############################################################################
# ADD THIS BLOCK TO REMOVE SVG FROM LIST OF "DISPLAYABLE_MIMES":
# pos = findfirst((x)->(x=="image/svg+xml"), VSCodeServer.DISPLAYABLE_MIMES)
# if !isnothing(pos)
#     popat!(VSCodeServer.DISPLAYABLE_MIMES, pos)
#     println("Popped!")
# end
###############################################################################





#test_elements, tad, a1= assemblydata(y_d)
dirichlet_n2f(y_d, dirichletnodes)

y_d.geo.vertices[dirichletnodes]
y_d.geo.vertices
#test_elements, tad, a1= assemblydata(X)
# y_d 116 fns => die (n,val) Einträge der tad.data Matrix haben n=(0),1,2,...,116 - 116 ist die Zahl der dirichletnodes
# also 284 Dreiecke
dirichletnodes[100]


tad.data[1,1,1]

length(y_d.fns)
length(realvertices(Γ))
length(Γ_c.faces)
length(Γ_nc.faces)
#length(tad.data)

max=0
min=1e30
for el in tad.data #! ist matrix!!!!
    @show el
    max <= el[1] && (max = el[1])
    min >= el[1] && (min = el[1])
end
max
min


##



#lhs = eq.equation.lhs
#rhs = eq.equation.rhs

@hilbertspace a1
@hilbertspace a2

eq1 = @varform BR_ΩΓ[a1,ntrc(a2)] + BR_ΓΓ[ntrc(a1),ntrc(a2)] + (BR_Ω+BR_ΩΩ)[a1,a2]
vie1 = @discretise eq1 a1∈X a2∈X 
tedic=vie1.test_space_dict
trdic=vie1.trial_space_dict
S1 = BEAST._spacedict_to_directproductspace(tedic)
S2 = BEAST._spacedict_to_directproductspace(trdic)
assemble(eq1, S1, S2)




# # LHS Assemble Test 
# assemble(TL_Γ, y, y) # ist 1/2 I
# assemble(TL_ΓΓ, y, y)

# assemble(TR_ΓΓ, y, ntrc(X))
# assemble(TR_ΓΩ, y, X)

# assemble(BL_ΓΓ, ntrc(X), y)
# assemble(BL_ΩΓ, X, y)
# #
# assemble(BR_Ω, X, X)
# assemble(BR_ΓΓ, ntrc(X), ntrc(X))
# assemble(BR_ΩΓ, X, ntrc(X))
# assemble(BR_ΓΩ, ntrc(X), X)
# assemble(BR_ΩΩ, X, X)

# # RHS Assemble Test, später noch * geg dirichlet spannungsvektor
# assemble(-TL_Γ, y, y_d)
# assemble(-TL_ΓΓ, y, y_d)

# assemble(-BL_ΓΓ, ntrc(X), y_d)
# assemble(-BL_ΩΓ, X, y_d)