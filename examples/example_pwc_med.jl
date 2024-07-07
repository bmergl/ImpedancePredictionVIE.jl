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


md = IP.setup(geoname = "cube.geo", meshname = "cube.msh", body = IP.cuboid(1.0, 1.0, 1.0), h = 0.18)
print("tehrahedrons: ", length(md.Ω.faces))
#Visu.mesh(md.Γ) 

y = md.y
w = md.w
X = md.X

ntrc = md.ntrc

##

# Operators row 1
B11_Γ = IPVIE1.B11_Γ()
assemble(B11_Γ, w, y)
B11_ΓΓ = IPVIE1.B11_ΓΓ(gammatype = Float64)
assemble(B11_ΓΓ, w, y)

UB12_ΓΓ = IPVIE1.UB12_ΓΓ(gammatype = Float64)
assemble(UB12_ΓΓ, w, w)

UB13_ΓΓn = IPVIE1.UB13_ΓΓn(gammatype = Float64)
#assemble(UB13_ΓΓn, w, NEW SPACE NEEDED!!!!)
UB13_ΓΩ = IPVIE1.UB13_ΓΩ(gammatype = Float64)
assemble(UB13_ΓΩ, w, X)


# Operators row 2
B21_ΓΓ = IPVIE1.B21_ΓΓ(gammatype = Float64)
assemble(B21_ΓΓ, y, y)

UB22_Γ = IPVIE1.UB22_Γ()
assemble(UB22_Γ, y, w)
UB22_ΓΓ = IPVIE1.UB22_ΓΓ(gammatype = Float64)
assemble(UB22_ΓΓ, y, w)

UB23_ΓΓn = IPVIE1.UB23_ΓΓn(gammatype = Float64)
#assemble(B23_ΓΓ, y, NEW SPACE NEEDED!!!!)
UB23_ΓΩ = IPVIE1.UB23_ΓΩ(gammatype = Float64)
assemble(UB23_ΓΩ, y, X)


# Operators row 3
B31_ΓΓ = IPVIE1.B31_ΓΓ(gammatype = Float64)
assemble(B31_ΓΓ, ntrc(X), y)
B31_ΩΓ = IPVIE1.B31_ΩΓ(gammatype = Float64)
assemble(B31_ΩΓ, X, y)

UB32_ΓΓ = IPVIE1.UB32_ΓΓ(gammatype = Float64)
assemble(UB32_ΓΓ, ntrc(X), w)
UB32_ΩΓ = IPVIE1.UB32_ΩΓ(gammatype = Float64)
assemble(UB32_ΩΓ, X, w)

UB33_Ω = IPVIE1.UB33_Ω()
assemble(UB33_Ω, X, X)
UB33_ΓΓn = IPVIE1.UB33_ΓΓn(gammatype = Float64)
#assemble(UB33_ΓΓ, ntrc(X), NEW SPACE NEEDED!!!!)
UB33_ΓΩ = IPVIE1.UB33_ΓΩ(gammatype = Float64)
assemble(UB33_ΓΩ, ntrc(X), X)
UB33_ΩΓn = IPVIE1.UB33_ΩΓn(gammatype = Float64)
#assemble(B33_ΩΓ, X, NEW SPACE NEEDED!!!!)
UB33_ΩΩ = IPVIE1.UB33_ΩΩ(gammatype = Float64)
assemble(UB33_ΩΩ, X, X)



# Schritt 1: ntrace artige function mit 
# input: Materialvektor(tetindex) = [30.0 131.6 123.0 ... 98.3]
# output 

X.geo.faces
elementsX, adX, cellsX = assemblydata(md.X)

elements
adX
for (n,b) in adX[860][4]
    @show n, b
end
cells
md.X.fns

# Schritt 2: erstelle X_mat basis für trial in den speziellen Fällen!




#center of a simplex:
el = elements[1135]
center = cartesian(CompScienceMeshes.center(el))
##



testmat = IP.constant_xsplit(100.0, nothing, 0.1, 50.0, nothing)
κ, ϵ = testmat() # hier noch Funktionen
κ0 = 1.0
ϵ0 = nothing
ω = nothing

τ, inv_τ, τ0, χ = gen_tau_chi(kappa = κ, kappa0 = κ0, epsilon = ϵ, epsilon0 = ϵ0, omega = ω) # hier noch Funktionen
cell2mat_inv_τ = Vector{typeof(inv_τ(SVector(0,0,0)))}(undef,length(elements))
cell2mat_χ  = Vector{typeof(χ(SVector(0,0,0)))}(undef,length(elements)) # Typenunterscheidung???
for (n,el) in enumerate(elements)
    center = cartesian(CompScienceMeshes.center(el))
    inv_τ_n = inv_τ(center)
    χ_n = χ(center)
    cell2mat_inv_τ[n] = inv_τ_n 
    cell2mat_χ[n] = χ_n   
end
cell2mat_inv_τ
cell2mat_χ


# erstelle X_mat aus X mittels cell2mat_χ: evtl sollte gen_tau_chi() den Typ T liefern
T = typeof(χ(SVector(0, 0, 0)))
newfns = Vector{Vector{BEAST.Shape{T}}}() 
for (i,shs) in enumerate(X.fns)
    newshs = Vector{BEAST.Shape{T}}()
    for (j,sh) in enumerate(shs)
        cellid = sh.cellid
        refid = sh.refid
        coeff = sh.coeff

        χ_cell = cell2mat_χ[cellid]
        new_coeff = coeff*χ_cell # can be ComplexF64
        push!(newshs, BEAST.Shape(cellid, refid, new_coeff))
    end
    push!(newfns, newshs)
end
X_mat = BEAST.NDLCDBasis(X.geo, newfns, X.pos)


# erstelle w_mat aus w mittels tri->tet, cell2mat_inv_τ:
D = connectivity(w.geo, X.geo)
@assert sum(D) == length(w.fns) == length(w.geo.faces)
rows, vals = rowvals(D), nonzeros(D)
# Bsp rows[i]=j vals[i]=a =>D[j,i]=a
face2cell = rows

T = typeof(inv_τ(SVector(0, 0, 0)))
newfns = Vector{Vector{BEAST.Shape{T}}}() 
for (i,shs) in enumerate(w.fns)
    newshs = Vector{BEAST.Shape{T}}()
    for (j,sh) in enumerate(shs)
        cellid = sh.cellid
        refid = sh.refid
        coeff = sh.coeff

        tricell = cellid
        tetcell = face2cell[tricell]

        inv_τ_tri = cell2mat_inv_τ[tetcell]

        new_coeff = coeff*inv_τ_tri # can be ComplexF64
        push!(newshs, BEAST.Shape(cellid, refid, new_coeff))
    end
    push!(newfns, newshs)
end
w_mat = BEAST.LagrangeBasis{0,-1,1}(w.geo, newfns, w.pos)

# erstelle  Xf_mat aus X und ...
swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)

D = connectivity(swg_faces_mesh, X.geo) # includes boundary tets
rows, vals = rowvals(D), nonzeros(D)
rows
vals

nzrange(D,2180)


r = 0
for k in BEAST.nzrange(D,2180)
    vals[k] == -1 && (r = rows[k]; break)
end
r
@assert r != 0



D
BEAST.nzrange(D,2180)

D[1063,2]

D

D[1063,3]
D[10]


sum(D)
sum(abs.(D))
@assert sum(D) == length(w.fns) == length(w.geo.faces)

# Bsp rows[i]=j vals[i]=a =>D[j,i]=a
face2cell = rows







Xf_mat = BEAST.LagrangeBasis{0,-1,1}(w.geo, newfns, w.pos) # kreuz und quer verteilte Dreiecke


# X.geo.faces
# elementsX, adX, cellsX = assemblydata(md.X)
# w.geo.faces
# elementsw, adw, cellsw = assemblydata(md.w)

#face2cell[170]
#plt = Visu.iplot()
#plt = Visu.points(elementsw[170].vertices)
#Visu.simplex(plt, elementsX[1025])



cnt = 0
for el in md.ntrcX.fns
    el == [] && (cnt += 1)
end
@show cnt



X.fns

t1 = BEAST.Shape(13,3,-1.0+13.3im)

#X.geo.vertices = realvertices(md.Γ)


X.geo.vertices

realvertices(md.Γ)

M=Matrix(connectivity(md.Γ_nc, md.Ω))

md.swg_faces


M=Matrix(connectivity(md.Γ_nc, md.Ω))

swg_faces_mesh = Mesh(md.Ω.vertices, md.swg_faces)

M = Matrix(connectivity(swg_faces_mesh, md.Ω))

maximum(M)
minimum(M)

(sum(abs.(M)) + 1*180 )/2

Visu.mesh(swg_faces_mesh)

asd




X = md.X

γ = swg_faces_mesh

#function ntrace(X::Space, swg_faces_mesh::Mesh)

    γ = swg_faces_mesh

    x = refspace(X)
    E, ad, P = assemblydata(X)
    igeo = geometry(X)
    @assert dimension(γ) == dimension(igeo)-1

    #ogeo = boundary(igeo) <--- ntrace klassisch
    ogeo = γ # ja...
    on_target = overlap_gpredicate(γ)
    ogeo = submesh(ogeo) do m,f
        ch = chart(m,f)
        on_target(ch)
    end

    D = connectivity(igeo, ogeo, abs) # nzrange(D,tetnummer) liefert die 4 bzw. seltener 3 Dreiecke
    rows, vals = rowvals(D), nonzeros(D)
    Dt = connectivity(ogeo, igeo, abs) # nzrange(Dt,trinummer) liefert die 2 bzw. seltener 1 Tetraeder
    rowst, valst = rowvals(Dt), nonzeros(Dt)


    T = typeof(χ(SVector(0, 0, 0)))
    S = BEAST.Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

   
    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))
            on_target(fc) || continue
            on_target(fc) == false && @warn "In this implementation every face should be on the target!"
            
            r = 0
            for k in nzrange(D,P[p]) # P ist Liste mit tet-Nummern (hier chronologisch aber nicht immer deshalb so wie hier...)
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0
            
            fc1 = chart(ogeo, r)
            Q = ntrace(x, el, q, fc1)

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]

                        v = a*Q[i,j]

                        # Berechne die 1-2 angrenzenden Tetraeder von fc1
                        tets = []
                        for k in nzrange(Dt,r) # r ist der index von fc1
                            push!(rowst[k],tets)
                        end
                        n̂_n = fc1.normals[1]
                        @show n̂_n

                        if length(tets) == 2
                            tet1 = tets[1]
                            tet2 = tets[2]
                            center1 = cartesian(CompScienceMeshes.center(tet1))
                            center2 = cartesian(CompScienceMeshes.center(tet2))
                            centerf = cartesian(CompScienceMeshes.center(fc1))
                            v_f1 = center1 - centerf
                            v_f2 = center2 - centerf

                            if dot(v_f1, n̂_n) > 0.0
                                tet_plus = tet1
                            elseif dot(v_f1, n̂_n) < 0.0

                            end

                        elseif length(tets) == 1

                        else
                            error("length(tets) problem!")
                        end



                        v_new 

                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        push!(fns[m], BEAST.Shape(r, i, v))
                    end
                end
            end

        end

    end

    ntrace(X, ogeo, fns)
#end

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




E[1]
a = faces(E[1])

on_target(a[1])

##


function ntrace(X::Space, γ)

    # on_target = overlap_gpredicate(γ)
    # ad = assemblydata(X)
    x = refspace(X)
    E, ad, P = assemblydata(X)
    igeo = geometry(X)
    @assert dimension(γ) == dimension(igeo)-1
    # Γ = geo
    # Dγ = dimension(γ)
    # Σ = skeleton(Γ,Dγ)

    ogeo = boundary(igeo)
    on_target = overlap_gpredicate(γ)
    ogeo = submesh(ogeo) do m,f
        ch = chart(m,f)
        on_target(ch)
    end

    # D = copy(transpose(connectivity(ogeo, igeo, abs)))
    D = connectivity(igeo, ogeo, abs)
    rows, vals = rowvals(D), nonzeros(D)

    T = scalartype(X)
    S = Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

    for (p,el) in enumerate(E)

        for (q,fc) in enumerate(faces(el))
            on_target(fc) || continue

            # print(Q)
            # @assert norm(Q,Inf) != 0
            
            r = 0
            for k in nzrange(D,P[p])
                vals[k] == q && (r = rows[k]; break)
            end
            @assert r != 0
            
            fc1 = chart(ogeo, r)
            Q = ntrace(x, el, q, fc1)

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[p,j]
                        # j == q && println("bingo",j,q)
                        v = a*Q[i,j]
                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        push!(fns[m], Shape(r, i, v))
                    end
                end
            end

        end

    end

    ntrace(X, ogeo, fns)
end





##






























mutable struct testop
    mf::Vector{Float64}
    v::Float64
end

a = Vector([0.1, 0.9, 1.3])
v_ini = 0.0
op1 = testop(a,v_ini)

function change!(op)
    #new_op = testop(op.mf,50.8)
    #op = new_op
    op.v = 50.8
    return nothing
end

change!(op1)
op1







# Quadstrat
qs3D = BEAST.SingleNumQStrat(6)
qs4D = BEAST.DoubleNumWiltonSauterQStrat(5,5,5,5,6,6,6,6) #BEAST.DoubleNumWiltonSauterQStrat(2,3,2,3,4,4,4,4)
qs5D6D = BEAST.SauterSchwab3DQStrat(5,5,6,6,6,6)

BEAST.defaultquadstrat(op::BEAST.LocalOperator, tfs, bfs) = qs3D
BEAST.defaultquadstrat(op::BEAST.Helmholtz3DOp, tfs, bfs) = qs4D
BEAST.defaultquadstrat(op::BEAST.VIEOperator, tfs, bfs) = qs5D6D