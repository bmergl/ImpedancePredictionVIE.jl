using SparseArrays



function geo2mesh(geopath::String, meshpath::String, meshparam::Float64; body::IP.geometric_body)


    h = meshparam
    if typeof(body) == IP.cuboid 
        L_x = body.L_x
        L_y = body.L_y
        L_z = body.L_z
        run(`gmsh $geopath -3 -clmax $h -format msh2 -setnumber lx $L_x -setnumber ly $L_y -setnumber lz $L_z -o $meshpath`)
    else
        #general formulation, use values from gmsh skript
        run(`gmsh $geopath -3 -clmax $h -format msh2 -setnumber var 55 -o $meshpath`) 
    end
    
    Ω = CompScienceMeshes.read_gmsh3d_mesh("$meshpath", physical="BodyVolume")

    Γ_c_t = CompScienceMeshes.read_gmsh_mesh("$meshpath", physical= "TopElectrode")
    Γ_c_b = CompScienceMeshes.read_gmsh_mesh("$meshpath", physical= "BottomElectrode")
    Γ_nc = CompScienceMeshes.read_gmsh_mesh("$meshpath", physical= "InsulatingContact")

    #Γ = weld(Γ_c_t, Γ_c_b, Γ_nc) !!!! hat die vertices nicht zusammengeführt!!!

    faces = deepcopy(Γ_c_t.faces) #weld in CSM als Alternatives! "verschweißen"
    append!(faces, Γ_c_b.faces)
    append!(faces, Γ_nc.faces)
    Γ = Mesh(Ω.vertices, faces)

    faces2 = deepcopy(Γ_c_t.faces)
    append!(faces2, Γ_c_b.faces)
    Γ_c = Mesh(Ω.vertices, faces2)

    @warn "Assumed that .geo has correct format: BodyVolume, TopE...!"

    return Ω, Γ, Γ_c, Γ_c_t, Γ_c_b, Γ_nc
end


function realvertices(mesh::Mesh) # mesh.vertices includes vertices that are not used but impotant to handle...
    active_nodes = realnodes(mesh)
    return mesh.vertices[active_nodes]
end

function realnodes(mesh::Mesh)
    active_nodes = Int64[]
    for n in skeleton(mesh, 0).faces
        push!(active_nodes, n[1])
    end
    return active_nodes
end

function pointlist2xyzlist(pointlist)

    L = length(pointlist)

    allx = Vector{Float64}(undef, L)
    ally = Vector{Float64}(undef, L)
    allz = Vector{Float64}(undef, L)

    for (i,p) in enumerate(pointlist)
        allx[i] = p[1]
        ally[i] = p[2]
        allz[i] = p[3]
    end

    return allx, ally, allz
end


# function dirichlet_n2f(y_d::BEAST.LagrangeBasis{1,0}, dirichletnodes)           #<--------SPEED? Distanz Octree...
#     # input: glob node nr. 12 => glob fns nr. = n2f[12]   !!! geht nur für knotenbasierte fns...
#     n2f = Int64[]
#     vertices = y_d.geo.vertices[dirichletnodes]
#     for mvert in vertices
#         for (j,Xvert) in enumerate(y_d.pos)
#             if norm(mvert - Xvert) < 1e-12
#                 push!(n2f, j)
#                 break
#             end
#         end
#     end
#     @assert length(dirichletnodes) == length(n2f)
#     return n2f # Vektor
# end

function swgfaces(volmesh::Mesh, ncbndmesh::Mesh; fast = true)           #<--------SPEED? Distanz Octree...
    
    if fast == false

        all_faces = skeleton(volmesh,2).faces
        excluded_faces = ncbndmesh.faces
    
        @assert length(all_faces[1]) == 3
        @assert length(excluded_faces[1]) == 3
        @assert volmesh.vertices == ncbndmesh.vertices
        @assert length(all_faces) > length(excluded_faces)

        faces = Vector{SVector{3, Int64}}()
        for f1 in all_faces
            addf1 = true 

            for f2 in excluded_faces
                ap = collect(BEAST.permutations(f2)) # circshift is not enough!!!

                if (f1 == ap[1])||(f1 == ap[2])||(f1 == ap[3])||(f1 == ap[4])||(f1 == ap[5])||(f1 == ap[6])
                    addf1 = false
                    break
                end
            end

            addf1 && push!(faces, f1)
        end


        @assert length(all_faces) == length(excluded_faces) + length(faces)
        return faces
    end


    all_faces_mesh = skeleton(volmesh,2)
    l = length(all_faces_mesh.faces)
    #all_charts = Vector{CompScienceMeshes.Simplex{3, 2, 1, 3, Float64}}(undef, l)
    all_charts_centers = Vector{SVector{3, Float64}}(undef, l)
    for i in 1:l
        chart = CompScienceMeshes.chart(all_faces_mesh, i)      #ÜBER assemblydata!!!
        #all_charts[i] = chart

        center = CompScienceMeshes.center(chart)
        all_charts_centers[i] = cartesian(center)
    end

    l_ = length(ncbndmesh.faces)
    excluded_charts = Vector{CompScienceMeshes.Simplex{3, 2, 1, 3, Float64}}(undef, l_)
    #excluded_charts_centers = Vector{SVector{3, Float64}}(undef, l_)
    for i in 1:l_
        chart = CompScienceMeshes.chart(ncbndmesh, i)   #ÜBER assemblydata!!!
        excluded_charts[i] = chart 

        #center = CompScienceMeshes.center(chart)
        #excluded_charts_centers[i] = cartesian(center)
    end

    swg_faces = Vector{SVector{3, Int64}}()
    excluded_charts_tree = BEAST.octree(excluded_charts)

    @assert length(all_charts_centers) >  length(excluded_charts)
    @assert length(all_charts_centers) ==  length(all_faces_mesh.faces)

    for (j,pos) in enumerate(all_charts_centers)

        i = CompScienceMeshes.findchart(excluded_charts, excluded_charts_tree, pos)

        if i === nothing    # nicht gefunden => face hinzufügen
            push!(swg_faces, all_faces_mesh.faces[j])
        else
            # gefunden! => face liegt auf Γ_nc und muss ausgeschlossen werden!
        end
    end

    @assert length(all_charts_centers) == length(excluded_charts) + length(swg_faces)
    
    return swg_faces
end













function gen_tau_chi(; kappa = nothing, kappa0 = nothing, epsilon = nothing, epsilon0 = nothing, omega = nothing)
    @warn "kappa(x), epsilon(x) must refer to the mesh, x must be a 3 dimensional vector!"

    
    if kappa !== nothing && kappa0 !== nothing && epsilon === nothing && epsilon0 === nothing && omega === nothing
        # :current
        tau = x -> kappa(x)
        tau0 = kappa0
        T = Float64
    
    elseif kappa === nothing && kappa0 === nothing && epsilon !== nothing && epsilon0 !== nothing && omega === nothing
        # :dielectic
        tau = x -> epsilon(x)
        tau0 = epsilon0 # this ist not ε0 !
        T = Float64

    elseif kappa !== nothing && kappa0 !== nothing && epsilon !== nothing && epsilon0 !== nothing && omega !== nothing
        # :general
        tau = x -> kappa(x) + im*omega*epsilon(x)
        tau0 = kappa0 + im*omega*epsilon0 # this ist not ε0 ! sollte omega da stehen???
        T = ComplexF64
        
    else
        error("Specify problem! current, dielectic or general")
    end

    function init_chi(tau0, tau)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
        function chi(x)
            taux = tau(x)
            taurx = taux/tau0
            return (taurx - 1.0) * inv(taux)
        end
        return chi
    end
    chi = init_chi(tau0, tau)

    function init_invtau(tau)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
        function invtau(x)
            taux = tau(x)
            return 1/taux
        end
        return invtau
    end
    invtau = init_invtau(tau)


    return tau, invtau, tau0, chi, T

end




#chart(mesh::Mesh, cell) ist simplex, davon kann man center machen
# s = simplex(volmesh.vertices[f1])
# cent = CompScienceMeshes.center(s)
# addf1 && @show cartesian(cent)    

# for f1 in excluded_faces
#     err = true
#     for f2 in all_faces
#         ap = collect(BEAST.permutations(f2))
#         @assert length(ap) == 6
#         if (f1 == ap[1]) || (f1 == ap[2]) ||(f1 == ap[3]) ||(f1 == ap[4]) ||(f1 == ap[5]) ||(f1 == ap[6])
#             err = false
#         end
#     end
#     err && @show f1
#     err && error("Diese excluded face existiert in all_faces gar nicht")
# end








######## PWC MATERIAL ################################

function gen_cell2mat(τ, inv_τ, τ0, χ, T, X::BEAST.NDLCDBasis)

    elementsX, adX, cellsX = assemblydata(X)
    cell2mat_inv_τ = Vector{T}(undef,length(elementsX))
    cell2mat_χ  = Vector{T}(undef,length(elementsX)) # Typenunterscheidung???

    for (n,el) in enumerate(elementsX)
        center = cartesian(CompScienceMeshes.center(el))
        inv_τ_n = inv_τ(center)
        χ_n = χ(center)
        cell2mat_inv_τ[n] = inv_τ_n 
        cell2mat_χ[n] = χ_n   
    end

    return cell2mat_inv_τ, cell2mat_χ
end


function gen_X_mat(X::BEAST.NDLCDBasis, cell2mat::Vector) # cell2mat_χ or cell2mat_inv_τ

    T = typeof(cell2mat[1])
    newfns = Vector{Vector{BEAST.Shape{T}}}() 
    for (i,shs) in enumerate(X.fns)
        newshs = Vector{BEAST.Shape{T}}()
        for (j,sh) in enumerate(shs)
            cellid = sh.cellid
            refid = sh.refid
            coeff = sh.coeff

            mat_of_cell = cell2mat[cellid]
            new_coeff = coeff*mat_of_cell # can be ComplexF64
            push!(newshs, BEAST.Shape(cellid, refid, new_coeff))
        end
        push!(newfns, newshs)
    end
    X_mat = BEAST.NDLCDBasis(X.geo, newfns, X.pos)

    return X_mat
end


function gen_ntrcX_mat(ntrcX::BEAST.LagrangeBasis{0, -1}, cell2mat_inv_τ::Vector)

    T = typeof(cell2mat_inv_τ[1])
    newfns = Vector{Vector{BEAST.Shape{T}}}()
    for (i,shs) in enumerate(ntrcX.fns)
        newshs = Vector{BEAST.Shape{T}}()
        shs == [] && (push!(newfns, newshs); continue)

        for (j,sh) in enumerate(shs)
            cellid = sh.cellid
            refid = sh.refid
            coeff = sh.coeff

            mat_of_cell = cell2mat_inv_τ[cellid]
            new_coeff = coeff*mat_of_cell # can be ComplexF64
            push!(newshs, BEAST.Shape(cellid, refid, new_coeff))
        end
        push!(newfns, newshs)
    end

    ntrcX_mat = BEAST.LagrangeBasis{0,-1,1}(ntrcX.geo, newfns, deepcopy(ntrcX.pos))

    return ntrcX_mat
    #ntrace(X::NDLCDBasis, geo, fns) = LagrangeBasis{0,-1,1}(geo, fns, deepcopy(X.pos))  
end


function gen_w_mat(w::BEAST.LagrangeBasis{0,-1}, X::BEAST.NDLCDBasis, cell2mat_inv_τ::Vector)

    # erstelle w_mat aus w mittels tri->tet, cell2mat_inv_τ:
    D = connectivity(w.geo, X.geo)
    rows, vals = SparseArrays.rowvals(D), SparseArrays.nonzeros(D)

    face2cell = rows # Bsp rows[i]=j vals[i]=a =>D[j,i]=a ... Nur hier so einfach!!!

    T = typeof(cell2mat_inv_τ[1])
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
    
    return w_mat
end


function inner_mat_ntrace(X::BEAST.NDLCDBasis, swg_faces_mesh::Mesh, cell2mat_χ::Vector)

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

    # D = connectivity(igeo, ogeo, abs) # nzrange(D,tetnummer) liefert die 4 bzw. seltener 3 Dreiecke
    # rows, vals = rowvals(D), nonzeros(D)
    Dt = connectivity(ogeo, igeo, abs) # nzrange(Dt,trinummer) liefert die 2 bzw. seltener 1 Tetraeder
    rowst, valst = rowvals(Dt), nonzeros(Dt)


    T = typeof(cell2mat_χ[1])
    S = BEAST.Shape{T}
    fns = [Vector{S}() for i in 1:numfunctions(X)]

    for s in 1:numfunctions(X)
        fc1 = chart(ogeo, s)

        # 1-2 angrenzende Tetraeder bestimmen
        tets = Vector{Int64}()
        for k in nzrange(Dt,s) # s ist der index von fc1
            push!(tets,rowst[k])
        end

        if length(tets) == 2

            i_tet1 = tets[1]
            i_tet2 = tets[2]
            tet1 = E[i_tet1]
            tet2 = E[i_tet2]
            center1 = cartesian(CompScienceMeshes.center(tet1))
            center2 = cartesian(CompScienceMeshes.center(tet2))
            fc1_center = cartesian(CompScienceMeshes.center(fc1))
            v_f1 = center1 - fc1_center
            v_f2 = center2 - fc1_center
            n̂_n = fc1.normals[1]
            if dot(v_f1, n̂_n) > 0.0
                i_tet_minus = i_tet1 # normale zeigt in tet1 => tet1=tet_minus
                i_tet_plus = i_tet2
                @assert dot(v_f2, n̂_n) < 0.0
            elseif dot(v_f1, n̂_n) < 0.0
                i_tet_minus = i_tet2
                i_tet_plus = i_tet1
                @assert dot(v_f2, n̂_n) > 0.0
            end
            χ_minus = cell2mat_χ[i_tet_minus]
            χ_plus = cell2mat_χ[i_tet_plus]
            δχ = χ_minus - χ_plus

            el = E[i_tet_plus] 
            fc1_center = cartesian(CompScienceMeshes.center(fc1))
            q = nothing
            fc = nothing
            for (k,fc_) in enumerate(faces(el))
                fc_center = cartesian(CompScienceMeshes.center(fc_))
                if isapprox(norm(fc_center-fc1_center),0,atol=sqrt(eps(real(T))))
                    q = k
                    fc = fc_
                    break
                end
            end
            Q = ntrace(x, el, q, fc1) # VORSICHT, wenn man dieses ntrace verwendet MUSS man mit dem tet arbeiten wo n̂_fc1 nach AUßEN zeigt!!!!

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[i_tet_plus,j] # das ist die assemblydata der Basis X
                        v = a*Q[i,j] # hier kommt das VZ mit rein vmtl. hätte man also auch auf dem anderen tet arbeiten können... aber nicht auf beiden!
                        v_new = v * δχ

                        isapprox(v,0,atol=sqrt(eps(real(T)))) && continue
                        isapprox(v_new,0,atol=sqrt(eps(real(T)))/100000) && continue
                        push!(fns[m], BEAST.Shape(s, i, v_new))
                    end
                end
            end

        elseif length(tets) == 1
            continue # !!! Muss vmlt verschwinden constantmaterial(2.0) => u_Jn Fehler viel niedriger

            i_tet = tets[1]
            tet = E[i_tet]

            center = cartesian(CompScienceMeshes.center(tet))
            fc1_center = cartesian(CompScienceMeshes.center(fc1))
            v_f = center - fc1_center
            n̂_n = fc1.normals[1]
            if dot(v_f, n̂_n) > 0.0
                i_tet_minus = i_tet # normale zeigt in tet1 => tet1=tet_minus
                χ_minus = cell2mat_χ[i_tet_minus]
                χ_plus = 0.0
            elseif dot(v_f, n̂_n) < 0.0
                i_tet_plus = i_tet
                χ_minus = 0.0
                χ_plus = cell2mat_χ[i_tet_plus]
            end
            δχ = χ_minus - χ_plus

            el = E[i_tet] 
            fc1_center = cartesian(CompScienceMeshes.center(fc1))
            q = nothing
            fc = nothing
            for (k,fc_) in enumerate(faces(el))
                fc_center = cartesian(CompScienceMeshes.center(fc_))
                if isapprox(norm(fc_center-fc1_center),0,atol=sqrt(eps(T)))
                    q = k
                    fc = fc_
                    break
                end
            end

            @show dot(fc.normals[1], n̂_n) 
            Q = ntrace(x, el, q, fc1)

            for i in 1:size(Q,1)
                for j in 1:size(Q,2)
                    for (m,a) in ad[i_tet,j] # das ist die assemblydata der Basis X
                        v = a*Q[i,j] # hier kommt das VZ mit rein vmtl. hätte man also auch auf dem anderen tet arbeiten können... aber nicht auf beiden!
                        v_new = v * δχ

                        isapprox(v,0,atol=sqrt(eps(T))) && continue
                        isapprox(v_new,0,atol=sqrt(eps(T))/100000) && continue
                        push!(fns[m], BEAST.Shape(s, i, v_new))
                    end
                end
            end



        else
            error()
        end

    end


    return BEAST.ntrace(X, ogeo, fns)
end





