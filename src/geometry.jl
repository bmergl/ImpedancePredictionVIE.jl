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

    p = point(0,0,0)

    @show kappa
    @show kappa0
    @show epsilon
    @show epsilon0
    @show omega
    
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
        tau0 = kappa0 + im*ω*epsilon0 # this ist not ε0 !
        T = ComplexF64

    # elseif problemtype == :dielectic
    #     error("NOT READY")
    #     epsilon === nothing && error("Add epsilon=FUNCTION IN X")
    #     epsilon0 === nothing && error("Add epsilon0=CONSTANT SCALAR")
    #     ((size(epsilon0, 1) != 1) || (size(epsilon0, 2) != 1)) && error("epsilon0 has to be a scalar value!")
    #     !(size(epsilon(p), 1) == size(epsilon(p), 2)) && error("Format error")
        
    #     s = size(epsilon(p), 1)
    #     !((s==1)||(s==3)) && error("epsilon must return skalar 1×1 or tensor 3×3")

    #     s==1 && (I = Float64(1.0)) 
    #     s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

    #     tau = x -> epsilon(x) # kläre jω ?


    # elseif problemtype == :general # epsilon0 must not be ε0! it should be 0.0 (recommended)

    #     error("NOT READY")

    #     kappa === nothing && error("")
    #     epsilon === nothing && error("")

    #     ((epsilon0 === nothing) && (kappa0 === nothing)) && error("Define either kappa0 or epsilon0")
    #     epsilon0 === nothing && (epsilon0 = 0.0)
    #     kappa0 === nothing && (kappa0 = 0.0)

    #     omega === nothing && error("Add omega=VALUE to gen_tau_chi in the general case!")

    #     sk = size(kappa(p), 1)
    #     se = size(epsilon(p), 1)

    #     sk >= se && (s = sk)
    #     se >= sk && (s = se)

    #     !((s==1)||(s==3)) && error("epsilon and kappa must return skalar 1×1 or tensor 3×3")

    #     s==1 && (I = Float64(1.0))
    #     s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

    #     sk > se && (tau = x -> kappa(x) + im*I*omega*epsilon(x))
    #     se > sk && (tau = x -> I*kappa(x) + im*omega*epsilon(x))
    #     se == sk && (tau = x -> kappa(x) + im*omega*epsilon(x))

    #     tau0 = kappa0 + im*ω*epsilon0
        
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



# function gen_tau_chi_const(; kappa = nothing, kappa0 = nothing, epsilon = nothing, epsilon0 = nothing, omega = nothing)

#     @show kappa
#     @show kappa0
#     @show epsilon
#     @show epsilon0
#     @show omega
    
#     if kappa !== nothing && kappa0 !== nothing && epsilon === nothing && epsilon0 === nothing && omega === nothing
#         # :current
#         tau = kappa
#         tau0 = kappa0
    
#     elseif kappa === nothing && kappa0 === nothing && epsilon !== nothing && epsilon0 !== nothing && omega === nothing
#         # :dielectic
#         tau = epsilon
#         tau0 = epsilon0 # this ist not ε0 !

#     elseif kappa !== nothing && kappa0 !== nothing && epsilon !== nothing && epsilon0 !== nothing && omega !== nothing
#         # :general
#         tau = kappa + im*omega*epsilon
#         tau0 = kappa0 + im*ω*epsilon0 # this ist not ε0 !
        
#     else
#         error("Specify problem! current, dielectic or general")
#     end

#     chi = (tau/tau0 - 1.0)*1/tau
#     invtau = 1/tau

#     return tau, invtau, tau0, chi

# end


#########################################################



function gen_tau_invtau(; problemtype = :current, omega = nothing, kappa = nothing, epsilon = nothing, arbitrary_meshpoint = SVector(0.0,0.0,0.0) )
    @warn "kappa(x), epsilon(x) must refer to the mesh, x must be a 3 dimensional vector!"

    p = arbitrary_meshpoint
    
    if problemtype == :current
        kappa === nothing && error("Add kappa=FUNCTION IN X")
        !(size(kappa(p), 1) == size(kappa(p), 2)) && error("Format error")

        s = size(kappa(p), 1)
        !((s==1)||(s==3)) && error("kappa must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0)) 
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))  #DAS I KANN ALLES WEG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        tau = x -> kappa(x)

    elseif problemtype == :dielectic
        error("NOT READY")
        epsilon === nothing && error("Add epsilon=FUNCTION IN X")
        epsilon0 === nothing && error("Add epsilon0=CONSTANT SCALAR")
        ((size(epsilon0, 1) != 1) || (size(epsilon0, 2) != 1)) && error("epsilon0 has to be a scalar value!")
        !(size(epsilon(p), 1) == size(epsilon(p), 2)) && error("Format error")
        
        s = size(epsilon(p), 1)
        !((s==1)||(s==3)) && error("epsilon must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0)) 
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

        tau = x -> epsilon(x) # kläre jω ?


    elseif problemtype == :general # epsilon0 must not be ε0! it should be 0.0 (recommended)

        error("NOT READY")

        kappa === nothing && error("")
        epsilon === nothing && error("")

        ((epsilon0 === nothing) && (kappa0 === nothing)) && error("Define either kappa0 or epsilon0")
        epsilon0 === nothing && (epsilon0 = 0.0)
        kappa0 === nothing && (kappa0 = 0.0)

        omega === nothing && error("Add omega=VALUE to gen_tau_chi in the general case!")

        sk = size(kappa(p), 1)
        se = size(epsilon(p), 1)

        sk >= se && (s = sk)
        se >= sk && (s = se)

        !((s==1)||(s==3)) && error("epsilon and kappa must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0))
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

        sk > se && (tau = x -> kappa(x) + im*I*omega*epsilon(x))
        se > sk && (tau = x -> I*kappa(x) + im*omega*epsilon(x))
        se == sk && (tau = x -> kappa(x) + im*omega*epsilon(x))

        tau0 = kappa0 + im*ω*epsilon0
        
    else
        error("Specify problem! problemtype=:current or :dielectic or :general")
    end

    #chi0 = 1/tau0

    # function chi(x)
    #     #taux = tau(x)
    #     #taurx = taux/tau0
    #     return x -> (kappa(x) + (-1)*I) * 1/kappa(x)
    # end

    # KATHASTROPHE!!!!!!!!!!!!!!!!!!!! der doppelte funktionsaufruf ......
    # function chi(x)
    #     taux = tau(x)
    #     taurx = taux/tau0
    #     return (taurx + (-1)*I) * inv(taux) 
    # end

    # function chi_chi0I(x)   # (χ - χ0*I)
    #     chix = chi(x)
    #     return chix + (-1)*chi0*I
    # end

    # function init_chi(tau0, tau, I)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
    #     function chi(x)
    #         taux = tau(x)
    #         taurx = taux/tau0
    #         return (taurx + (-1.0)*I) * inv(taux)
    #     end
    #     return chi
    # end

    # chi = init_chi(tau0, tau, I)


    # function init_chi_chi0I(chi0, chi, I)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
    #     function chi_chi0I(x)
    #         chix = chi(x)
    #         return chix + (-1.0)*chi0*I         #(χ - χ0*I)
    #     end
    #     return chi_chi0I
    # end
    
    # chi_chi0I = init_chi_chi0I(chi0, chi, I)


    function init_invtau(tau)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
        function invtau(x)
            taux = tau(x)
            return inv(taux)
        end
        return invtau
    end
    invtau = init_invtau(tau)



    
    #chi = x -> kappa(x) + (-1)*chi0 # + (-1)*I) #* 1/kappa(x)
    #chi_chi0I = x -> (kappa(x) + (-1)*I) * 1/kappa(x) + (-1)*chi0*I # 

    #return tau, tau0, chi, chi0, chi_chi0I
    return tau, invtau

end


# Tensor ... aber muss eh überarbeitet werden....
function gen_tau_chi_old(; problemtype = :current, omega = nothing, kappa = nothing, epsilon = nothing, kappa0 = nothing, epsilon0 = nothing, arbitrary_meshpoint = SVector(0.0,0.0,0.0) )
    @warn "kappa(x), epsilon(x) must refer to the mesh, x must be a 3 dimensional vector!"

    p = arbitrary_meshpoint
    
    if problemtype == :current
        kappa === nothing && error("Add kappa=FUNCTION IN X")
        kappa0 === nothing && error("Add kappa0=CONSTANT SCALAR")
        ((size(kappa0, 1) != 1) || (size(kappa0, 2) != 1)) && error("kappa0 has to be a scalar value!")
        !(size(kappa(p), 1) == size(kappa(p), 2)) && error("Format error")

        s = size(kappa(p), 1)
        !((s==1)||(s==3)) && error("kappa must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0)) 
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))  #DAS I KANN ALLES WEG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        tau = x -> kappa(x)
        tau0 = kappa0   #????

        size(kappa(p),1) == 3 && error("NOT READY...ntrace of Tensor*Vector problem...")

    elseif problemtype == :dielectic
        error("NOT READY")
        epsilon === nothing && error("Add epsilon=FUNCTION IN X")
        epsilon0 === nothing && error("Add epsilon0=CONSTANT SCALAR")
        ((size(epsilon0, 1) != 1) || (size(epsilon0, 2) != 1)) && error("epsilon0 has to be a scalar value!")
        !(size(epsilon(p), 1) == size(epsilon(p), 2)) && error("Format error")
        
        s = size(epsilon(p), 1)
        !((s==1)||(s==3)) && error("epsilon must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0)) 
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

        tau = x -> epsilon(x) # kläre jω ?


    elseif problemtype == :general # epsilon0 must not be ε0! it should be 0.0 (recommended)

        error("NOT READY")

        kappa === nothing && error("")
        epsilon === nothing && error("")

        ((epsilon0 === nothing) && (kappa0 === nothing)) && error("Define either kappa0 or epsilon0")
        epsilon0 === nothing && (epsilon0 = 0.0)
        kappa0 === nothing && (kappa0 = 0.0)

        omega === nothing && error("Add omega=VALUE to gen_tau_chi in the general case!")

        sk = size(kappa(p), 1)
        se = size(epsilon(p), 1)

        sk >= se && (s = sk)
        se >= sk && (s = se)

        !((s==1)||(s==3)) && error("epsilon and kappa must return skalar 1×1 or tensor 3×3")

        s==1 && (I = Float64(1.0))
        s==3 && (I = SMatrix{3,3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))

        sk > se && (tau = x -> kappa(x) + im*I*omega*epsilon(x))
        se > sk && (tau = x -> I*kappa(x) + im*omega*epsilon(x))
        se == sk && (tau = x -> kappa(x) + im*omega*epsilon(x))

        tau0 = kappa0 + im*ω*epsilon0
        
    else
        error("Specify problem! problemtype=:current or :dielectic or :general")
    end

    function init_chi(tau0, tau, I)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
        function chi(x)
            taux = tau(x)
            taurx = taux/tau0
            return (taurx + (-1.0)*I) * inv(taux)
        end
        return chi
    end
    chi = init_chi(tau0, tau, I)

    function init_invtau(tau)     #SO IST ES SCHNELL!!!!!!!!! WIESO?????? FOLGLICH ALLE EINZELN...
        function invtau(x)
            taux = tau(x)
            return inv(taux)
        end
        return invtau
    end
    invtau = init_invtau(tau)


    return tau, invtau, tau0, chi

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


function gen_w_mat(w::BEAST.LagrangeBasis{0,-1}, X::BEAST.NDLCDBasis, cell2mat_inv_τ::Vector)

    # erstelle w_mat aus w mittels tri->tet, cell2mat_inv_τ:
    D = connectivity(w.geo, X.geo)
    @assert sum(D) == length(w.fns) == length(w.geo.faces)
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

    D = connectivity(igeo, ogeo, abs) # nzrange(D,tetnummer) liefert die 4 bzw. seltener 3 Dreiecke
    rows, vals = rowvals(D), nonzeros(D)
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
        
        # Wähle einen Tetraeder aus von dem aus die innere ntrace betrachtet wird -  Achtung n̂ ! 
        
        i_el = tets[1] # der zweite wird gar nicht betrachtet für den Fall dass es existiert
        el = E[i_el] 
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
        Q = ntrace(x, el, q, fc1)

        for i in 1:size(Q,1)
            for j in 1:size(Q,2)
                for (m,a) in ad[i_el,j] # das ist die assemblydata der Basis X

                    v = a*Q[i,j]

                    n̂_n = fc1.normals[1]

                    # Betrachte die 1-2 angrenzenden Tetraeder von fc1
                    if length(tets) == 2
                        i_tet1 = tets[1]
                        i_tet2 = tets[2]
                        tet1 = E[i_tet1]
                        tet2 = E[i_tet2]
                        center1 = cartesian(CompScienceMeshes.center(tet1))
                        center2 = cartesian(CompScienceMeshes.center(tet2))
                        v_f1 = center1 - fc1_center
                        v_f2 = center2 - fc1_center

                        if dot(v_f1, n̂_n) > 0.0
                            tet_minus = tet1
                            tet_plus = tet2
                            i_tet_minus = i_tet1 
                            i_tet_plus = i_tet2
                            @assert dot(v_f2, n̂_n) < 0.0
                        elseif dot(v_f1, n̂_n) < 0.0
                            tet_minus = tet2
                            tet_plus = tet1
                            i_tet_minus = i_tet2
                            i_tet_plus = i_tet1
                            @assert dot(v_f2, n̂_n) > 0.0
                        end
                        χ_minus = cell2mat_χ[i_tet_minus]
                        χ_plus = cell2mat_χ[i_tet_plus]
                        δχ = χ_minus - χ_plus

                    elseif length(tets) == 1
                        i_tet1 = tets[1]
                        tet1 = E[i_tet1]
                        center1 = cartesian(CompScienceMeshes.center(tet1))
                        v_f1 = center1 - fc1_center
                        
                        if dot(v_f1, n̂_n) > 0.0
                            tet_minus = tet1
                            tet_plus = nothing
                            i_tet_minus = i_tet1 
                            i_tet_plus = nothing
                            χ_minus = cell2mat_χ[i_tet_minus]
                            χ_plus = 0.0
                        elseif dot(v_f1, n̂_n) < 0.0
                            tet_minus = nothing
                            tet_plus = tet1
                            i_tet_minus = nothing
                            i_tet_plus = i_tet1
                            χ_minus = 0.0
                            χ_plus = cell2mat_χ[i_tet_plus]
                        end
                        δχ = χ_minus - χ_plus

                    else
                        error("length(tets) problem!")
                    end
                    #@show δχ
                    v_new = v * δχ
                    isapprox(v,0,atol=sqrt(eps(T))) && continue # WICHTIG sonst unnötig viele Nullen aber ohne δχ, sonst für kleine ... ggf Probleme
                    push!(fns[m], BEAST.Shape(s, i, v_new))
                end
            end
        end

    end

    return BEAST.ntrace(X, ogeo, fns)
end





