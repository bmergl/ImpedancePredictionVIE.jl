# function simpletest()
#     @show "Hallo"
#     return nothing
# end


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

    
    if kappa !== nothing && kappa0 !== nothing && epsilon === nothing && epsilon0 === nothing && omega === nothing
        # :current
        tau = x -> kappa(x)
        tau0 = kappa0
    
    else if kappa === nothing && kappa0 === nothing && epsilon !== nothing && epsilon0 !== nothing && omega === nothing
        # :dielectic
        tau = x -> epsilon(x)
        tau0 = epsilon0 # this ist not ε0 !

    else if kappa !== nothing && kappa0 !== nothing && epsilon !== nothing && epsilon0 !== nothing && omega !== nothing
        # :general
        tau = x -> kappa(x) + im*omega*epsilon(x)
        tau0 = kappa0 + im*ω*epsilon0 # this ist not ε0 !

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


    return tau, invtau, tau0, chi

end




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