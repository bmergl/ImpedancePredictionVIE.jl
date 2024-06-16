
####################################################################################
# SauterSchwab3D and the corresponding BEAST.momintegrals! work for ∫∫_∂Ω ∫∫∫_Ω and ∫∫∫_Ω ∫∫∫_Ω but NOT for  ∫∫∫_Ω ∫∫_∂Ω
# This can be fixed by split BEAST.momintegrals! into three parts, one for every new operator type.  

# strat.sing.T für Testkörper, strat.sing.S für Basiskörper,
# reorder(sing::Singularity5D...) hat IMMER return TET, TRI => Bei sing. kann es vertauscht sein!
####################################################################################



#kernelvals(igd.op,tgeo,bgeo) hier kann man für einen spezifischen operator
# ... die Dyade einbauen um die dann in integrand() zu benutzen!


abstract type BoundaryOperatorΩΓ <: BEAST.BoundaryOperator end
abstract type BoundaryOperatorΓΩ <: BEAST.BoundaryOperator end
abstract type VolumeOperatorΩΩ <: BEAST.VolumeOperator end





#####  kernelvals für Dyadische VIE Operatoren  ######################################################################
struct KernelValsVIEdyad{T,U,P,Q,K} # ÄNDERN!!!!!!
    gamma::U
    vect::P
    dist::T
    dyadgreen::Q
    tau::K
end

function kernelvalsdyad(viop, p ,q) # p=r_vec, q=r'_vec
    # Achtung! Speziell auf gamma=0 zugeschnitten, gamma für typ wichtig
    Y = viop.gamma  #ComplexF64/Float64 unterscheidung läuft normalerweise über Gamma...
    r = cartesian(p)-cartesian(q)
    R = norm(r)
    Rsq = R^2

    p_ = cartesian(p)
    q_ = cartesian(q)
    xd = p_[1]-q_[1]
    yd = p_[2]-q_[2]
    zd = p_[3]-q_[3]

    dyadgreen =  @SMatrix [3*xd^2-Rsq xd*yd xd*zd;  
    yd*xd 3*yd^2-Rsq yd*zd;
    zd*xd zd*yd 3*zd^2-Rsq]
    dyadgreen = dyadgreen * (1/(4*pi*R^5)) # = ∇'∇'G


    tau = viop.tau(cartesian(q))

    KernelValsVIEdyad(Y,r,R, dyadgreen, tau)
end






##### momintegrals!  ######################################################################

###### Sonderintegranden für den CommonFace6D Fall
struct VIEIntegrandΩΩ_cf6d{S,T,M,O,K,L}
    test_tetrahedron_element::S
    trial_tetrahedron_element::T
    trial_tetrahedron_element_inv::M
    op::O
    test_local_space::K
    trial_local_space::L
end
function (igd::VIEIntegrandΩΩ_cf6d)(u,v) # spezielle reoder_dof????

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,u)
    bgeo = neighborhood(igd.trial_tetrahedron_element,v)
    bgeo_refspace = neighborhood(igd.trial_tetrahedron_element_inv,v)

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo)
    #bval = igd.trial_local_space(bgeo) <---- geht nicht, denn bgeo ist MP mit non-CSM tet
    bval = igd.trial_local_space(bgeo_refspace)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,bgeo) * j
end



#######

# 5D: ∫∫∫_Ω ∫∫_Γ      -> Achtung, hier haben wir keinen Unit test!!!
struct VIEIntegrandΩΓ{S,T,O,K,L}
    test_tetrahedron_element::S
    trial_tetrahedron_element::T
    op::O
    test_local_space::K
    trial_local_space::L
end

function (igd::VIEIntegrandΩΓ)(u,v)     # wird nur für den SS3D-Fall verwendet, d.h. Kontakt zw. test/trial simplex

    @assert length(igd.test_tetrahedron_element.vertices) == 4
    @assert length(igd.trial_tetrahedron_element.vertices) == 3
    @assert length(u) == 3
    @assert length(v) == 2

    #@show u, v
    #error("Ja, das ist VIEIntegrandΩΓ")

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,u) #!!!
    bgeo = neighborhood(igd.trial_tetrahedron_element,v) #!!!

    @assert length(tgeo.patch.vertices) == 4
    @assert length(bgeo.patch.vertices) == 3

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo) 
    bval = igd.trial_local_space(bgeo)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,bgeo) * j
end

function BEAST.qr_boundary(op::BoundaryOperatorΩΓ, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
    qs::BEAST.SauterSchwab3DQStrat)
    
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    hits = 0
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    D = dimension(τ)+dimension(σ)
    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    @assert hits <= 3
    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   
    # Ohne die simplex Methode!
    # hits == 3 && return SauterSchwab3D.CommonFace5D(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(SauterSchwab3D._legendre(acc,0.0,1.0)))
    # hits == 2 && return SauterSchwab3D.CommonEdge5D(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(SauterSchwab3D._legendre(acc,0.0,1.0)))
    # hits == 1 && return SauterSchwab3D.CommonVertex5D(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(SauterSchwab3D._legendre(acc,0.0,1.0)))
    #hits == 3 && return SauterSchwab3D.CommonFace5D(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(qd.sing_qp[1]))
    #hits == 2 && return SauterSchwab3D.CommonEdge5D(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(qd.sing_qp[1]))
    #hits == 1 && return SauterSchwab3D.CommonVertex5D(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(qd.sing_qp[1]))


    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_t,idx_s),(qd.sing_qp[3],qd.sing_qp[2]))


    return BEAST.DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end

function BEAST.momintegrals!(op::BoundaryOperatorΩΓ,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_tetrahedron_element, trial_tetrahedron_element, out, strat::SauterSchwab3DStrategy)

    #rand()< 0.05 && @show "momintegrals!(op::BoundaryOperatorΩΓ,...) aufgerufen!!!"

    #Find permutation of vertices to match location of singularity to SauterSchwab
    I, J = SauterSchwab3D.reorder(strat.sing)
      
    #Get permutation and rel. orientatio of DoFs 
    K,O1 = BEAST.reorder_dof(test_local_space, I)
    L,O2 = BEAST.reorder_dof(trial_local_space, J)


    #Apply permuation to elements
    test_tetrahedron_element  = simplex(
            test_tetrahedron_element.vertices[I[1]],
            test_tetrahedron_element.vertices[I[2]],
            test_tetrahedron_element.vertices[I[3]],
            test_tetrahedron_element.vertices[I[4]])
    trial_tetrahedron_element  = simplex(
            trial_tetrahedron_element.vertices[J[1]],
            trial_tetrahedron_element.vertices[J[2]],
            trial_tetrahedron_element.vertices[J[3]])


    #Define integral (returns a function that only needs barycentric coordinates)
    igd = VIEIntegrandΩΓ(test_tetrahedron_element, trial_tetrahedron_element,
    op, test_local_space, trial_local_space)
    
    #Evaluate integral
    Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat) 
    
    # Es wird erwartet f((1-x1,x1-x1*x2,x1*x2*x3),(1-x1*x4,x1*x4*x5)) für die igd function d.h. igd(3D,2D)
    # Zudem nutzt  sauterschwab_parameterized(integrand, method) nur noch qps = method.qps
    # ... also ist die Reihenfolge strat.sing.T bzw. .S egal.

    #Undo permuation on DoFs
    for j in 1 : length(L)
        for i  in 1 : length(K)
            out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
        end
    end
    nothing
end







# 5D: ∫∫_Γ ∫∫∫_Ω 

function BEAST.qr_boundary(op::BoundaryOperatorΓΩ, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
    qs::BEAST.SauterSchwab3DQStrat)  # T <-> S
    #error("Ja das richtige qr boundary")
    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    hits = 0
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    D = dimension(τ)+dimension(σ)
    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            d2 = LinearAlgebra.norm_sqr(t-s)
            dmin2 = min(dmin2, d2)
            if d2 < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
   
    # T <-> S  : 
    hits == 3 && return SauterSchwab3D.CommonFace5D_S(SauterSchwab3D.Singularity5DFace(idx_s,idx_t),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 2 && return SauterSchwab3D.CommonEdge5D_S(SauterSchwab3D.Singularity5DEdge(idx_s,idx_t),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    hits == 1 && return SauterSchwab3D.CommonVertex5D_S(SauterSchwab3D.Singularity5DPoint(idx_s,idx_t),(qd.sing_qp[3],qd.sing_qp[2]))


    return BEAST.DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end

function BEAST.momintegrals!(op::BoundaryOperatorΓΩ,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_tetrahedron_element, trial_tetrahedron_element, out, strat::SauterSchwab3DStrategy)
    
    #error("Ja, das ist momintegrals für BoundaryOperatorΓΩ")

    #Find permutation of vertices to match location of singularity to SauterSchwab
    J, I = SauterSchwab3D.reorder(strat.sing) #Singularity...(idx_s,idx_t)
      
    #Get permutation and rel. orientatio of DoFs 
    K,O1 = BEAST.reorder_dof(test_local_space, I)
    L,O2 = BEAST.reorder_dof(trial_local_space, J)

    #Apply permuation to elements
    test_tetrahedron_element  = simplex(
            test_tetrahedron_element.vertices[I[1]],
            test_tetrahedron_element.vertices[I[2]],
            test_tetrahedron_element.vertices[I[3]])
    trial_tetrahedron_element  = simplex(
            trial_tetrahedron_element.vertices[J[1]],
            trial_tetrahedron_element.vertices[J[2]],
            trial_tetrahedron_element.vertices[J[3]],
            trial_tetrahedron_element.vertices[J[4]])

    #Define integral (returns a function that only needs barycentric coordinates)
    igd = BEAST.VIEIntegrand(test_tetrahedron_element, trial_tetrahedron_element,
        op, test_local_space, trial_local_space) # standard
    
    #Evaluate integral
    Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat) 
  
    #Undo permuation on DoFs
    for j in 1 : length(L)
        for i  in 1 : length(K)
            out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
        end
    end
    nothing
end


# 6D: ∫∫∫_Ω ∫∫∫_Ω
struct VIEIntegrandΩΩ{S,T,O,K,L}
    test_tetrahedron_element::S
    trial_tetrahedron_element::T
    op::O
    test_local_space::K
    trial_local_space::L
end
function (igd::VIEIntegrandΩΩ)(u,v)     # wird nur für den SS3D-Fall verwendet, d.h. Kontakt zw. test/trial simplex

    @assert length(igd.test_tetrahedron_element.vertices) == 4
    @assert length(igd.trial_tetrahedron_element.vertices) == 4
    @assert length(u) == 3
    @assert length(v) == 3

    #@show u, v
    #error("JA, das ist 6D: ∫∫∫_Ω ∫∫∫_Ω")

    #mesh points
    tgeo = neighborhood(igd.test_tetrahedron_element,u) #!!! auch mal für vol teil hier tauschen...
    bgeo = neighborhood(igd.trial_tetrahedron_element,v) #!!!

    @assert length(tgeo.patch.vertices) == 4
    @assert length(bgeo.patch.vertices) == 4

    #kernel values
    kerneldata = kernelvals(igd.op,tgeo,bgeo)

    #values & grad/div/curl of local shape functions
    tval = igd.test_local_space(tgeo) 
    bval = igd.trial_local_space(bgeo)

    #jacobian
    j = jacobian(tgeo) * jacobian(bgeo)
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,bgeo) * j
end
function BEAST.qr_volume(op::VolumeOperatorΩΩ, g::RefSpace, f::RefSpace, i, τ, j, σ, qd,
    qs::BEAST.SauterSchwab3DQStrat)

    dtol = 1.0e3 * eps(eltype(eltype(τ.vertices)))

    hits = 0
    idx_t = Int64[]
    idx_s = Int64[]
    sizehint!(idx_t,4)
    sizehint!(idx_s,4)
    dmin2 = floatmax(eltype(eltype(τ.vertices)))
    D = dimension(τ)+dimension(σ)
    for (i,t) in enumerate(τ.vertices)
        for (j,s) in enumerate(σ.vertices)
            d2 = LinearAlgebra.norm_sqr(t-s)
            d = norm(t-s)
            dmin2 = min(dmin2, d2)
            # if d2 < dtol
            if d < dtol
                push!(idx_t,i)
                push!(idx_s,j)
                hits +=1
                break
            end
        end
    end

    #singData = SauterSchwab3D.Singularity{D,hits}(idx_t, idx_s )
    @assert hits <= 4
    hits == 4 && return SauterSchwab3D.CommonVolume6D_S(SauterSchwab3D.Singularity6DVolume(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[4]))
    hits == 3 && return SauterSchwab3D.CommonFace6D_S(SauterSchwab3D.Singularity6DFace(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3]))
    #hits == 2 && return SauterSchwab3D.CommonEdge6D_S(SauterSchwab3D.Singularity6DEdge(idx_t,idx_s),(qd.sing_qp[1],qd.sing_qp[2],qd.sing_qp[3],qd.sing_qp[4]))
    #hits == 2 && return SauterSchwab3D.CommonEdge6D(SauterSchwab3D.Singularity6DEdge(idx_t,idx_s),(SauterSchwab3D._legendre(???,0.0,1.0)))
    hits == 2 && return SauterSchwab3D.CommonEdge6D(SauterSchwab3D.Singularity6DEdge(idx_t,idx_s),(qd.sing_qp[1])) # KEIN _S !!!!!!
    hits == 1 && return SauterSchwab3D.CommonVertex6D_S(SauterSchwab3D.Singularity6DPoint(idx_t,idx_s),qd.sing_qp[3])


    return BEAST.DoubleQuadRule( # auch bei hits = 2 => wähle andere Quadraturordnung!!! nicht nötig in dem Fall
        qd[1][1,i],
        qd[2][1,j])

end
function BEAST.momintegrals!(op::VolumeOperatorΩΩ,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_tetrahedron_element, trial_tetrahedron_element, out, strat::SauterSchwab3DStrategy)

    #Find permutation of vertices to match location of singularity to SauterSchwab
    I, J = SauterSchwab3D.reorder(strat.sing)
      
    #Get permutation and rel. orientatio of DoFs 
    K,O1 = BEAST.reorder_dof(test_local_space, I)
    L,O2 = BEAST.reorder_dof(trial_local_space, J)

    #Apply permuation to elements
    test_tetrahedron_element  = simplex(
            test_tetrahedron_element.vertices[I[1]],
            test_tetrahedron_element.vertices[I[2]],
            test_tetrahedron_element.vertices[I[3]],
            test_tetrahedron_element.vertices[I[4]])
    trial_tetrahedron_element  = simplex(
            trial_tetrahedron_element.vertices[J[1]],
            trial_tetrahedron_element.vertices[J[2]],
            trial_tetrahedron_element.vertices[J[3]],
            trial_tetrahedron_element.vertices[J[4]])

    #Define integral (returns a function that only needs barycentric coordinates)
    if strat.sing isa SauterSchwab3D.Singularity6DFace  # cf6d    !!! Achtung nur weil es für lag:value/gradient passt heißt nicth dass vect auch passt volumen invertieren???
        
        tet = trial_tetrahedron_element
        tangs = SVector{3,SVector{3,Float64}}(-tet.tangents[1],-tet.tangents[2],-tet.tangents[3])
        vol = -tet.volume
        trial_tetrahedron_element_inv = CompScienceMeshes.Simplex(tet.vertices,tangs,tet.normals,vol) # <--- MINUS volumen???

        igd = VIEIntegrandΩΩ_cf6d(test_tetrahedron_element, trial_tetrahedron_element, trial_tetrahedron_element_inv,
        op, test_local_space, trial_local_space)
    else
       igd = VIEIntegrandΩΩ(test_tetrahedron_element, trial_tetrahedron_element,
       op, test_local_space, trial_local_space)
    end


    # if strat.sing isa SauterSchwab3D.Singularity6DEdge
    #     for j in 1 : length(L)
    #         for i  in 1 : length(K)
    #             out[i,j] += 0.0
    #         end
    #     end 
    #     return nothing
    # end


    #Evaluate integral
    Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat) 
    
    #Undo permuation on DoFs
    for j in 1 : length(L)
        for i  in 1 : length(K)
            out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
        end
    end
    nothing
end












#########################################################################

abstract type MaterialIdentity <: BEAST.LocalOperator end

struct KernelValsMaterialIdentity{U}
    tau::U # nehmen gleiche Notation wie in VIE Part
end

function BEAST.kernelvals(localop::MaterialIdentity, p)

    tau = localop.tau(cartesian(p)) #skalare oder tensorielle Funktion der Ortes

    return KernelValsMaterialIdentity(tau)
end

BEAST.scalartype(localop::MaterialIdentity) = typeof(localop.α)  # typeof(tau) geht ja schlecht weil tau function ist


















    #display(test_tetrahedron_element.vertices)
    #display(trial_tetrahedron_element.vertices)

    # function tet_circ_LHS(s) # simplex als input
    #     is_circ_lhs = false
    #     #s = simplex(mesh.vertices[face])
    #     t41=s.vertices[1]-s.vertices[4]
    #     t42=s.vertices[2]-s.vertices[4]
    #     t43=s.vertices[3]-s.vertices[4]
    
    #     @assert t41 == s.tangents[1]
    #     @assert t42 == s.tangents[2]
    #     @assert t43 == s.tangents[3]
    
    #     c = cross(t41,t42)
    #     dot(c,t43) > 0.0 && (is_circ_lhs = true)
    
    
    #     # t12=s.vertices[1]-s.vertices[2]
    #     # t13=s.vertices[1]-s.vertices[3]
    #     # t14=s.vertices[4]-s.vertices[1]
    #     # c2=cross(t12,t13)
    #     # @assert dot(c2,t14)<0
    
    
    #     """
    #     Ergebnis:
    
    #     1 2 3 CircRechteHand => n1
    #     => 4 ist in -n1  Richtung
        
    #     Alternativ direkt Linke hand
    #     """
    
    #     return is_circ_lhs
    # end

    # tet = test_tetrahedron_element
    # tri = trial_tetrahedron_element
    # @assert tet_circ_LHS(tet) == true
    # c_tri = cartesian(CompScienceMeshes.center(tri))
    # d = dot(tri.normals[1], tet.vertices[1] - c_tri) + dot(tri.normals[1], tet.vertices[2] - c_tri) + dot(tri.normals[1], tet.vertices[3] - c_tri) + dot(tri.normals[1], tet.vertices[4] - c_tri)
    # # drei sind immer null ... aber nur so bekommt man den nicht-tri knoten....
    # @assert d < 0.0







    # tet = test_tetrahedron_element
    # tri = trial_tetrahedron_element
    # @assert tet_circ_LHS(tet) == true
    # c_tri = cartesian(CompScienceMeshes.center(tri))
    # d = dot(tri.normals[1], tet.vertices[1] - c_tri) + dot(tri.normals[1], tet.vertices[2] - c_tri) + dot(tri.normals[1], tet.vertices[3] - c_tri) + dot(tri.normals[1], tet.vertices[4] - c_tri)
    # # drei sind immer null ... aber nur so bekommt man den nicht-tri knoten....
    # @assert d < 0.0


    # if typeof(strat) <: SauterSchwab3D.CommonFace5D_S
    #     @show typeof(strat)
    #     tet = test_tetrahedron_element
    #     tri = trial_tetrahedron_element # ja... eigentlich dreieck...
    #     # Das ist die Darstellung die ANGEBLICH durch die reorder Funktion erreicht wird!
    #     # @show norm(tet[1] - tri[1]) #< 1.0e-14
    #     # @show norm(tet[2] - tri[2]) #< 1.0e-14
    #     # @show norm(tet[4] - tri[3]) #< 1.0e-14


    #     # Das ist die Darstellung die laut example_cf_2.5d.jl nötig ist
    #     # const P = simplex(pI,pII,pIV,pIII)
    #     # const Q = simplex(pI,pIII,pII)
    #     @show norm(tet[1] - tri[1]) < 1.0e-14
    #     @show norm(tet[2] - tri[3]) < 1.0e-14
    #     @show norm(tet[4] - tri[2]) < 1.0e-14
        

    #     display(tet.vertices)
    #     display(tri.vertices)
    #     display("-------------------------------------")
    # end
    #@show strat.sing.T
    #@show strat.sing.S
    # function reorder(sing::Singularity5DFace)
    # Find the permutation P of t and s that make
    # Pt = [P1, P2, A1, P3]
    # Ps = [P1, P2, P3]



