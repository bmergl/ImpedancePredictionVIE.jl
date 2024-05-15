
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


##### momintegrals!  ######################################################################

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
    error("")

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

# function qr_boundary(op::BoundaryOperatorΩΓ, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
#    qs::SauterSchwab3DQStrat) no changes => take existing op::BoundaryOperator version

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
function qr_boundary(op::BoundaryOperatorΓΩ, g::RefSpace, f::RefSpace, i, τ, j,  σ, qd,
    qs::BEAST.SauterSchwab3DQStrat)  # T <-> S

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


    return DoubleQuadRule(
        qd[1][1,i],
        qd[2][1,j])

end

function BEAST.momintegrals!(op::BoundaryOperatorΓΩ,
    test_local_space::RefSpace, trial_local_space::RefSpace,
    test_tetrahedron_element, trial_tetrahedron_element, out, strat::SauterSchwab3DStrategy)

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
    error("")

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
    igd = VIEIntegrandΩΩ(test_tetrahedron_element, trial_tetrahedron_element,
        op, test_local_space, trial_local_space) # should work...

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




