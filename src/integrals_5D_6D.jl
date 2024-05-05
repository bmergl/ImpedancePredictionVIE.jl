
####################################################################################
# SauterSchwab3D and the corresponding BEAST.momintegrals! work for ∫∫_∂Ω ∫∫∫_Ω and ∫∫∫_Ω ∫∫∫_Ω but NOT for  ∫∫∫_Ω ∫∫_∂Ω
# This can be fixed by split BEAST.momintegrals! into three parts, one for every new operator type.  

# strat.sing.T für Testkörper, strat.sing.S für Basiskörper,
# reorder(sing::Singularity5D...) hat IMMER return TET, TRI => Bei sing. kann es vertauscht sein!
####################################################################################



#abstract type NAME <: BEAST.VolumeOperator end #sicher nötig? haben wird  doch schon!!!

#abstract type BoundaryOperatorΓΩ <: BEAST.BoundaryOperator end



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
    #error("")

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
    
    integrand(igd.op, kerneldata,tval,tgeo,bval,bgeo) * j # HIER WAR ZENTRALER FEHLER!!!!
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



    @assert length(trial_tetrahedron_element.vertices) == 3
    @assert length(test_tetrahedron_element.vertices) == 4

    #Define integral (returns a function that only needs barycentric coordinates)
    igd = VIEIntegrandΩΓ(test_tetrahedron_element, trial_tetrahedron_element,
        op, test_local_space, trial_local_space)

    
    
    #Evaluate integral
    Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat) 
    
    # Es wird erwartet f((1-x1,x1-x1*x2,x1*x2*x3),(1-x1*x4,x1*x4*x5)) für die igd function d.h. igd(3D,2D)
  
    #Undo permuation on DoFs
    for j in 1 : length(L)
        for i  in 1 : length(K)
            out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
        end
    end
    nothing
end







#Find permutation of vertices to match location of singularity to SauterSchwab
    #T = deepcopy(strat.sing.T) # length 1 to 3
    #S = deepcopy(strat.sing.S) # same length         

    #strat.sing.T = S
    #strat.sing.S = T

    # Ist das mutable?





############################### REF
# function momintegrals!(op::VIEOperator,
#     test_local_space::RefSpace, trial_local_space::RefSpace,
#     test_tetrahedron_element, trial_tetrahedron_element, out, strat::SauterSchwab3DStrategy)

#     #Find permutation of vertices to match location of singularity to SauterSchwab
#     J, I= SauterSchwab3D.reorder(strat.sing)
      
#     #Get permutation and rel. orientatio of DoFs 
#     K,O1 = reorder_dof(test_local_space, I)
#     L,O2 = reorder_dof(trial_local_space, J)
#     #Apply permuation to elements
 
#     if length(I) == 4
#         test_tetrahedron_element  = simplex(
#             test_tetrahedron_element.vertices[I[1]],
#             test_tetrahedron_element.vertices[I[2]],
#             test_tetrahedron_element.vertices[I[3]],
#             test_tetrahedron_element.vertices[I[4]])
#     elseif  length(I) == 3
#         test_tetrahedron_element  = simplex(
#             test_tetrahedron_element.vertices[I[1]],
#             test_tetrahedron_element.vertices[I[2]],
#             test_tetrahedron_element.vertices[I[3]])
#     end

#     #test_tetrahedron_element  = simplex(test_tetrahedron_element.vertices[I]...)

#     if length(J) == 4
#     trial_tetrahedron_element  = simplex(
#         trial_tetrahedron_element.vertices[J[1]],
#         trial_tetrahedron_element.vertices[J[2]],
#         trial_tetrahedron_element.vertices[J[3]],
#         trial_tetrahedron_element.vertices[J[4]])
#     elseif  length(J) == 3
#         trial_tetrahedron_element  = simplex(
#         trial_tetrahedron_element.vertices[J[1]],
#         trial_tetrahedron_element.vertices[J[2]],
#         trial_tetrahedron_element.vertices[J[3]])
#     end

#     #trial_tetrahedron_element = simplex(trial_tetrahedron_element.vertices[J]...)

#     #Define integral (returns a function that only needs barycentric coordinates)
#     igd = VIEIntegrand(test_tetrahedron_element, trial_tetrahedron_element,
#         op, test_local_space, trial_local_space)

#     #Evaluate integral
#     Q = SauterSchwab3D.sauterschwab_parameterized(igd, strat)
  
#     #Undo permuation on DoFs
#     for j in 1 : length(L)
#         for i  in 1 : length(K)
#             out[i,j] += Q[K[i],L[j]]*O1[i]*O2[j]
#         end
#     end
#     nothing
# end











