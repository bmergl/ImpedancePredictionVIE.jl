#das ist  per include auf top level => using BEAST nicht nötig
#scalartype(op::VIEOperator) = typeof(op.gamma) In Beast d.h. (ε,κ) -> gamma = 0.0 + 0.0im wenn nur eines geg dann gamma = 0.0




##### SauterSchwab3D 5D/6D ######################################################################

# top right (tr) block ##########################################################################

struct tr_ΓΩ{T,U,P} <: BEAST.BoundaryOperator
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::tr_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(gx[i] * gradG, Ty*fy[j]) for i in 1:3, j in 1:4]
end


# bottom left (bl) block ##########################################################################

struct bl_ΩΓ{T,U,P} <: BoundaryOperatorΩΓ # 5D
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::bl_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)
    @assert length(bgeo.patch.vertices) == 3

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:3]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    α = viop.α

    # pnt = cartesian(bgeo) #Ortsvektor zu einem Punkt der Dreiecksfläche ausgegen vom Zentrum des Würfels => n * Ortsvektor > 0 IMMER
    # n = bgeo.patch.normals[1]
    # dot(pnt,n) < 0.0 && error("n̂ points in wrong direction!")

    return @SMatrix[α * dgx[i] * dot(bgeo.patch.normals[1], gradG * fy[j]) for i in 1:4, j in 1:3]
end


# bottom right (br) block #########################################################################

struct br_ΩΓ{T,U,P} <: BoundaryOperatorΩΓ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::br_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)
    @assert length(bgeo.patch.vertices) == 3

    dgx = @SVector[tvals[i].divergence for i in 1:4]  
    fy = @SVector[bvals[i].value for i in 1:1]  #ntrace => nur ein Beitrag: 1 PWC  
    G = kerneldata.green

    α = viop.α

    return @SMatrix[α * dgx[i] * G * fy[j] for i in 1:4, j in 1:1]
end

struct br_ΓΩ{T,U,P} <: BEAST.BoundaryOperator
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::br_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:1] #ntrace => nur ein Beitrag: 1 PWC  
    fy = @SVector[bvals[i].value for i in 1:4] 

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * gx[i] * dot(gradG, Ty*fy[j]) for i in 1:1, j in 1:4]
end

struct br_ΩΩ{T,U,P} <: BEAST.VolumeOperator
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::br_ΩΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dgx[i] * dot(gradG, Ty * fy[j]) for i in 1:4, j in 1:4]
end

###################################################################################################






##### ∫_Γ and ∫_Ω -> LocalOperator ################################################################


struct br_Ω{T,U} <: MaterialIdentity 
    α::T
    tau::U
end
function BEAST.integrand(localop::MaterialIdentity, kerneldata, x, g, f) # ist ja für beide gleich...

    gx = g.value
    fx = f.value

    Tx = kerneldata.tau


    α = localop.α

    #Tx = kerneldata.tau(cartesian(x)) # anderes als bei VIE... Denkfehler!
    #@show x
    # @show α
    #@show Tx
    # @show gx
    #@show fx
    # @show Tx * fx
    #error()

    return α * dot(gx, Tx * fx) # geht auch wenn g und f skalar sind
end






