



###### MaterialSL ##############################################################
struct MaterialSL{T,K,Q} <: BEAST.Helmholtz3DOp{T,K}
    gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
    alpha::T
    tau::Q
end

function (igd::BEAST.Integrand{<:MaterialSL})(x,y,f,g)
    α = igd.operator.alpha
    γ = BEAST.gamma(igd.operator)

    Ty = igd.operator.tau(y) #!!!

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1 / R
    green = exp(-γ*R)*(BEAST.i4pi*iR)

    αGTy = α * green * Ty

    BEAST._integrands(f,g) do fi, gi
        dot(gi.value, αGTy*fi.value)
    end
end

###### MaterialDL ##############################################################
struct MaterialDL{T,K,Q} <: BEAST.Helmholtz3DOp{T,K}
    gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
    alpha::T
    tau::Q
end

function (igd::BEAST.Integrand{<:MaterialDL})(x,y,f,g)
    γ = BEAST.gamma(igd.operator)
    α = igd.operator.alpha

    Ty = igd.operator.tau(y) #!!!

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*BEAST.i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)
    αgradgreen = α * gradgreen
    n = normal(y)
    fvalue = BEAST.getvalue(f)
    gvalue = BEAST.getvalue(g)

    return BEAST._krondot(fvalue,gvalue) * dot(n, -αgradgreen*Ty) # "-" siehe ∇' ... 
end

###### MaterialADL ##############################################################
struct MaterialADL{T,K,Q} <: BEAST.Helmholtz3DOp{T,K}
    gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
    alpha::T
    tau::Q
end

function (igd::BEAST.Integrand{<:MaterialADL})(x,y,f,g)
    γ = BEAST.gamma(igd.operator)
    α = igd.operator.alpha

    Ty = igd.operator.tau(y) #!!!

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*BEAST.i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)
    αgradgreen = α * gradgreen
    n = normal(x)           # Unterschied zu DL ist 1. normal(x) statt y 
    fvalue = BEAST.getvalue(f)
    gvalue = BEAST.getvalue(g)

    return BEAST._krondot(fvalue,gvalue) * dot(n, αgradgreen*Ty) # 2. + statt -
end

###### Hypersingular Dyadic: KONVERGIERT  NICHT!!!! ##############################################################

# struct HyperSingularDyadic{T,K} <: BEAST.Helmholtz3DOp{T,K}
#     gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
#     alpha::T
# end

# function (igd::BEAST.Integrand{<:HyperSingularDyadic})(x,y,f,g)
#     α = igd.operator.alpha
#     γ = BEAST.gamma(igd.operator)

#     r = cartesian(x) - cartesian(y)
#     R = norm(r)

#     Rsq = R^2

#     p_ = cartesian(x)
#     q_ = cartesian(y)
#     xd = p_[1]-q_[1]
#     yd = p_[2]-q_[2]
#     zd = p_[3]-q_[3]

#     dyadgreen = (1/(4*pi*R^5)) * @SMatrix [3*xd^2-Rsq xd*yd xd*zd;
#     yd*xd 3*yd^2-Rsq yd*zd;
#     zd*xd zd*yd 3*zd^2-Rsq;
#     ]

#    # rand() < 0.00001 && display(dyadgreen)
#     @assert norm(dyadgreen-transpose(dyadgreen))<1e-13

#     αdaydG = α * dyadgreen

#     nx = x.patch.normals[1]
#     ny = y.patch.normals[1] 

#     # retunrn BEAST._integrands(f,g) do fi, gi
#     #     dot(gi.value*nx, αdaydG*(fi.value*ny))
#     # end

#     fvalue = BEAST.getvalue(f)
#     gvalue = BEAST.getvalue(g)
#     return BEAST._krondot(fvalue,gvalue) * dot(nx, αdaydG*ny) 
# end

###### Hypersingular curl=0 because of pwc ... problems??? #############################################
# struct HyperSingular0curl{T,K} <: BEAST.Helmholtz3DOp{T,K}
#     gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
#     alpha::T
# end

# function (igd::BEAST.Integrand{<:HyperSingular0curl})(x,y,f,g)
#     α = igd.operator.alpha
#     #β = igd.operator.beta
#     γ = BEAST.gamma(igd.operator)

#     r = cartesian(x) - cartesian(y)
#     R = norm(r)
#     iR = 1 / R
#     green = exp(-γ*R)*(BEAST.i4pi*iR)
#     nx = normal(x)
#     ny = normal(y)

#     BEAST._integrands(f,g) do fi, gi
#         α*dot(nx,ny)*gi.value*fi.value*green #+ β*dot(gi.curl,fi.curl)*green  # <- zero for pwc
#     end
# end



##### SauterSchwab3D 5D/6D ######################################################################

struct gradG_ΓΩ{T,U,P} <: BoundaryOperatorΓΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::gradG_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gx = @SVector[tvals[i].value for i in 1:1] #!!!!!!!!!!!!!!!!!!!!!!!!
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(gx[i] * gradG, Ty*fy[j]) for i in 1:1, j in 1:4]
end
struct n_gradG_ΓΩ{T,U,P} <: BoundaryOperatorΓΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::n_gradG_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:1] #ntrace => nur ein Beitrag: 1 PWC  
    fy = @SVector[bvals[i].value for i in 1:4] 

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * gx[i] * dot(gradG, Ty*fy[j]) for i in 1:1, j in 1:4]
end

struct n_dyadG_ΓΩ{T,U,P} <: BoundaryOperatorΓΩ
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::n_dyadG_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4] 

    dyadG = -kerneldata.dyadgreen # is ∇'∇'G

    Ty = kerneldata.tau

    α = viop.α

    nx = tgeo.patch.normals[1]
    
    return @SMatrix[α * gx[i] * dot(nx, dyadG * (Ty*fy[j])) for i in 1:3, j in 1:4]
end
struct KernelValsVIEdyad{T,U,P,Q,K} # ÄNDERN!!!!!!
    gamma::U
    vect::P
    dist::T
    dyadgreen::Q
    tau::K
end
function BEAST.kernelvals(viop::n_dyadG_ΓΩ, p ,q) # p=r_vec, q=r'_vec
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

    dyadgreen = (1/(4*pi*R^5)) * @SMatrix [3*xd^2-Rsq xd*yd xd*zd;
    yd*xd 3*yd^2-Rsq yd*zd;
    zd*xd zd*yd 3*zd^2-Rsq;
    ]

    tau = viop.tau(cartesian(q))

    KernelValsVIEdyad(Y,r,R, dyadgreen, tau)
end


struct div_ngradG_ΩΓ{T,U,P} <: BoundaryOperatorΩΓ # 5D
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::div_ngradG_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)
    @assert length(bgeo.patch.vertices) == 3

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:3]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α
    

    # pnt = cartesian(bgeo) #Ortsvektor zu einem Punkt der Dreiecksfläche ausgegen vom Zentrum des Würfels => n * Ortsvektor > 0 IMMER
    # n = bgeo.patch.normals[1]
    # dot(pnt,n) < 0.0 && error("n̂ points in wrong direction!")

    return @SMatrix[α * dgx[i] * dot(bgeo.patch.normals[1], gradG * Ty * fy[j]) for i in 1:4, j in 1:3]
end

struct div_G_ΩΓ{T,U,P} <: BoundaryOperatorΩΓ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::div_G_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)
    @assert length(bgeo.patch.vertices) == 3

    dgx = @SVector[tvals[i].divergence for i in 1:4]  
    fy = @SVector[bvals[i].value for i in 1:1]  #ntrace => nur ein Beitrag: 1 PWC  
    G = kerneldata.green

    Ty = kerneldata.tau # KEIN TENSOR ERLAUBT!

    α = viop.α

    return @SMatrix[α * dgx[i] * G * Ty * fy[j] for i in 1:4, j in 1:1]
end


struct div_gradG_ΩΩ{T,U,P} <: VolumeOperatorΩΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::div_gradG_ΩΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dgx[i] * dot(gradG, Ty * fy[j]) for i in 1:4, j in 1:4]
end



# DYADE...








struct MatId{T,U} <: MaterialIdentity 
    α::T
    tau::U
end

function BEAST.integrand(localop::MatId, kerneldata, x, g, f)

    gx = g.value
    fx = f.value

    Tx = kerneldata.tau

    α = localop.α

    return α * dot(gx, Tx * fx) # geht auch wenn g und f skalar sind
end



