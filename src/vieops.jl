

###### MaterialSL ##############################################################
struct MaterialSL{T,K,Q} <: BEAST.Helmholtz3DOp{T,K}
    gamma::K # Reihenfolge im VIE Teil ist gamma,alpha -> im HH3D Teil alpha,gamma
    alpha::T
    tau::Q
end

function (igd::BEAST.Integrand{<:MaterialSL})(x,y,f,g)
    α = igd.operator.alpha
    γ = BEAST.gamma(igd.operator)

    Ty = igd.operator.tau(cartesian(y)) #!!!

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

    Ty = igd.operator.tau(cartesian(y)) #!!!

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

    Ty = igd.operator.tau(cartesian(y)) #!!!

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = exp(-γ*R)*(iR*BEAST.i4pi)
    gradgreen = -(γ + iR) * green * (iR * r)
    αgradgreen = α * gradgreen
    n = normal(x)           # Unterschied zu DL ist 1. normal(x) statt y 
    fvalue = BEAST.getvalue(f)
    gvalue = BEAST.getvalue(g)


    #αgradgreen = α * 4*R^2*(cartesian(x) - cartesian(y)) # <---------------------------------------------


    return BEAST._krondot(fvalue,gvalue) * dot(n, αgradgreen*Ty) # 2. + statt -
end




##### SauterSchwab3D 5D/6D ######################################################################
struct n_gradGdiv_ΓΩ{T,U,P} <: BEAST.BoundaryOperator#ΓΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::n_gradGdiv_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gx = @SVector[tvals[i].value for i in 1:3]
    dfy = @SVector[bvals[i].divergence for i in 1:4]


    G = kerneldata.green
    gradG = kerneldata.gradgreen # "+" to get ∇G(r,r') (without ')

    Ty = kerneldata.tau

    α = viop.α

    nx = tgeo.patch.normals[1]

    return @SMatrix[α * dot(gx[i] * nx,  gradG * Ty * dfy[j]) for i in 1:3, j in 1:4]
end


struct gradG_ΓΩ{T,U,P} <: BEAST.BoundaryOperator
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::gradG_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')
    Ty = kerneldata.tau
    α = viop.α

    #αgradGTy = α * gradG * Ty # VORSICHT! Wenn Ty ComplexF64 darf es NICHT im ersten Argument des Skalarprodukt stehen!!!
    αgradG = α * gradG

    return BEAST._integrands(tvals,bvals) do fi, gj # f ist jetzt test, g ist jetzt basis!
       dot(fi.value*αgradG, Ty*gj.value)
    end
end

struct gradG_Γ1Ω{T,U,P} <: BEAST.BoundaryOperator#ΓΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::gradG_Γ1Ω, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gx = @SVector[tvals[i].value for i in 1:1] #!!!!!!!!!!!!!!!!!!!!!!!!
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(gx[i] * gradG, Ty*fy[j]) for i in 1:1, j in 1:4]
end
struct gradG_Γ3Ω{T,U,P} <: BEAST.BoundaryOperator#ΓΩ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::gradG_Γ3Ω, kerneldata, tvals, tgeo, bvals, bgeo) # 5D

    gx = @SVector[tvals[i].value for i in 1:3] #!!!!!!!!!!!!!!!!!!!!!!!!
    fy = @SVector[bvals[i].value for i in 1:4]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(gx[i] * gradG, Ty*fy[j]) for i in 1:3, j in 1:4]
end

struct n_gradG_ΓΩ{T,U,P} <: BEAST.BoundaryOperator#ΓΩ
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

struct n_dyadG_ΓΩ{T,U,P} <: BEAST.BoundaryOperator#ΓΩ
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::n_dyadG_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]
    fy = @SVector[bvals[i].value for i in 1:4] 

    dyadG = kerneldata.dyadgreen # is ∇'∇'G 

    Ty = kerneldata.tau

    α = viop.α

    nx = tgeo.patch.normals[1]
    
    return @SMatrix[α * gx[i] * dot(nx, dyadG * (Ty*fy[j])) for i in 1:3, j in 1:4]
end
BEAST.kernelvals(viop::n_dyadG_ΓΩ, p ,q) = kernelvalsdyad(viop, p, q)

struct div_ngradG_ΩΓ{T,U,P} <: BEAST.BoundaryOperator#ΩΓ # 5D
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::div_ngradG_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    fy = @SVector[bvals[i].value for i in 1:3]

    G = kerneldata.green
    gradG = -kerneldata.gradgreen # "-" to get nabla'G(r,r')

    Ty = kerneldata.tau

    α = viop.α
    
    ny = bgeo.patch.normals[1]

    gradGTy = gradG * Ty

    return @SMatrix[α * dgx[i] * dot(ny, gradGTy * fy[j]) for i in 1:4, j in 1:3]
end

struct div_G_ΩΓ{T,U,P} <: BEAST.BoundaryOperator#ΩΓ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::div_G_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:4]  
    fy = @SVector[bvals[i].value for i in 1:1] 

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dgx[i] * G * Ty * fy[j] for i in 1:4, j in 1:1]
end


struct div_gradG_ΩΩ{T,U,P} <: BEAST.VolumeOperator#ΩΩ
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


struct gradG_ΩΓ{T,U,P} <: BEAST.BoundaryOperator#ΩΓ
    gamma::T
    α::U
    tau::P
end
function BEAST.integrand(viop::gradG_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:4]  
    fy = @SVector[bvals[i].value for i in 1:1] 

    gradG = kerneldata.gradgreen # "+" to get ∇G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(gx[i], gradG * Ty * fy[j]) for i in 1:4, j in 1:1]
end



struct ncgrad_gradGc_ΓΩ{T,U,P} <: BEAST.BoundaryOperator#ΩΓ
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::ncgrad_gradGc_ΓΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].curl for i in 1:3] # = n̂ × ∇a(r_vec), a(r_vec) is 2D lin. Lagrange
    fy = @SVector[bvals[i].value for i in 1:4] 

    
    gradG = -kerneldata.gradgreen # "-" to get ∇'G(r,r')

    Ty = kerneldata.tau

    α = viop.α

    return @SMatrix[α * dot(-gx[i], cross(gradG, Ty*fy[j])) for i in 1:3, j in 1:4]

end

struct Gdiv_Γ1Ω{T,U,P} <: BEAST.BoundaryOperator
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::Gdiv_Γ1Ω, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:1]   
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α

    αGTy = α*G*Ty

    return @SMatrix[gx[i]*αGTy*dfy[j] for i in 1:1, j in 1:4]

end
struct Gdiv_Γ3Ω{T,U,P} <: BEAST.BoundaryOperator
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::Gdiv_Γ3Ω, kerneldata, tvals, tgeo, bvals, bgeo)

    gx = @SVector[tvals[i].value for i in 1:3]   
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α

    αGTy = α*G*Ty

    return @SMatrix[gx[i]*αGTy*dfy[j] for i in 1:3, j in 1:4]

end

struct div_Gdiv_ΩΩ{T,U,P} <: BEAST.VolumeOperator
    gamma::T
    α::U 
    tau::P
end
function BEAST.integrand(viop::div_Gdiv_ΩΩ, kerneldata, tvals, tgeo, bvals, bgeo)

    dgx = @SVector[tvals[i].divergence for i in 1:4]
    dfy = @SVector[bvals[i].divergence for i in 1:4]

    G = kerneldata.green

    Ty = kerneldata.tau

    α = viop.α

    αGTy = α*G*Ty

    return @SMatrix[dgx[i]*αGTy*dfy[j] for i in 1:4, j in 1:4]

end







#### Material Local Operators #######################################################

abstract type MaterialLocalOp <: BEAST.LocalOperator end

struct KernelValsMaterialLocalOp{U}
    tau::U # nehmen gleiche Notation wie in VIE Part
end
function BEAST.kernelvals(localop::MaterialLocalOp, p)

    tau = localop.tau(cartesian(p)) #skalare oder tensorielle Funktion der Ortes

    return KernelValsMaterialLocalOp(tau)
end
BEAST.scalartype(localop::MaterialLocalOp) = typeof(localop.α)  # typeof(tau) geht ja schlecht weil tau function ist



struct MatId{T,U} <: MaterialLocalOp # ja...muss ja keine identity mehr sein...besser unten....
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

struct MatLoc{T,U} <: MaterialLocalOp
    α::T
    tau::U
end
function BEAST.integrand(localop::MatLoc, kerneldata, x, g, f)

    gx = g.value
    fx = f.value

    Tx = kerneldata.tau

    α = localop.α

    return α * dot(gx, Tx * fx) # geht auch wenn g und f skalar sind
end

struct grad__Ω{T,U} <: MaterialLocalOp
    α::T
    tau::U
end
function BEAST.integrand(localop::grad__Ω, kerneldata, x, g, f)

    dgx = g.gradient
    fx = f.value


    Tx = kerneldata.tau

    α = localop.α

    return α * dot(dgx, Tx * fx)
end
struct grad_grad_Ω{T,U} <: MaterialLocalOp
    α::T
    tau::U
end
function BEAST.integrand(localop::grad_grad_Ω, kerneldata, x, g, f)

    dgx = g.gradient
    dfx = f.gradient


    Tx = kerneldata.tau

    α = localop.α

    return α * dot(dgx, Tx * dfx)
end

struct _grad_Ω{T,U} <: MaterialLocalOp
    α::T
    tau::U
end
function BEAST.integrand(localop::_grad_Ω, kerneldata, x, g, f)

    gx = g.value
    dfx = f.gradient


    Tx = kerneldata.tau

    α = localop.α

    return α * dot(gx, Tx * dfx)
end

struct _div_Ω{T,U} <: MaterialLocalOp
    α::T
    tau::U
end
function BEAST.integrand(localop::_div_Ω, kerneldata, x, g, f)

    gx = g.value
    dfx = f.divergence

    Tx = kerneldata.tau

    α = localop.α

    return α * gx * Tx * dfx
end

BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.NDLCDRefSpace{T}, ::BEAST.LagrangeRefSpace{T,D1,4}) where {T,D1} = BEAST.SingleNumQStrat(6)
function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.NDLCDRefSpace{T},
    f::BEAST.LagrangeRefSpace{T,Deg,4}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
    # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end

BEAST.defaultquadstrat(::BEAST.LocalOperator, ::BEAST.LagrangeRefSpace{T,D1,4}, ::BEAST.NDLCDRefSpace{T}) where {T,D1} = BEAST.SingleNumQStrat(6)
function BEAST.quaddata(op::BEAST.LocalOperator, g::BEAST.LagrangeRefSpace{T,Deg,4},
    f::BEAST.NDLCDRefSpace{T}, tels::Vector, bels::Vector, qs::BEAST.SingleNumQStrat) where {T,Deg} 
    # besser: quaddata(op::LocalOperator, g::LinearRefSpaceTetr, f::LinearRefSpaceTetr... gibt es aber nicht mit lag... 

    o, x, y, z = CompScienceMeshes.euclidianbasis(3)
    reftet = simplex(x,y,z,o)
    qps = quadpoints(reftet, qs.quad_rule)
    qd = [(w, parametric(p)) for (p,w) in qps]
    A = BEAST._alloc_workspace(qd, g, f, tels, bels)
    return qd, A
end




# J_ana ....
mutable struct Jana{P} <: BEAST.Functional
    J::P #function
end

BEAST.scalartype(x::Jana{P}) where {P} = Float64

function setup_Jana(; J = error("J=...as Function...")) #generator
    P = Function
    return Jana{P}(J)
end

function (e::Jana)(p) # p is CSM MP
    return e.J(BEAST.cartesian(p))
end

#*(a::Number, e::LinearPotentialVIE) = LinearPotentialVIE_(e.direction, a*e.amplitude)

BEAST.integrand(::Jana, test_vals, field_val) = test_vals[1] ⋅ field_val








#####  kernelvals für Dyadische VIE Operatoren  ######################################################################
# struct KernelValsVIEdyad{T,U,P,Q,K} # ÄNDERN!!!!!!
#     gamma::U
#     vect::P
#     dist::T
#     dyadgreen::Q
#     tau::K
# end

# function kernelvalsdyad(viop, p, q) # p=r_vec, q=r'_vec
#     # Achtung! Speziell auf gamma=0 zugeschnitten, gamma für typ wichtig
#     Y = viop.gamma  #ComplexF64/Float64 unterscheidung läuft normalerweise über Gamma...
#     r = cartesian(p)-cartesian(q)
#     R = norm(r)
#     Rsq = R^2

#     p_ = cartesian(p)
#     q_ = cartesian(q)
#     xd = p_[1]-q_[1]
#     yd = p_[2]-q_[2]
#     zd = p_[3]-q_[3]
    
#     xd_yd_3 = 3*xd*yd
#     xd_zd_3 = 3*xd*zd
#     yd_zd_3 = 3*yd*zd

#     dyadgreen =  @SMatrix [3*xd^2-Rsq   xd_yd_3      xd_zd_3;  
#                             xd_yd_3     3*yd^2-Rsq   yd_zd_3;
#                             xd_zd_3     yd_zd_3      3*zd^2-Rsq]
#     dyadgreen = dyadgreen/(4*pi*R^5) # = ∇'∇'G

#     #rand()<0.0001 && @show norm(dyadgreen)

#     # dyadgreen =  @SMatrix [8*xd^2+4*Rsq   8*xd*yd        8*xd*zd;
#     #                        8*xd*yd        8*yd^2+4*Rsq   8*yd*zd;
#     #                        8*xd*zd        8*yd*zd        8*zd^2+4*Rsq]

#     tau = viop.tau(q_)

#     KernelValsVIEdyad(Y,r,R, dyadgreen, tau)
# end


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


##### SauterSchwab3D 5D/6D dyad (problems...) ##########################################################
# struct ndyadG_ΩΓ{T,U,P} <: BEAST.BoundaryOperator#ΩΓ
#     gamma::T
#     α::U 
#     tau::P
# end
# function BEAST.integrand(viop::ndyadG_ΩΓ, kerneldata, tvals, tgeo, bvals, bgeo)

#     gx = @SVector[tvals[i].value for i in 1:4]
#     fy = @SVector[bvals[i].value for i in 1:3] 

#     dyadG = -kerneldata.dyadgreen # is ∇∇'G

#     Ty = kerneldata.tau     # ja hier ist tau eh 1.0 ....
#     Ty != 1.0 && error("Integrand ndyadG_ΩΓ not ready for tau different from 1.0")

#     α = viop.α

#     ny = bgeo.patch.normals[1]

#     dyadGny = dyadG * ny #GEHT NICHT!!!!!!!!!!!

#     # r1 = dyadG[1,1]*ny[1] + dyadG[1,2]*ny[2] + dyadG[1,3]*ny[3]
#     # r2 = dyadG[2,1]*ny[1] + dyadG[2,2]*ny[2] + dyadG[2,3]*ny[3]
#     # r3 = dyadG[3,1]*ny[1] + dyadG[3,2]*ny[2] + dyadG[3,3]*ny[3]
#     # dyadGny = @SVector [r1, r2, r3]

#     #error("STOP")
    
#     return @SMatrix[α * dot(gx[i], dyadGny) * fy[j] for i in 1:4, j in 1:3]
# end
# BEAST.kernelvals(viop::ndyadG_ΩΓ, p ,q) = kernelvalsdyad(viop, p, q) # !!!!!!!!!!!!!!! Entscheidend!!!!
