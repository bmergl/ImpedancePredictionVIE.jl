module IPVIE3 # 2 × 2 arbitrary Material Formulation (non HC)
    using ..ImpedancePredictionVIE # <----- geht nur wenn das auch das PARENT MODUL ist
    Mod = ImpedancePredictionVIE
    using BEAST


    # ACHTUNG! Im VIE Teil bei BEAST ist die Materialfunktion immer "tau" d.h.operatoren
    # haben die Felder gamma, alpha, tau! Hier kann trotzdem tau=different_tau stehen!


    # B11 Block
    function B11_Γ(; alpha = 1.0) # 3D

        return alpha*Identity() # Später hier - 1/2....
    end
    function B11_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end

    # B12 Block
    function B12_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, invtau0 = nothing) # 4D
        gamma = gammatype(0.0)
        invtau0 === nothing && error()

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha*invtau0)
    end
    function B12_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) # 5D
        gamma = gammatype(0.0)
        chi === nothing && error()

        return Mod.gradG_Γ3Ω(gamma, alpha, chi) 
    end

    # B21
    function B21_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B21_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)

        return Mod.div_ngradG_ΩΓ(gamma, alpha, x -> 1.0)
    end

    # B22
    function B22_Ω(; alpha = -1.0, invtau = nothing)
        invtau === nothing && error()

        return Mod.MatId(alpha, invtau)
    end

    function B22_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, invtau0 = nothing) # 
        gamma = gammatype(0.0)
        invtau0 === nothing && error()

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha*invtau0)
    end
    function B22_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error()

        return Mod.gradG_Γ1Ω(gamma, alpha, chi) 
    end
    function B22_ΩΓ(; gammatype = ComplexF64, alpha = -1.0, invtau0 = nothing) #
        gamma = gammatype(0.0)
        invtau0 === nothing && error()

        return Mod.div_G_ΩΓ(gamma, alpha, x -> invtau0)
    end
    function B22_ΩΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error()

        return Mod.div_gradG_ΩΩ(gamma, alpha, chi)
    end

 

    

end