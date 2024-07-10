module IPVIE2 # 2 × 2 HC Formulation
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
    function UB12_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        #return Mod.MaterialSL(gamma,alpha,x -> 1.0)
        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function UB12_ΓΩ(; gammatype = ComplexF64, alpha = 1.0) # 5D
        gamma = gammatype(0.0)

        return Mod.Gdiv_Γ3Ω(gamma, alpha, x -> 1.0) 
    end
    function UB12_ΓΓn(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        #return Mod.MaterialSL(gamma,alpha,x -> 1.0)
        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end

    # B21
    function B21_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B21_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)

        return Mod.div_ngradG_ΩΓ(gamma, alpha, x -> 1.0)
    end

    # B22
    function UB22_Ω(; alpha = -1.0)

        return alpha*Identity()
    end
    function UB22_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 
        gamma = gammatype(0.0)

        #return Mod.MaterialSL(gamma,alpha,x -> 1.0)
        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function UB22_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) #
        gamma = gammatype(0.0)

        return Mod.div_G_ΩΓ(gamma, alpha, x -> 1.0)
    end
    function UB22_ΓΩ(; gammatype = ComplexF64, alpha = 1.0) #
        gamma = gammatype(0.0)

        return Mod.Gdiv_Γ1Ω(gamma, alpha, x -> 1.0)
    end
    function UB22_ΩΩ(; gammatype = ComplexF64, alpha = -1.0) #
        gamma = gammatype(0.0)

        return Mod.div_Gdiv_ΩΩ(gamma, alpha, x -> 1.0)
    end
    function UB22_ΓΓn(; gammatype = ComplexF64, alpha = 1.0) # 
        gamma = gammatype(0.0)

        #return Mod.MaterialSL(gamma,alpha,x -> 1.0)
        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function UB22_ΩΓn(; gammatype = ComplexF64, alpha = -1.0) #
        gamma = gammatype(0.0)

        return Mod.div_G_ΩΓ(gamma, alpha, x -> 1.0)
    end

    

end