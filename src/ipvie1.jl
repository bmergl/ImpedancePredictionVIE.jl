module IPVIE1    # HAUPTMODUL: Konstruktor für Operatoren der Version 1
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

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end

    # B13 Block
    function UB13_ΓΩ(; gammatype = ComplexF64, alpha = 1.0) # 5D
        gamma = gammatype(0.0)

        return Mod.Gdiv_ΓΩ(gamma, alpha, x->1.0)
    end
    function UB13_ΓΓn(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end


    # B21 Block
    function B21_ΓΓ(; gammatype = ComplexF64, beta = -1.0) # 4D 
        gamma = gammatype(0.0)
 
        return Helmholtz3D.hypersingular(gamma = gamma, beta = beta) # das versteckt alpha std. Null für gamma =0.0!!! !ja ...-1 für beta
    end

    # B22 Block
    function UB22_Γ(; alpha = -1.0) # 3D

        return alpha*Identity()
    end
    function UB22_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer_transposed(gamma = gamma, alpha = alpha)
    end

    # B23 Block
    function UB23_ΓΩ(; gammatype = ComplexF64, alpha = 1.0) # 5D
        gamma = gammatype(0.0)
        
        return Mod.n_gradGdiv_ΓΩ(gamma, alpha, x->1.0)
    end
    function UB23_ΓΓn(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer_transposed(gamma = gamma, alpha = alpha)
    end



    # B31 Block
    function B31_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B31_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)

        return Mod.div_ngradG_ΩΓ(gamma, alpha, x -> 1.0) 
    end

    # B32 Block
    function UB32_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function UB32_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)

        return Mod.div_G_ΩΓ(gamma, alpha, x -> 1.0)
    end

    # B33 Block
    function UB33_Ω(; alpha = -1.0)

        return alpha*Identity()
    end
    function UB33_ΓΩ(; gammatype = ComplexF64, alpha = 1.0) #
        gamma = gammatype(0.0)

        return Mod.Gdiv_ΓΩ(gamma, alpha, x -> 1.0)
    end
    function UB33_ΩΩ(; gammatype = ComplexF64, alpha = -1.0) #
        gamma = gammatype(0.0)

        return Mod.div_Gdiv_ΩΩ(gamma, alpha, x -> 1.0)
    end
    function UB33_ΓΓn(; gammatype = ComplexF64, alpha = 1.0) # 
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function UB33_ΩΓn(; gammatype = ComplexF64, alpha = -1.0) #
        gamma = gammatype(0.0)

        return Mod.div_G_ΩΓ(gamma, alpha, x -> 1.0)
    end



end