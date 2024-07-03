module IPVIE    # HAUPTMODUL: Konstruktor für Operatoren der Version 2
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
    function B12_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing) # 4D
        gamma = gammatype(0.0)
        invtau === nothing && error("")

        return Mod.MaterialSL(gamma,alpha,invtau)
    end

    # B13 Block
    function B13_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 4D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.MaterialSL(gamma, alpha, chi)
    end
    function B13_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 5D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.gradG_ΓΩ(gamma, alpha, chi)
    end

    # B21 Block
    function B21_ΓΓ(; gammatype = ComplexF64, beta = 1.0) # 4D
        gamma = gammatype(0.0)
 
        return Helmholtz3D.hypersingular(gamma = gamma, beta=beta) # das versteckt alpha std. Null für gamma =0.0!!!
    end

    # B22 Block
    function B22_Γ(; alpha = 1.0, invtau = nothing) # 3D
        invtau === nothing && error("")

        return Mod.MatId(alpha, invtau)
    end
    function B22_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing) # 4D
        gamma = gammatype(0.0)
        invtau === nothing && error("")

        return Mod.MaterialADL(gamma, alpha, invtau) #+ (1/2)*Mod.MatId(alpha, invtau) # hinterer 0 immer STOPP evtl. falsch! ganz weg...
    end

    # B23 Block
    function B23_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 4D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.MaterialADL(gamma, alpha, chi) #+ (1/2)*Mod.MatId(alpha, chi) # hinterer 0 immer STOPP evtl. falsch! ganz weg...
    end
    function B23_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 5D
        
        gamma = gammatype(0.0)
        chi === nothing && error("")
        
        return Mod.ncgrad_gradGc_ΓΩ(gamma, alpha, chi)
    end

    function B23_constmed(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 5D
        @warn "χ const. needed for B23_constmed"
        gamma = gammatype(0.0)
        chi === nothing && error("")
        
        return Mod.n_gradGdiv_ΓΩ(gamma, alpha, chi)
    end

    function B23_dyad(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 5D
        @warn "B23_ΓΓ needed for B23_dyad"
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.n_dyadG_ΓΩ(gamma, alpha, chi)
    end




    # B31 Block
    function B31_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, dyad = false) # 4D
        gamma = gammatype(0.0)

        if dyad
            return 0.0*Identity()
        end

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B31_ΩΓ(; gammatype = ComplexF64, alpha = 1.0, dyad = false) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals needs tau

        if dyad
            alpha = -1.0*alpha
            return Mod.ndyadG_ΩΓ(gamma, alpha, tau)
        end

        return Mod.div_ngradG_ΩΓ(gamma, alpha, tau) 
    end

    # B32 Block
    function B32_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing, dyad = false) # 4D
        gamma = gammatype(0.0)
        invtau === nothing && error("")

        if dyad
            return 0.0*Identity()
        end

        return Mod.MaterialSL(gamma, alpha, invtau)
    end
    function B32_ΩΓ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing, dyad = false) # 5D
        gamma = gammatype(0.0)
        invtau === nothing && error("")
        
        if dyad
            alpha = -1.0*alpha
            return Mod.gradG_ΩΓ(gamma, alpha, invtau)
        end

        return Mod.div_G_ΩΓ(gamma, alpha, invtau)
    end

    # B33 Block
    function B33_Ω(; alpha = 1.0, invtau = nothing) # 3D (Material Identity)
        invtau === nothing && error("invtau=FUNCTION missing in B33_Ω")

        return Mod.MatId(alpha, invtau)
    end
    function B33_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) 
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.MaterialSL(gamma,alpha,chi)
    end
    function B33_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.n_gradG_ΓΩ(gamma, alpha, chi)
    end
    function B33_ΩΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.div_G_ΩΓ(gamma, alpha, chi)
    end
    function B33_ΩΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.div_gradG_ΩΩ(gamma, alpha, chi)
    end
    




    # Testops
    function genMatSL(;gamma,alpha = 1.0,tau = nothing)
        tau === nothing && error()
        return Mod.MaterialSL(gamma,alpha,tau)
    end
    function genMatDL(;gamma,alpha = 1.0,tau = nothing)
        tau === nothing && error()
        return Mod.MaterialDL(gamma,alpha,tau)
    end
    function genMatADL(;gamma,alpha = 1.0,tau = nothing)
        tau === nothing && error()
        return Mod.MaterialADL(gamma,alpha,tau)
    end


end