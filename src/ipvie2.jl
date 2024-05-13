module IPVIE2    # HAUPTMODUL: Konstruktor für Operatoren der Version 2
    using ..ImpedancePredictionVIE # <----- geht nur wenn das auch das PARENT MODUL ist
    Mod = ImpedancePredictionVIE
    using BEAST


    # ACHTUNG! Im VIE Teil bei BEAST ist die Materialfunktion immer "tau" d.h.operatoren
    # haben die Felder gamma, alpha, tau! Hier kann trotzdem tau=different_tau stehen!


    # B11 Block
    function B11_Γ(;alpha = 1.0) # 3D

        return alpha*Identity() # Später hier - 1/2....
    end
    function B11_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end

    # B12 Block
    function B12_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end

    # B13 Block
    function B13_ΓΓ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) # 5D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.MaterialSL(gamma, alpha, chi)
    end
    function B13_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) # 5D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.gradG_ΓΩ(gamma, alpha, chi)
    end

    # B21 Block
    function B21_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)
 
        return Mod.HyperSingularDyadic(gamma, alpha) #Helmholtz3D.hypersingular(gamma = gamma, alpha = alpha) 
    end

    # B22 Block
    function B22_Γ(; alpha = 1.0) # 3D

        return alpha*Identity()
    end
    function B22_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer_transposed(gamma = gamma, alpha = alpha) + alpha*(1/2)*Identity()#!!!! + oder - Identity oder gar nichts???:
    end

    # B23 Block
    function B23_ΓΓ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) # 4D
        gamma = gammatype(0.0)
        chi === nothing && error("")
        # Testintegral ∫_dS ist skalar! n̂ kommt vom Operator => Identity besteht aus skalarer Testfunktion und skalarer Basis (ntrace) 
        return Mod.MaterialADL(gamma, alpha, chi) + alpha*(1/2)*Identity() #!!!!! + oder - Identity oder gar nichts???
    end
    function B23_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) # 5D
        gamma = gammatype(0.0)
        chi === nothing && error("")

        return Mod.n_dyadG_ΓΩ(gamma, alpha, chi) # + Zusatzterm??????
    end

    # B31 Block
    function B31_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B31_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals needs tau

        return Mod.div_ngradG_ΩΓ(gamma, alpha, tau) 
    end

    # B32 Block
    function B32_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function B32_ΩΓ(; gammatype = ComplexF64, alpha = 1.0) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals needs tau

        return Mod.div_G_ΩΓ(gamma, alpha, tau)
    end

    # B33 Block
    function B33_Ω(; alpha = -1.0, invtau = nothing) # 3D (Material Identity)
        invtau === nothing && error("invtau=FUNCTION missing in B33_Ω")

        return Mod.MatIdΩ(alpha, invtau)  
    end
    function B33_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.MaterialSL(gamma,alpha,chi)
    end
    function B33_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.n_gradG_ΓΩ(gamma, alpha, chi)#... FALSCH!!!??
    end
    function B33_ΩΓ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.div_G_ΩΓ(gamma, alpha, chi)#...
    end
    function B33_ΩΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return Mod.div_gradG_ΩΩ(gamma, alpha, chi)#...
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