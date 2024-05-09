module IPVIE2    # HAUPTMODUL: Konstruktor für Operatoren der Version 2
    using ..ImpedancePredictionVIE # <----- geht nur wenn das auch das PARENT MODUL ist
    Mod = ImpedancePredictionVIE
    using BEAST


    # ACHTUNG! Im VIE Teil bei BEAST ist die Materialfunktion immer "tau" d.h.operatoren
    # haben die Felder gamma, alpha, tau! Hier kann trotzdem tau=different_tau stehen!


    # B11 Block
    function B11_Γ() # 3D

        return Identity() # Später hier - 1/2....
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
    # MaterialSL nötig!!! ############
    function B13_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, invtau = nothing) # 5D
        gamma = gammatype(0.0)

        return #evtl gibt es den op in beast schon !mainMod.tr_ΓΩ(gamma, alpha, invtau)
    end

    # B21 Block
    function B21_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!! + oder - Identity oder gar nichts???:
        return #Helmholtz3D. Hypersingular????? oder was anderes?
    end

    # B22 Block
    function B22_Γ() # 3D

        return Identity()
    end
    function B22_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!! + oder - Identity oder gar nichts???:
        return Helmholtz3D.doublelayer_transposed(gamma = gamma, alpha = alpha) + alpha*(1/2)*Identity()
    end

    # B23 Block
    function B23_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!! + oder - Identity oder gar nichts???:
        return # MaterialDL...Problem: Brauchen trace des Tensors... + seltsame Identity
    end
    function B23_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing) # 5D
        gamma = gammatype(0.0)

        return #Dyade... + Zusatzterm???
    end

    # B31 Block
    function B31_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end
    function B31_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals needs tau

        return # DEF THIS OF GENERAL... does not ex in BEAST ... mainMod.bl_ΩΓ(gamma, alpha, tau)
    end

    # B32 Block
    function B32_ΓΓ(; gammatype = ComplexF64, alpha = -1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha)
    end
    function B32_ΩΓ(; gammatype = ComplexF64, alpha = 1.0) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals needs tau

        return # DEF THIS OF GENERAL... does not ex in BEAST ... mainMod.bl_ΩΓ(gamma, alpha, tau)
    end

    # B33 Block
    function B33_Ω(; alpha = -1.0, invtau = nothing) # 3D (Material Identity)
        invtau === nothing && error("invtau=FUNCTION missing in B33_Ω")

        return mainMod.MatIdΩ(alpha, invtau)  
    end
    function B33_ΓΓ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return #...
    end
    function B33_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return #...
    end
    function B33_ΩΓ(; gammatype = ComplexF64, alpha = -1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return #...
    end
    function B33_ΩΩ(; gammatype = ComplexF64, alpha = 1.0, chi = nothing) #
        gamma = gammatype(0.0)
        chi === nothing && error("chi missing")

        return #...
    end




end