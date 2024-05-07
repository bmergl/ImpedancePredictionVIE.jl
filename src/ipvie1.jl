module IPVIE1    # HAUPTMODUL: Konstruktor für Operatoren der Version 1
    using ..ImpedancePredictionVIE # <----- geht nur wenn das auch das PARENT MODUL ist
    mainMod = ImpedancePredictionVIE
    using BEAST



    # ACHTUNG! Im VIE Teil bei BEAST ist die Materialfunktion immer "tau" d.h.operatoren
    # haben die Felder gamma, alpha, tau! Hier kann trotzdem chi gefordert werden!

    # TL-Block
    
    function tl_Γ() # 3D
        
        return Identity() # Später hier - 1/2....
    end

    function tl_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)
        #alpha != 1.0 && error("Can't change tl_ΓΓ, see tl_Γ connection!")

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end

    
    # TR-Block
    
    function tr_ΓΩ(; gammatype = ComplexF64, alpha = -1.0, invtau = nothing) # 5D
        gamma = gammatype(0.0)

        return mainMod.tr_ΓΩ(gamma, alpha, invtau)
    end

    function tr_ΓΓ(; gammatype = ComplexF64, alpha = 2.0, invtau = nothing) # 4D
        gamma = gammatype(0.0)
        invtau === nothing && error("Add invtau=SCALARVALUE to tr_ΓΓ")

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha*invtau)
    end


    # BL-Block

    function bl_ΓΓ(; gammatype = ComplexF64, alpha = 1.0) # 4D
        gamma = gammatype(0.0)

        return Helmholtz3D.doublelayer(gamma = gamma, alpha = alpha) + alpha*(-1/2)*Identity()
    end

    function bl_ΩΓ(; gammatype = ComplexF64, alpha = -1.0) # 5D
        gamma = gammatype(0.0)
        tau = x -> 1.0  # VIE-kernelvals fordert tau

        return mainMod.bl_ΩΓ(gamma, alpha, tau)
    end


    # BR-Block

    function br_Ω(; alpha = -1.0, invtau = nothing) # 3D (Material Identity)

        invtau === nothing && error("invtau=FUNCTION missing in br_Ω")

        return mainMod.br_Ω(alpha, invtau)  
    end

    function br_ΓΩ(; gammatype = ComplexF64, alpha = 1.0, invtau = nothing) # 5D
        gamma = gammatype(0.0)

        return mainMod.br_ΓΩ(gamma, alpha, invtau)
    end

    function br_ΩΩ(; gammatype = ComplexF64, alpha = -1.0, invtau = nothing) # 6D
        gamma = gammatype(0.0)

        return mainMod.br_ΩΩ(gamma, alpha, invtau)
    end



    function br_ΓΓ(; gammatype = ComplexF64, alpha = 2.0, invtau = nothing) # 4D
        gamma = gammatype(0.0)
        invtau === nothing && error("Add invtau=SCALARVALUE to br_ΓΓ")

        return Helmholtz3D.singlelayer(gamma = gamma, alpha = alpha * invtau)
    end
    
    function br_ΩΓ(; gammatype = ComplexF64, alpha = -2.0, invtau = nothing) # 5D
        gamma = gammatype(0.0)
        #tau = x -> 1.0 # VIE-kernelvals fordert tau

        invtau === nothing && error("Add invtau=SCALARVALUE to br_ΩΓ")

        return mainMod.br_ΩΓ(gamma, alpha, invtau)
    end













    
    # function bl_Ω(; alpha = -1.0) # <---   -1.0 !     IDEE DENN ∇ Volumen-.----- Identity....

    #     tau = x -> 1.0

    #     return mainMod.br_Ω(alpha, tau)
    # end


end