"""
Method configuration for the Photosynthesis model
"""
Base.@kwdef mutable struct PhotosynthesisMethods{FT<:AbstractFloat}
    # constants
    # photosynthesis rate constants
    "Rate constants for PSI and PSII combined (most for PSII?)"
    PS_RATE_CONSTANTS::PSRateConstants{FT} = PSRateConstants{FT}()
    "Rate constants for PSI"
    PSI_RATE_CONSTANTS::PSIRateConstants{FT} = PSIRateConstants{FT}()
    "Rate constants for PSII"
    PSII_RATE_CONSTANTS::PSIIRateConstants{FT} = PSIIRateConstants{FT}()

    # methods
    # photosynthesis models - Ac, Aj, Ap methods
    "C3 Ac method"
    C3_AC_METHOD::AcMethodC3Vcmax = AcMethodC3Vcmax()
    "C3 Aj method"
    C3_AJ_METHOD::Union{AjMethodC3Jmax,AjMethodC3Vqmax} = AjMethodC3Jmax()
    "C3 Ap method"
    C3_AP_METHOD::Union{ApMethodC3Inf,ApMethodC3Vcmax} = ApMethodC3Vcmax()
    "C4 Ac method"
    C4_AC_METHOD::Union{AcMethodC4Vcmax,AcMethodC4Vqmax} = AcMethodC4Vcmax()
    "C4 Aj method"
    C4_AJ_METHOD::Union{AjMethodC4JPSII,AjMethodC4Vqmax} = AjMethodC4JPSII()
    "C4 Ap method"
    C4_AP_METHOD::Union{ApMethodC4Vcmax,ApMethodC4Vpmax,ApMethodC4Vqmax} = ApMethodC4Vcmax()

    # photosynthesis models - colimitation methods
    "Colimitation method for Ac and Aj => Ai"
    COLIMIT_CJ::UnionColimit{FT} = MinimumColimit()
    "Colimitation method for Ai and Ap => Ag"
    COLIMIT_IP::UnionColimit{FT} = MinimumColimit()
    "Colimitation method for J (for C3 only)"
    COLIMIT_J::UnionColimit{FT} = ColimitJCLM(FT)

    # fluorescence methods
    "Fluorescence method for C3 plants"
    FLUORESCENCE_METHOD_C3::UnionFluorescence{FT} = KNFluorescenceModel{FT}()
    "Fluorescence method for C4 plants"
    FLUORESCENCE_METHOD_C4::UnionFluorescence{FT} = KNFluorescenceModel{FT}()

    # photosynthesis models - temperature dependency methods
    "Jmax temperature dependency"
    TD_JMAX::UnionTD{FT} = JmaxTDCLM(FT)
    "Kc temperature dependency"
    TD_KC::UnionTD{FT} = KcTDCLM(FT)
    "Ko temperature dependency"
    TD_KO::UnionTD{FT} = KoTDCLM(FT)
    "Kpep temperature dependency"
    TD_KPEP::UnionTD{FT} = KpepTDBoyd(FT)
    "Kpep temperature dependency to use with C4CLM method"
    TD_KPEP_CLM::UnionTD{FT} = Q10TDKpepCLM(FT)
    "Kq temperature dependency"
    TD_KQ::UnionTD{FT} = KqTDJohnson(FT)
    "Respiration temperature dependency for C3 plants"
    TD_R_C3::UnionTD{FT} = RespirationTDCLMC3(FT)
    "Respiration temperature dependency for C4 plants"
    TD_R_C4::UnionTD{FT} = RespirationTDCLMC4(FT)
    "Vcmax temperature dependency for C3 plants"
    TD_VCMAX_C3::UnionTD{FT} = VcmaxTDCLMC3(FT)
    "Vcmax temperature dependency for C4 plants"
    TD_VCMAX_C4::UnionTD{FT} = VcmaxTDCLMC4(FT)
    "Vpmax temperature dependency"
    TD_VPMAX::UnionTD{FT} = VpmaxTDBoyd(FT)
    "Γ* temperature dependency"
    TD_Γ::UnionTD{FT} = ΓStarTDCLM(FT)
    "η_C temperature dependency"
    TD_ηC::UnionTD{FT} = ηCTDWang(FT)
    "η_L temperature dependency"
    TD_ηL::UnionTD{FT} = ηLTDWang(FT)

    # features
    # optional features
    "Fix the TD of η for Johnson-Berry model"
    FIX_ETA_TD::Bool = true
end;
