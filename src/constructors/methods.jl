"""Use the C3 Jmax method for the given PhotosynthesisMethods instance."""
function use_c3_jmax_model!(methods::PhotosynthesisMethods{FT}) where {FT}
    methods.C3_AC_METHOD = AcMethodC3Vcmax();
    methods.C3_AJ_METHOD = AjMethodC3Jmax();
    methods.C3_AP_METHOD = ApMethodC3Vcmax();
    methods.COLIMIT_J = ColimitJCLM(FT);
    methods.FLUORESCENCE_METHOD_C3 = KNFluorescenceModel{FT}();

    return nothing
end;


"""use the C3 Vqmax method for the given PhotosynthesisMethods instance."""
function use_c3_vqmax_model!(methods::PhotosynthesisMethods{FT}) where {FT}
    methods.C3_AC_METHOD = AcMethodC3Vcmax();
    methods.C3_AJ_METHOD = AjMethodC3Vqmax();
    methods.C3_AP_METHOD = ApMethodC3Vcmax();
    methods.COLIMIT_J = SerialColimit();
    methods.FLUORESCENCE_METHOD_C3 = CytochromeFluorescenceModel();

    return nothing
end;


"""Use the C4 Vcmax method for the given PhotosynthesisMethods instance."""
function use_c4_vcmax_model!(methods::PhotosynthesisMethods{FT}) where {FT}
    methods.C4_AC_METHOD = AcMethodC4Vcmax();
    methods.C4_AJ_METHOD = AjMethodC4JPSII();
    methods.C4_AP_METHOD = ApMethodC4Vcmax();
    methods.FLUORESCENCE_METHOD_C4 = KNFluorescenceModel{FT}();

    return nothing
end;


"""Use the C4 Vpmax method for the given PhotosynthesisMethods instance."""
function use_c4_vpmax_model!(methods::PhotosynthesisMethods{FT}) where {FT}
    methods.C4_AC_METHOD = AcMethodC4Vcmax();
    methods.C4_AJ_METHOD = AjMethodC4JPSII();
    methods.C4_AP_METHOD = ApMethodC4Vpmax();
    methods.FLUORESCENCE_METHOD_C4 = KNFluorescenceModel{FT}();

    return nothing
end;


"""Use the C4 Vqmax method for the given PhotosynthesisMethods instance."""
function use_c4_vqmax_model!(methods::PhotosynthesisMethods{FT}) where {FT}
    methods.C4_AC_METHOD = AcMethodC4Vqmax();
    methods.C4_AJ_METHOD = AjMethodC4Vqmax();
    methods.C4_AP_METHOD = ApMethodC4Vqmax();
    methods.COLIMIT_J = SerialColimit();
    methods.FLUORESCENCE_METHOD_C4 = CytochromeFluorescenceModel();

    return nothing
end;
