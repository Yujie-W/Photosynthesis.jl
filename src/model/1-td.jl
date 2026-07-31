"""

    temperature_correction(
                td::AbstractTemperatureDependency{FT},
                t::FT;
                t_ref::FT = td.T_REF
        ) where {FT}

Return the correction ratio for a temperature dependent variable, given
- `td` `Arrhenius`, `ArrheniusPeak`, `ArrheniusPeak2`, `Q10`, `Q10Peak`, `Q10PeakHT`, or `Q10PeakLTHT` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

"""
function temperature_correction end;

temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = exp( td.ΔHA / GAS_R(FT) * (1/t_ref - 1/t) );

temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHA, ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = exp( ΔHA / GAS_R(FT) * (1 / t_ref - 1 / t) );
    f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return f_a * f_b
);

temperature_correction(td::ArrheniusPeak2{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHA, ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = exp( ΔHA / GAS_R(FT) * (1 / t_ref - 1 / t) );
    f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return min(1,f_a) * min(1,f_b)
);

temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = td.Q_10 ^ ( (t - t_ref) / 10 );

temperature_correction(td::Q10Peak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHD, ΔSV) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = td.Q_10 ^ ( (t - t_ref) / 10 );
    f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return f_a * f_b
);

temperature_correction(td::Q10PeakHT{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔT_REF, ΔT_SLOPE) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = td.Q_10 ^ ( (t - t_ref) / 10 );
    f_b = 1 / (1 + exp(ΔT_SLOPE * (t - ΔT_REF)));

    return f_a * f_b
);

temperature_correction(td::Q10PeakLTHT{FT}, t::FT; t_ref::FT = td.T_REF) where {FT} = (
    (; ΔHT_REF, ΔHT_SLOPE, ΔLT_REF, ΔLT_SLOPE) = td;

    # f_a: activation correction, f_b: de-activation correction
    f_a = td.Q_10 ^ ( (t - t_ref) / 10 );
    f_b = 1 / (1 + exp(ΔHT_SLOPE * (t - ΔHT_REF))) / (1 + exp(ΔLT_SLOPE * (ΔLT_REF - t)));

    return f_a * f_b
);


"""

    temperature_corrected_value(td::AbstractTemperatureDependency{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}

Return the temperature corrected value, given
- `td` `Arrhenius`, `ArrheniusPeak`, `ArrheniusPeak2`, `Q10`, `Q10Peak`, `Q10PeakHT`, or `Q10PeakLTHT` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

"""
function temperature_corrected_value(td::UnionTD{FT}, t::FT; t_ref::FT = td.T_REF) where {FT}
    return td.VAL_REF * temperature_correction(td, t; t_ref=t_ref)
end;


"""

    photosystem_temperature_dependence!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                p_air::FT,
                t::FT) where {FT}

Update the temperature dependencies of C3 photosynthesis model, given
- `methods` `PhotosynthesisMethods` structure that contains rate constants, temperature dependency methods, colimitation methods, and model features
- `ps` `LeafPhotosystem` structure
- `p_air` Surface atmosphere pressure in `Pa`
- `t` Target temperature in `K`

"""
function photosystem_temperature_dependence! end;

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_air::FT,
            t::FT) where {FT} = photosystem_temperature_dependence!(methods, ps.trait, ps.auxil, p_air, t);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            t::FT) where {FT} = photosystem_temperature_dependence!(methods, pst, psa, methods.C3_AC_METHOD, methods.C3_AJ_METHOD, methods.C3_AP_METHOD, p_air, t);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            t::FT) where {FT} = photosystem_temperature_dependence!(methods, pst, psa, methods.C4_AC_METHOD, methods.C4_AJ_METHOD, methods.C4_AP_METHOD, p_air, t);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3JmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax},
            p_air::FT,
            t::FT) where {FT} = (
    (; PS_RATE_CONSTANTS, TD_JMAX, TD_KC, TD_KO, TD_R_C3, TD_VCMAX_C3, TD_Γ) = methods;

    psa.j_max  = pst.j_max25 * temperature_correction(TD_JMAX, t);
    psa.k_c    = temperature_corrected_value(TD_KC, t);
    psa.k_o    = temperature_corrected_value(TD_KO, t);
    psa.r_d    = pst.r_d25 * temperature_correction(TD_R_C3, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(TD_VCMAX_C3, t);
    psa.γ_star = temperature_corrected_value(TD_Γ, t);
    psa.k_m    = psa.k_c * (1 + p_air * F_O₂(FT) / psa.k_o);

    # TODO: make this a method?
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC3VcmaxPi,
            ajm::AjMethodC3VqmaxPi,
            apm::Union{ApMethodC3Inf, ApMethodC3Vcmax},
            p_air::FT,
            t::FT) where {FT} = (
    (; FIX_ETA_TD, PSI_RATE_CONSTANTS, TD_KC, TD_KO, TD_KQ, TD_R_C3, TD_VCMAX_C3, TD_Γ, TD_ηC, TD_ηL) = methods;

    psa.k_c    = temperature_corrected_value(TD_KC, t);
    psa.k_o    = temperature_corrected_value(TD_KO, t);
    psa.k_q    = temperature_corrected_value(TD_KQ, t);
    psa.r_d    = pst.r_d25 * temperature_correction(TD_R_C3, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(TD_VCMAX_C3, t);
    psa.γ_star = temperature_corrected_value(TD_Γ, t);
    psa.η_c    = temperature_corrected_value(TD_ηC, t);
    psa.η_l    = temperature_corrected_value(TD_ηL, t);
    psa.k_m    = psa.k_c * (1 + p_air * F_O₂(FT) / psa.k_o);
    psa.v_qmax = pst.b₆f * psa.k_q;

    if FIX_ETA_TD
        psa.η_c = min(TD_ηC.VAL_REF, psa.η_c);
        psa.η_l = min(TD_ηL.VAL_REF, psa.η_l);
    end;

    psa.ϕ_psi_max = PSI_RATE_CONSTANTS.K_P / (PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VcmaxPi,
            p_air::FT,
            t::FT) where {FT} = (
    (; PS_RATE_CONSTANTS, TD_KPEP_CLM, TD_R_C4, TD_VCMAX_C4) = methods;

    psa.k_pep_clm = temperature_corrected_value(TD_KPEP_CLM, t);
    psa.r_d       = pst.r_d25 * temperature_correction(TD_R_C4, t);
    psa.v_cmax    = pst.v_cmax25 * temperature_correction(TD_VCMAX_C4, t);

    # TODO: make this a method?
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);

photosystem_temperature_dependence!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            acm::AcMethodC4Vcmax,
            ajm::AjMethodC4JPSII,
            apm::ApMethodC4VpmaxPi,
            p_air::FT,
            t::FT) where {FT} = (
    (; PS_RATE_CONSTANTS, TD_KPEP, TD_R_C4, TD_VCMAX_C4, TD_VPMAX) = methods;

    psa.k_pep  = temperature_corrected_value(TD_KPEP, t);
    psa.r_d    = pst.r_d25 * temperature_correction(TD_R_C4, t);
    psa.v_cmax = pst.v_cmax25 * temperature_correction(TD_VCMAX_C4, t);
    psa.v_pmax = pst.v_pmax25 * temperature_correction(TD_VPMAX, t);

    # TODO: make this a method?
    psa.k_d = max(0.8738, 0.0301 * (t - 273.15) + 0.0773);
    psa.ϕ_psii_max = PS_RATE_CONSTANTS.K_P / (psa.k_d + PS_RATE_CONSTANTS.K_F + PS_RATE_CONSTANTS.K_P);

    return nothing
);
