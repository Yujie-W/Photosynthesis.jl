"""

    photosystem_electron_transport!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                ppar::Vector{FT},
                p_i::Union{FT,Vector{FT}};
                β::FT = FT(1)
        ) where {FT}

Update the electron transport rates, given
- `methods` `PhotosynthesisMethods` structure that contains rate constants, temperature dependency methods, colimitation methods, and model features
- `ps` `LeafPhotosystem` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `p_i` Internal CO₂ partial pressure in `Pa`, used to compute e_to_c
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_electron_transport! end;

photosystem_electron_transport!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            ppar::Vector{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport!(methods, ps.trait, ps.state, ps.auxil, ppar, p_i; β = β);

photosystem_electron_transport!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport_c3!(methods.C3_AJ_METHOD, methods.COLIMIT_J, pss, psa, ppar, p_i; β = β);

photosystem_electron_transport!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = photosystem_electron_transport_c4!(methods.C4_AJ_METHOD, methods.COLIMIT_J, pst, psa, ppar; β = β);

photosystem_electron_transport_c3!(
            ajm::AjMethodC3Jmax,
            colimj::UnionColimit{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    @. psa.e2c = (p_i - psa.γ_star) / (pss.EFF_1 * p_i + pss.EFF_2 * psa.γ_star);
    @. psa.j_psii = psa.f_psii * psa.ϕ_psii_max * ppar;
    colimited_rate!(β * psa.j_max, psa.j_psii, psa.j, colimj);
    for i in eachindex(psa.e2c)
        isnan(psa.e2c[i]) ? (psa.e2c[i] = 1 / pss.EFF_1) : nothing;
    end;

    return nothing
);

photosystem_electron_transport_c3!(
            ajm::AjMethodC3Vqmax,
            colimj::SerialColimit,
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = (
    _j = psa.cache1;
    @. _j = ppar * (1 - psa.f_psii) * psa.ϕ_psi_max;
    colimited_rate!(β * psa.v_qmax, _j, psa.j_psi, colimj);

    @. psa.η = 1 - psa.η_l / psa.η_c + (3 * p_i + 7 * psa.γ_star) / (4 * p_i + 8 * psa.γ_star) / psa.η_c;
    for i in eachindex(psa.η)
        isnan(psa.η[i]) ? (psa.η[i] = 1 - psa.η_l / psa.η_c + 3 / 4 / psa.η_c) : nothing;
    end;
    @. psa.j_psii = psa.j_psi / psa.η;
    @. psa.j = psa.j_psii;
    @. psa.e2c = (p_i - psa.γ_star) / (4 * p_i + 8 * psa.γ_star);
    for i in eachindex(psa.e2c)
        isnan(psa.e2c[i]) ? (psa.e2c[i] = 1 / 4) : nothing;
    end;

    return nothing
);

photosystem_electron_transport_c4!(
            ajm::AjMethodC4JPSII,
            colimj::UnionColimit{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    @. psa.e2c = 1 / 6;
    @. psa.j_psii = psa.f_psii * psa.ϕ_psii_max * ppar;
    @. psa.j = psa.j_psii;

    return nothing
);

photosystem_electron_transport_c4!(
            ajm::AjMethodC4Vqmax,
            colimj::SerialColimit,
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # partitioning the PPAR to PSI and PSII in mesophyll and bundle sheath
    ppar_psi_m = psa.cache1;
    ppar_psi_bs = psa.cache2;
    ppar_psii_m = psa.cache3;
    ppar_psii_bs = psa.cache4;
    @. ppar_psi_m = ppar * (1 - pst.f_b6f_bs) * (1 - psa.f_psii_m);
    @. ppar_psii_m = ppar * (1 - pst.f_b6f_bs) * psa.f_psii_m;
    @. ppar_psi_bs = ppar * pst.f_b6f_bs * (1 - psa.f_psii_bs);
    @. ppar_psii_bs = ppar * pst.f_b6f_bs * psa.f_psii_bs;

    # the potential J in PSI
    j1_ppar_m = psa.cache5;
    j1_ppar_bs = psa.cache6;
    @. j1_ppar_m = ppar_psi_m * psa.ϕ_psi_max;
    @. j1_ppar_bs = ppar_psi_bs * psa.ϕ_psi_max;
    colimited_rate!(β * psa.v_qmax_m, j1_ppar_m, psa.j_psi_m, colimj);
    colimited_rate!(β * psa.v_qmax_bs, j1_ppar_bs, psa.j_psi_bs, colimj);

    # the potential J in PSII
    @. psa.η_m = 1 - psa.η_l / psa.η_c + 1 / psa.η_c;
    @. psa.η_bs = 1 - psa.η_l / psa.η_c + (3 * psa.p_c_bs + 7 * psa.γ_star) / (4 * psa.p_c_bs + 8 * psa.γ_star) / psa.η_c;
    for i in eachindex(psa.η_bs)
        isnan(psa.η_bs[i]) ? (psa.η_bs[i] = 1 - psa.η_l / psa.η_c + 3 / 4 / psa.η_c) : nothing;
    end;
    @. psa.j_psii_m = psa.j_psi_m / psa.η_m;
    @. psa.j_psii_bs = psa.j_psi_bs / psa.η_bs;
    @. psa.j = psa.j_psii_bs;

    # update the e2c (only in the bundle sheath)
    @. psa.e2c = (psa.p_c_bs - psa.γ_star) / (4 * psa.p_c_bs + 8 * psa.γ_star);
    for i in eachindex(psa.e2c)
        isnan(psa.e2c[i]) ? (psa.e2c[i] = 1 / 4) : nothing;
    end;

    return nothing
);
