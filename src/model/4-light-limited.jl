"""

    light_limited_rate!(
                ps::LeafPhotosystem{FT},
                air::AirLayer{FT},
                g_lc::Vector{FT};
                β::FT = FT(1)) where {FT}

Update the electron transport limited photosynthetic rate, given
- `ps` `LeafPhotosystem` struct
- `air` `AirLayer` struct for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function light_limited_rate! end;

light_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate!(methods, ps.state, ps.auxil, p_air, p_co2, g_lc; β = β);

light_limited_rate!(ps::LeafPhotosystem{FT}) where {FT} = light_limited_rate!(ps.auxil);

light_limited_rate!(psa::LeafPhotosystemAuxil{FT}) where {FT} = (@. psa.a_j = psa.j * psa.e2c; return nothing);

light_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate_c3!(pss, psa, methods.C3_AJ_METHOD, p_air, p_co2, g_lc; β = β);

light_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = light_limited_rate!(psa);

light_limited_rate_c3!(
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::AjMethodC3JmaxPi,
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j[1] == 0 && psa.j[end] == 0
        @. psa.a_j = 0;

        return nothing
    end;

    # unpack the cache variables
    b  = psa.cache1;
    f  = psa.cache2;
    qa = psa.cache3;
    qb = psa.cache4;
    qc = psa.cache5;
    an = psa.cache6;

    a = psa.j;
    c = pss.EFF_1;
    d = pss.EFF_2 * psa.γ_star;
    p = p_co2;
    r = β * psa.r_d;
    @. b = psa.j * psa.γ_star;
    @. f = p_air / g_lc * FT(1e-6);

    @. qa = c * f;
    @. qb = c * f * r - c * p - d - a * f;
    @. qc = a * p - b - r * (c * p + d);
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_j[i] = r;
        else
            psa.a_j[i] = an[i] + r;
        end;
    end;

    return nothing
);

light_limited_rate_c3!(
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::AjMethodC3VqmaxPi,
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if psa.j_psi[1] == 0 && psa.j_psi[end] == 0
        @. psa.a_j = 0;

        return nothing
    end;

    eff_a = 1 - psa.η_l / psa.η_c;
    eff_b = 1 / psa.η_c;
    eff_1 = eff_a * pss.EFF_1 + 3 * eff_b;
    eff_2 = eff_a * pss.EFF_2 + 7 * eff_b;

    # unpack the cache variables
    b  = psa.cache1;
    f  = psa.cache2;
    qa = psa.cache3;
    qb = psa.cache4;
    qc = psa.cache5;
    an = psa.cache6;

    a  = psa.j_psi;
    c  = eff_1;
    d  = eff_2 * psa.γ_star;
    p  = p_co2;
    r  = β * psa.r_d;
    @. b = psa.j_psi * psa.γ_star;
    @. f = p_air / g_lc * FT(1e-6);

    @. qa = c * f;
    @. qb = c * f * r - c * p - d - a * f;
    @. qc = a * p - b - r * (c * p + d);
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_j[i] = r;
        else
            psa.a_j[i] = an[i] + r;
        end;
    end;

    return nothing
);
