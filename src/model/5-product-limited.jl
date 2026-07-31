"""

    product_limited_rate!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                p_air::FT,
                p_co2::FT,
                g_lc::Vector{FT};
                β::FT = FT(1)
        ) where {FT}
    product_limited_rate!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                p_i::Union{FT,Vector{FT}};
                β::FT = FT(1)
        ) where {FT}

Update the product limited photosynthetic rate, given
- `methods` `PhotosynthesisMethods` structure that contains rate constants, temperature dependency methods, colimitation methods, and model features
- `ps` `LeafPhotosystem` struct
- `p_air` Surface atmosphere pressure in `Pa`
- `p_co2` Surface CO₂ partial pressure in `Pa`
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]` (if appears in the function signature, it is meant for conductance mode)
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd
- `p_i` Internal CO₂ partial pressure in `Pa`

"""
function product_limited_rate! end;

product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate!(methods, ps.trait, ps.auxil, p_air, p_co2, g_lc; β = β);

product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c3!(psa, methods.C3_AP_METHOD; β = β);

product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = product_limited_rate_c4!(pst, psa, methods.C4_AP_METHOD, p_air, p_co2, g_lc; β = β);

product_limited_rate_c3!(
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC3Inf;
            β::FT = FT(1)) where {FT} = (@. psa.a_p = FT(Inf); return nothing);

product_limited_rate_c3!(
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC3Vcmax;
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.v_cmax / 2; return nothing);

product_limited_rate_c4!(
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VcmaxPi,
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack variables from the cache
    g   = psa.cache1;
    p_i = psa.cache2;

    a = p_air;
    k = β * psa.k_pep_clm * pst.v_cmax25;
    p = p_co2;
    r = β * psa.r_d;
    @. g = FT(1e6) * g_lc;

    @. p_i = (g * p + a * r) / (a * k + g);
    @. psa.a_p = k * p_i;

    return nothing
);

product_limited_rate_c4!(
            ::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VpmaxPi,
            p_air::FT,
            p_co2::FT,
            g_lc::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    # unpack the variables from the cache
    f  = psa.cache1;
    qb = psa.cache2;
    an = psa.cache3;

    a = β * psa.v_pmax;
    d = psa.k_pep;
    p = p_co2;
    r = β * psa.r_d;
    @. f = p_air / g_lc * FT(1e-6);

    qa = f;
    qc = a * p - r * (p + d);
    @. qb = f * r - p - d - a * f;
    @. an = lower_quadratic(qa, qb, qc);

    for i in eachindex(g_lc)
        if g_lc[i] == 0
            psa.a_p[i] = r;
        else
            psa.a_p[i] = an[i] + r;
        end;
    end;

    return nothing
);

# Pressure mode for A-Ci fiting
product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = product_limited_rate!(methods, ps.trait, ps.auxil, p_i; β = β);

product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            ::C3Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = product_limited_rate_c3!(psa, methods.C3_AP_METHOD; β = β);

product_limited_rate!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = product_limited_rate_c4!(pst, psa, methods.C4_AP_METHOD, p_i; β = β);

product_limited_rate_c4!(
            pst::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VcmaxPi,
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.k_pep_clm * pst.v_cmax25 * p_i; return nothing);

product_limited_rate_c4!(
            ::C4Trait{FT},
            psa::LeafPhotosystemAuxil{FT},
            ::ApMethodC4VpmaxPi,
            p_i::Union{FT,Vector{FT}};
            β::FT = FT(1)) where {FT} = (@. psa.a_p = β * psa.v_pmax * p_i / (p_i + psa.k_pep); return nothing);
