"""

    photosynthesis!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                p_air::FT,
                p_co2::FT,
                ppar::Vector{FT},
                t::FT,
                g_lc::Vector{FT}
        ) where {FT}
    photosynthesis!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                p_air::FT,
                ppar::Vector{FT},
                t::FT,
                p_i::Union{FT,Vector{FT}}
        ) where {FT}

Update the photosynthesis rate, given
- `methods` `PhotosynthesisMethods` structure that contains rate constants, temperature dependency methods, colimitation methods, and model features
- `ps` `LeafPhotosystem` struct
- `p_air` Surface atmosphere pressure in `Pa`
- `p_co2` Surface CO₂ partial pressure in `Pa`
- `ppar` Vector of photosynthetically active radiation
- `t` Leaf temperature in `K`
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]` (if appears in the function signature, it is meant for conductance mode)
- `p_i` Vector of intercellular CO₂ partial pressures

"""
function photosynthesis! end;

photosynthesis!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_air::FT,
            p_co2::FT,
            ppar::Vector{FT},
            t::FT,
            g_lc::Vector{FT}) where {FT} = (
    photosystem_temperature_dependence!(methods, ps, p_air, t);
    photosystem_electron_transport!(methods, ps, ppar, p_co2 / 2);
    rubisco_limited_rate!(methods, ps, p_air, p_co2, g_lc);
    light_limited_rate!(methods, ps, p_air, p_co2, g_lc);
    product_limited_rate!(methods, ps, p_air, p_co2, g_lc);
    colimit_photosynthesis!(methods, ps);
    p_i = @. p_co2 - ps.auxil.a_n / g_lc * p_air .* FT(1e-6);
    photosystem_electron_transport!(methods, ps, ppar, p_i);
    photosystem_coefficients!(methods, ps, ppar);

    return nothing
);

photosynthesis!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            p_air::FT,
            ppar::Vector{FT},
            t::FT,
            p_i::Union{FT,Vector{FT}};
            f_o2::Number = F_O₂()) where {FT} = (
    photosystem_temperature_dependence!(methods, ps, p_air, t);
    photosystem_electron_transport!(methods, ps, ppar, p_i);
    rubisco_limited_rate!(methods, ps, p_i);
    light_limited_rate!(ps);
    product_limited_rate!(methods, ps, p_i);
    colimit_photosynthesis!(methods, ps);
    photosystem_coefficients!(methods, ps, ppar);

    # update p_c_bs and p_o_bs for steady state solution
    if ps.trait isa C4Trait{FT} && methods.C4_AC_METHOD isa AcMethodC4Vqmax
        p_c = ps.auxil.cache1;
        p_o = ps.auxil.cache2;
        @. p_c = p_i + (ps.auxil.a_p - ps.auxil.a_n) / 1000000 / ps.trait.g_bs_c * p_air;
        @. p_o = p_air * f_o2 + (ps.auxil.a_n - ps.auxil.a_p / 2) / 1000000 / ps.trait.g_bs_o * p_air;
        @. ps.auxil.p_c_bs += (p_c - ps.auxil.p_c_bs) / 10;
        @. ps.auxil.p_o_bs += (p_o - ps.auxil.p_o_bs) / 10;
    end;

    return nothing
);
