"""

    photosynthesis!(config::SPACConfig{FT}, leaf::Leaf{FT}, air::AirLayer{FT}; use_glc::Bool = false) where {FT}
    photosynthesis!(config::SPACConfig{FT}, ps::LeafPhotosystem{FT}, air::AirLayer{FT}, p_i::Vector{FT}, ppar::Vector{FT}, t::FT) where {FT}

Update the photosynthesis rate, given
- `config` `SPACConfig` struct
- `leaf` `Leaf` struct
- `air` `AirLayer` struct
- `use_glc` Whether to use leaf conductance to CO₂ to compute photosynthesis rate (otherwise use intercellular CO₂ partial pressure)
- `ps` `LeafPhotosystem` struct
- `p_i` Vector of intercellular CO₂ partial pressures
- `ppar` Vector of photosynthetically active radiation
- `t` Leaf temperature in Kelvin

"""
function photosynthesis! end;

photosynthesis!(methods::PhotosynthesisMethods{FT}, ps::LeafPhotosystem{FT}, p_air::FT, ppar::Vector{FT}, t::FT, p_i::Union{FT,Vector{FT}}) where {FT} = (
    photosystem_temperature_dependence!(methods, ps, p_air, t);
    photosystem_electron_transport!(methods, ps, ppar, p_i);
    rubisco_limited_rate!(methods, ps, p_i);
    light_limited_rate!(ps);
    product_limited_rate!(methods, ps, p_i);
    colimit_photosynthesis!(methods, ps);
    photosystem_coefficients!(methods, ps, ppar);

    return nothing
);

photosynthesis_glc!(methods::PhotosynthesisMethods{FT}, ps::LeafPhotosystem{FT}, p_air::FT, p_co2::FT, ppar::Vector{FT}, t::FT, g_lc::Vector{FT}) where {FT} = (
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
