"""

$(TYPEDEF)

Struct that contains the auxiliary variables for leaf photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystemAuxil{FT}
    # photosynthetic rates
    "RubisCO limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_c::Vector{FT}
    "Gross photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_g::Vector{FT}
    "Intermediate photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_i::Vector{FT}
    "Light limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_j::Vector{FT}
    "Net photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_n::Vector{FT}
    "Product limited photosynthetic rate `[μmol m⁻² s⁻¹]`"
    a_p::Vector{FT}

    # electron transport rates
    "Electron to CO₂ coefficient"
    e2c::Vector{FT}
    "Fraction of absorbed light used by PSII ETR"
    f_psii::FT = 0.5
    "Fraction of absorbed light used by PSII ETR in the mesophyll `[-]`"
    f_psii_m::FT = 0.47
    "Fraction of absorbed light used by PSII ETR in the bundle sheath `[-]`"
    f_psii_bs::FT = 0.47
    "Electron transport `[μmol m⁻² s⁻¹]`"
    j::Vector{FT}
    "Maximal electron transport rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    j_max::FT = 0
    "PSI electron transport rate after colimitation"
    j_psi::Vector{FT}
    "PSI electron transport rate after colimitation in the mesophyll"
    j_psi_m::Vector{FT}
    "PSI electron transport rate after colimitation in the bundle sheath"
    j_psi_bs::Vector{FT}
    "Potential Electron Transport Rate in PSII `[μmol m⁻² s⁻¹]`"
    j_psii::Vector{FT}
    "Potential Electron Transport Rate in PSII in the mesophyll `[μmol m⁻² s⁻¹]`"
    j_psii_m::Vector{FT}
    "Potential Electron Transport Rate in PSII in the bundle sheath `[μmol m⁻² s⁻¹]`"
    j_psii_bs::Vector{FT}

    # photosynthesis rate coefficients
    "RubisCO coefficient Kc `[Pa]`"
    k_c::FT = 0
    "Michaelis-Menten's coefficient `[Pa]`"
    k_m::Vector{FT}
    "RubisCO coefficient Ko `[Pa]`"
    k_o::FT = 0
    "PEP coefficient Kpep `[Pa]`"
    k_pep::FT = 0
    "PEP coefficient Kpep fro CLM (different algorithm) `[Pa]`"
    k_pep_clm::FT = 0
    "Maximal turnover rate of Cytochrome b₆f `[e⁻ s⁻¹]`"
    k_q::FT = 0
    "CO₂ compensation point with the absence of Rd `[Pa]`"
    γ_star::FT = 0

    # respiration and carboxylation
    "Respiration rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    r_d::FT = 0
    "Respiration rate at leaf temperature in the mesophyll `[μmol m⁻² s⁻¹]`"
    r_d_m::FT = 0
    "Respiration rate at leaf temperature in the bundle sheath `[μmol m⁻² s⁻¹]`"
    r_d_bs::FT = 0
    "Maximal carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_cmax::FT = 0
    "Maximal PEP carboxylation rate at leaf temperature `[μmol m⁻² s⁻¹]`"
    v_pmax::FT = 0
    "Maximal Cytochrome b₆f activity `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax::FT = 0
    "Maximal Cytochrome b₆f activity in the mesophyll `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax_m::FT = 0
    "Maximal Cytochrome b₆f activity in the bundle sheath `[μmol e⁻ m⁻² s⁻¹]`"
    v_qmax_bs::FT = 0

    # C3 Cytochrome model variables
    "Ratio between J_P700 and J_P680"
    η::Vector{FT}
    "Ratio between J_P700 and J_P680 in the mesophyll"
    η_m::Vector{FT}
    "Ratio between J_P700 and J_P680 in the bundle sheath"
    η_bs::Vector{FT}
    "Coupling efficiency of cyclic electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_c::FT = 0
    "Coupling efficiency of linear electron flow `[mol ATP mol⁻¹ e⁻]`"
    η_l::FT = 0

    # yield variables
    "Heat dissipation yield"
    ϕ_d::Vector{FT}
    "Fluorescence yield"
    ϕ_f::Vector{FT}
    "Non-photochemical yield"
    ϕ_n::Vector{FT}
    "Photochemical yield"
    ϕ_p::Vector{FT}

    # fluorescence yeolds of two photosystems
    "Heat dissipation yield of PSI"
    ϕ_d1::Vector{FT}
    "Heat dissipation yield of PSI in the mesophyll"
    ϕ_d1_m::Vector{FT}
    "Heat dissipation yield of PSI in the bundle sheath"
    ϕ_d1_bs::Vector{FT}
    "Heat dissipation yield of PSII"
    ϕ_d2::Vector{FT}
    "Heat dissipation yield of PSII in the mesophyll"
    ϕ_d2_m::Vector{FT}
    "Heat dissipation yield of PSII in the bundle sheath"
    ϕ_d2_bs::Vector{FT}
    "Fluorescence yield of PSI"
    ϕ_f1::Vector{FT}
    "Fluorescence yield of PSI in the mesophyll"
    ϕ_f1_m::Vector{FT}
    "Fluorescence yield of PSI in the bundle sheath"
    ϕ_f1_bs::Vector{FT}
    "Fluorescence yield of PSII"
    ϕ_f2::Vector{FT}
    "Fluorescence yield of PSII in the mesophyll"
    ϕ_f2_m::Vector{FT}
    "Fluorescence yield of PSII in the bundle sheath"
    ϕ_f2_bs::Vector{FT}
    "Non-photochemical yield of PSI"
    ϕ_n1::Vector{FT}
    "Non-photochemical yield of PSI in the mesophyll"
    ϕ_n1_m::Vector{FT}
    "Non-photochemical yield of PSI in the bundle sheath"
    ϕ_n1_bs::Vector{FT}
    "Non-photochemical yield of PSII"
    ϕ_n2::Vector{FT}
    "Non-photochemical yield of PSII in the mesophyll"
    ϕ_n2_m::Vector{FT}
    "Non-photochemical yield of PSII in the bundle sheath"
    ϕ_n2_bs::Vector{FT}
    "Photochemical yield of PSI"
    ϕ_p1::Vector{FT}
    "Photochemical yield of PSI in the mesophyll"
    ϕ_p1_m::Vector{FT}
    "Photochemical yield of PSI in the bundle sheath"
    ϕ_p1_bs::Vector{FT}
    "Photochemical yield of PSII"
    ϕ_p2::Vector{FT}
    "Photochemical yield of PSII in the mesophyll"
    ϕ_p2_m::Vector{FT}
    "Photochemical yield of PSII in the bundle sheath"
    ϕ_p2_bs::Vector{FT}

    # fluorescence variables
    "Dark adapted yield (`Kp=0`)"
    f_m::FT = 0
    "Light adapted yield (`Kp=0`)"
    f_m′::Vector{FT}
    "Dark-adapted fluorescence yield (`Kp=max`)"
    f_o::FT = 0
    "Light-adapted fluorescence yield in the dark (`Kp=max`)"
    f_o′::Vector{FT}
    "Non-Photochemical quenching "
    npq::Vector{FT}
    "Energy quenching"
    q_e::Vector{FT}
    "Photochemical quenching"
    q_p::Vector{FT}

    # fluorescence rate coefficients
    "Rate constant for thermal dissipation"
    k_d::FT = 0
    "Reversible NPQ rate constant (initially zero)"
    k_n::Vector{FT}
    "Reversible NPQ rate constant (initially zero) in the mesophyll"
    k_n_m::Vector{FT}
    "Reversible NPQ rate constant (initially zero) in the bundle sheath"
    k_n_bs::Vector{FT}
    "Rate constant for photochemistry"
    k_p::Vector{FT}
    "Rate constant for photochemistry in the mesophyll"
    k_p_m::Vector{FT}
    "Rate constant for photochemistry in the bundle sheath"
    k_p_bs::Vector{FT}
    "Maximal PS I photochemical yield"
    ϕ_psi_max::FT = 0
    "max PSII yield (_k_npq_rev = 0, all RC open)"
    ϕ_psii_max::FT = 0

    # internal partial pressures
    "CO₂ partial pressure in the bundle sheath `[Pa]`"
    p_c_bs::Vector{FT}
    "O₂ partial pressure in the bundle sheath `[Pa]`"
    p_o_bs::Vector{FT}

    # cache variables with the dimension of the number of leaves
    cache1::Vector{FT}
    cache2::Vector{FT}
    cache3::Vector{FT}
    cache4::Vector{FT}
    cache5::Vector{FT}
    cache6::Vector{FT}
    cache7::Vector{FT}
    cache8::Vector{FT}
    cache9::Vector{FT}
    cache0::Vector{FT}
end;

LeafPhotosystemAuxil{FT}(dim::Int) where {FT} = (
    return LeafPhotosystemAuxil{FT}(
                a_c       = zeros(FT, dim),
                a_g       = zeros(FT, dim),
                a_i       = zeros(FT, dim),
                a_j       = zeros(FT, dim),
                a_n       = zeros(FT, dim),
                a_p       = zeros(FT, dim),
                e2c       = zeros(FT, dim),
                j         = zeros(FT, dim),
                j_psi     = zeros(FT, dim),
                j_psi_m   = zeros(FT, dim),
                j_psi_bs  = zeros(FT, dim),
                j_psii    = zeros(FT, dim),
                j_psii_m  = zeros(FT, dim),
                j_psii_bs = zeros(FT, dim),
                k_m       = zeros(FT, dim),
                η         = zeros(FT, dim),
                η_m       = zeros(FT, dim),
                η_bs      = zeros(FT, dim),
                ϕ_d       = zeros(FT, dim),
                ϕ_f       = zeros(FT, dim),
                ϕ_n       = zeros(FT, dim),
                ϕ_p       = zeros(FT, dim),
                ϕ_d1      = zeros(FT, dim),
                ϕ_d1_m    = zeros(FT, dim),
                ϕ_d1_bs   = zeros(FT, dim),
                ϕ_d2      = zeros(FT, dim),
                ϕ_d2_m    = zeros(FT, dim),
                ϕ_d2_bs   = zeros(FT, dim),
                ϕ_f1      = zeros(FT, dim),
                ϕ_f1_m    = zeros(FT, dim),
                ϕ_f1_bs   = zeros(FT, dim),
                ϕ_f2      = zeros(FT, dim),
                ϕ_f2_m    = zeros(FT, dim),
                ϕ_f2_bs   = zeros(FT, dim),
                ϕ_n1      = zeros(FT, dim),
                ϕ_n1_m    = zeros(FT, dim),
                ϕ_n1_bs   = zeros(FT, dim),
                ϕ_n2      = zeros(FT, dim),
                ϕ_n2_m    = zeros(FT, dim),
                ϕ_n2_bs   = zeros(FT, dim),
                ϕ_p1      = zeros(FT, dim),
                ϕ_p1_m    = zeros(FT, dim),
                ϕ_p1_bs   = zeros(FT, dim),
                ϕ_p2      = zeros(FT, dim),
                ϕ_p2_m    = zeros(FT, dim),
                ϕ_p2_bs   = zeros(FT, dim),
                f_m′      = zeros(FT, dim),
                f_o′      = zeros(FT, dim),
                npq       = zeros(FT, dim),
                q_e       = zeros(FT, dim),
                q_p       = zeros(FT, dim),
                k_n       = zeros(FT, dim),
                k_n_m     = zeros(FT, dim),
                k_n_bs    = zeros(FT, dim),
                k_p       = zeros(FT, dim),
                k_p_m     = zeros(FT, dim),
                k_p_bs    = zeros(FT, dim),
                p_c_bs    = zeros(FT, dim) .+ 20,
                p_o_bs    = zeros(FT, dim) .+ F_O₂() .* P_ATM(),
                cache1    = zeros(FT, dim),
                cache2    = zeros(FT, dim),
                cache3    = zeros(FT, dim),
                cache4    = zeros(FT, dim),
                cache5    = zeros(FT, dim),
                cache6    = zeros(FT, dim),
                cache7    = zeros(FT, dim),
                cache8    = zeros(FT, dim),
                cache9    = zeros(FT, dim),
                cache0    = zeros(FT, dim),
    )
);
