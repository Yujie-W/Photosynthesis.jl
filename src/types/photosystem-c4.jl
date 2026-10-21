"""

$(TYPEDEF)

Struct that contains the trait variables for C4 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4Trait{FT}
    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 8862 * 1e-6

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 2
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 100
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 100

    # fraction of b₆f and respiration in the mesophyll and bundle sheath
    "Fraction of b₆f in the bundle sheath (analogous to the antenna and pigments distribution) `[-]`"
    f_b6f_bs::FT = 0.6
    "Fraction of respiration in the bundle sheath `[-]`"
    f_rd_bs::FT = 0.4

    # Bundle sheath conductance to CO₂ and O₂
    "Bundle sheath conductance to CO₂ `[mol m⁻² s⁻¹]`"
    g_bs_c::FT = 0.003
    "Bundle sheath conductance to O₂ `[mol m⁻² s⁻¹]`"
    g_bs_o::FT = 0.003 * 0.047
end;


"""

$(TYPEDEF)

Struct that contains the state variables for C4 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C4State{FT}
    # Prognostic variables
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;
