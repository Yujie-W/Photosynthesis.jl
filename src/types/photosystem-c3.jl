"""

$(TYPEDEF)

Struct that contains the trait variables for C3 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3Trait{FT}
    # Related to OCS uptake
    "Multiplier to derive internal conductance for OCS `[mol μmol⁻¹]`"
    K_OCS::FT = 1400 * 1e-6

    # Prognostic variables
    "Total concentration of Cytochrome b₆f `[μmol m⁻²]`"
    b₆f::FT = 350 / 300
    "Maximal electron transport rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    j_max25::FT = 83.5
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
end;


"""

$(TYPEDEF)

Struct that contains the state variables for C3 photosynthesis (VJP model)

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct C3State{FT}
    # General model information
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_1::FT = 4
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    EFF_2::FT = 8

    # Prognostic variables (for VJP model)
    "Sustained NPQ rate constant (for seasonal changes, default is zero)"
    k_npq_sus::FT = 0
end;
