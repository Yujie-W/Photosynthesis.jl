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
    "Respiration rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    r_d25::FT = 0.75
    "Maximal carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_cmax25::FT = 50
    "Maximal PEP carboxylation rate at 298.15 K `[μmol m⁻² s⁻¹]`"
    v_pmax25::FT = 50
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
