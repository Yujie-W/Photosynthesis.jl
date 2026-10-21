"""

$(TYPEDEF)

Struct that contains the fields for C3 photosynthesis

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct LeafPhotosystem{FT}
    "Trait variables"
    trait::Union{C3Trait{FT}, C4Trait{FT}}
    "State variables"
    state::Union{C3State{FT}, C4State{FT}}
    "Auxilary variables"
    auxil::LeafPhotosystemAuxil{FT}
end;

LeafPhotosystem(FT, c3c4::String = "C3"; dim::Int = 1) = (
    @assert c3c4 in ["C3", "C4"] "The model string should be either C3 or C4!";

    return if c3c4 == "C3"
        LeafPhotosystem{FT}(trait = C3Trait{FT}(), state = C3State{FT}(), auxil = LeafPhotosystemAuxil{FT}(dim));
    else
        ps = LeafPhotosystem{FT}(trait = C4Trait{FT}(), state = C4State{FT}(), auxil = LeafPhotosystemAuxil{FT}(dim));
        ps.auxil.f_psii = 0.41;
        ps
    end;
);
