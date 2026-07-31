module Photosynthesis

using PkgUtility.MathTools: lower_quadratic
using PkgUtility.UniversalConstants: F_O₂, GAS_R, T₀, T₂₅


# export public functions
export LeafPhotosystem, PhotosynthesisMethods
export photosynthesis!, photosynthesis_glc!


# types and structs
include("types/colimitation-method.jl");
include("types/fluorescence-model.jl");
include("types/photosynthesis-model.jl");
include("types/rate-constant.jl");
include("types/temperature-dependency.jl");

include("types/photosystem-auxil.jl");
include("types/photosystem-c3.jl");
include("types/photosystem-c4.jl");
include("types/photosystem.jl");

include("constructors/colimitation-method-settings.jl");
include("constructors/fluorescence-model-settings.jl");
include("constructors/temperature-dependency-settings.jl");
include("types/methods.jl");


# core model
include("model/1-td.jl");
include("model/2-etr.jl");
include("model/3-rubisco-limited.jl");
include("model/4-light-limited.jl");
include("model/5-product-limited.jl");
include("model/6-colimit.jl");
include("model/7-fluorescence.jl");
include("model/pipeline.jl");


end; # module
