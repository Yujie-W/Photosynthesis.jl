module Photosynthesis

using PkgUtility.UniversalConstants: F_O₂, GAS_R, T₀, T₂₅


# types and structs
include("types/colimitation-method.jl");
include("types/colimitation-method-settings.jl");
include("types/fluorescence-model.jl");
include("types/fluorescence-model-settings.jl");
include("types/photosynthesis-model.jl");
include("types/photosystem-auxil.jl");
include("types/photosystem-c3.jl");
include("types/photosystem-c4.jl");
include("types/photosystem.jl");
include("types/rate-constant.jl");
include("types/temperature-dependency.jl");
include("types/temperature-dependency-settings.jl");

include("types/methods.jl");


# core model
include("model/1-td.jl");
include("model/2-etr.jl");
include("model/3-rubisco-limited.jl");
include("model/colimit.jl");
include("model/fluorescence.jl");
include("model/light-limited.jl");
include("model/photosystem.jl");
include("model/product-limited.jl");


end; # module
