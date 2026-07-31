module Photosynthesis

using PkgUtility.UniversalConstants: T₀, T₂₅


# types and structs
include("types/rate-constant.jl");
include("types/temperature-dependency.jl");
include("types/temperature-dependency-settings.jl");
include("types/photosystem-auxil.jl");
include("types/photosystem-c3.jl");
include("types/photosystem-c4.jl");
include("types/photosystem.jl");

# core model
include("model/colimit.jl");
include("model/etr.jl");
include("model/fluorescence.jl");
include("model/light-limited.jl");
include("model/photosystem.jl");
include("model/product-limited.jl");
include("model/rubisco-limited.jl");
include("model/td.jl");


end; # module
