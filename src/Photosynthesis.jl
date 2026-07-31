module Photosynthesis

using PkgUtility.UniversalConstants: T₀, T₂₅


# types and structs
include("types/photosystem-rate-constant.jl");
include("types/temperature-dependency.jl");
include("types/temperature-dependency-settings.jl");


end; # module
