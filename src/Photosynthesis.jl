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


end; # module
