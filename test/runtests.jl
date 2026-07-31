using Test


@testset verbose = true "Photosynthesis" begin
    @testset "Tutorials" verbose = true begin
        include("tutorials/photosynthesis-methods.jl");
        include("tutorials/fluorescence-methods.jl");
    end;
end;
