using Test

import Photosynthesis as PS


@testset "Fluorescence Methods" verbose = true begin
    @testset "C3 - Vcmax - KN" begin
        # these config settings are the default methods, I am just being explicit here
        methods = PS.PhotosynthesisMethods{Float64}();
        PS.use_c3_jmax_model!(methods);
        lps = PS.LeafPhotosystem(Float64, "C3");
        p_air = 101325.0;   # Pa
        p_co2 = 40.0;       # Pa
        p_i = 20.0;         # Pa
        ppar = 1000.0;      # μmol m⁻² s⁻¹
        t = 298.15;         # K
        g_lc = [0.1];       # mol m⁻² s⁻¹
        PS.photosynthesis!(methods, lps, p_air, [ppar,], t, p_i);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
        PS.photosynthesis!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C3 - Vcmax - QL" begin
        # these config settings are the default methods, I am just being explicit here
        methods = PS.PhotosynthesisMethods{Float64}();
        PS.use_c3_jmax_model!(methods);
        methods.FLUORESCENCE_METHOD_C3 = PS.QLFluorescenceModel{Float64}();
        lps = PS.LeafPhotosystem(Float64, "C3");
        p_air = 101325.0;   # Pa
        p_co2 = 40.0;       # Pa
        p_i = 20.0;         # Pa
        ppar = 1000.0;      # μmol m⁻² s⁻¹
        t = 298.15;         # K
        g_lc = [0.1];       # mol m⁻² s⁻¹
        PS.photosynthesis!(methods, lps, p_air, [ppar,], t, p_i);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
        PS.photosynthesis!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C3 - Vqmax - B6F" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        PS.use_c3_vqmax_model!(methods);
        lps = PS.LeafPhotosystem(Float64, "C3");
        p_air = 101325.0;   # Pa
        p_co2 = 40.0;       # Pa
        p_i = 20.0;         # Pa
        ppar = 1000.0;      # μmol m⁻² s⁻¹
        t = 298.15;         # K
        g_lc = [0.1];       # mol m⁻² s⁻¹
        PS.photosynthesis!(methods, lps, p_air, [ppar,], t, p_i);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
        PS.photosynthesis!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C4 - Vcmax - KN" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        PS.use_c4_vcmax_model!(methods);
        lps = PS.LeafPhotosystem(Float64, "C4");
        p_air = 101325.0;   # Pa
        p_co2 = 40.0;       # Pa
        p_i = 10.0;         # Pa
        ppar = 1000.0;      # μmol m⁻² s⁻¹
        t = 298.15;         # K
        g_lc = [0.1];       # mol m⁻² s⁻¹
        PS.photosynthesis!(methods, lps, p_air, [ppar,], t, p_i);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
        PS.photosynthesis!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C4 - Vcmax - QL" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        PS.use_c4_vcmax_model!(methods);
        methods.FLUORESCENCE_METHOD_C4 = PS.QLFluorescenceModel{Float64}();
        lps = PS.LeafPhotosystem(Float64, "C4");
        p_air = 101325.0;   # Pa
        p_co2 = 40.0;       # Pa
        p_i = 10.0;         # Pa
        ppar = 1000.0;      # μmol m⁻² s⁻¹
        t = 298.15;         # K
        g_lc = [0.1];       # mol m⁻² s⁻¹
        PS.photosynthesis!(methods, lps, p_air, [ppar,], t, p_i);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
        PS.photosynthesis!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;
end;
