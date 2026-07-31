using Test

import Photosynthesis as PS


@testset "Photosynthesis Methods" verbose = true begin
    @testset "C3 - Jmax" begin
        # these config settings are the default methods, I am just being explicit here
        methods = PS.PhotosynthesisMethods{Float64}();
        methods.C3_AC_METHOD = PS.AcMethodC3VcmaxPi();
        methods.C3_AJ_METHOD = PS.AjMethodC3JmaxPi();
        methods.C3_AP_METHOD = PS.ApMethodC3Vcmax();
        methods.COLIMIT_J = PS.ColimitJCLM(Float64);
        methods.FLUORESCENCE_METHOD_C3 = PS.KNFluorescenceModel{Float64}();
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
        PS.photosynthesis_glc!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C3 - Vqmax" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        methods.C3_AC_METHOD = PS.AcMethodC3VcmaxPi();
        methods.C3_AJ_METHOD = PS.AjMethodC3VqmaxPi();
        methods.C3_AP_METHOD = PS.ApMethodC3Vcmax();
        methods.COLIMIT_J = PS.SerialColimit();
        methods.FLUORESCENCE_METHOD_C3 = PS.CytochromeFluorescenceModel();
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
        PS.photosynthesis_glc!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C4 - Vcmax" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        methods.C4_AC_METHOD = PS.AcMethodC4Vcmax();
        methods.C4_AJ_METHOD = PS.AjMethodC4JPSII();
        methods.C4_AP_METHOD = PS.ApMethodC4VcmaxPi();
        methods.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
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
        PS.photosynthesis_glc!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;

    @testset "C4 - Vpmax" begin
        methods = PS.PhotosynthesisMethods{Float64}();
        methods.C4_AC_METHOD = PS.AcMethodC4Vcmax();
        methods.C4_AJ_METHOD = PS.AjMethodC4JPSII();
        methods.C4_AP_METHOD = PS.ApMethodC4VpmaxPi();
        methods.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
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
        PS.photosynthesis_glc!(methods, lps, p_air, p_co2, [ppar,], t, g_lc);
        @test lps.auxil.a_n[1] > 0;
        @test lps.auxil.ϕ_f[1] > 0;
    end;
end;
