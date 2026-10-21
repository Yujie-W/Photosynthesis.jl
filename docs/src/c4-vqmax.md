# C4-Vqmax model

## CO₂ Response
```@example
import Photosynthesis as PS
using Plots

methods = PS.PhotosynthesisMethods{Float64}();
methods.C4_AC_METHOD = PS.AcMethodC4Vqmax();
methods.C4_AJ_METHOD = PS.AjMethodC4Vqmax();
methods.C4_AP_METHOD = PS.ApMethodC4Vqmax();
methods.COLIMIT_J = PS.SerialColimit();
methods.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
p_air = 101325.0;                   # Pa
p_i = collect(0.5:0.5:40);            # Pa
ppar = ones(length(p_i)) * 1000;    # μmol m⁻² s⁻¹
t = 298.15;                         # K
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
# run N times to reach steady state solution
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(p_i, lps.auxil.a_c, label="Ac")
plot!(p_i, lps.auxil.a_j, label="Aj")
plot!(p_i, lps.auxil.a_p, label="Ap")
```

## Light Responses
```@example
import Photosynthesis as PS
using Plots

methods = PS.PhotosynthesisMethods{Float64}();
methods.C4_AC_METHOD = PS.AcMethodC4Vqmax();
methods.C4_AJ_METHOD = PS.AjMethodC4Vqmax();
methods.C4_AP_METHOD = PS.ApMethodC4Vqmax();
methods.COLIMIT_J = PS.SerialColimit();
methods.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
p_air = 101325.0;                   # Pa
ppar = collect(0.0:50:2000);       # μmol m⁻² s⁻¹
p_i = ones(length(ppar)) .* 20;     # Pa
t = 298.15;                         # K
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
# run N times to reach steady state solution
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(ppar, lps.auxil.a_c, label="Ac")
plot!(ppar, lps.auxil.a_j, label="Aj")
plot!(ppar, lps.auxil.a_p, label="Ap")
```

## Temperature Impacts
```@example
import Photosynthesis as PS
using Plots

methods = PS.PhotosynthesisMethods{Float64}();
methods.C4_AC_METHOD = PS.AcMethodC4Vqmax();
methods.C4_AJ_METHOD = PS.AjMethodC4Vqmax();
methods.C4_AP_METHOD = PS.ApMethodC4Vqmax();
methods.COLIMIT_J = PS.SerialColimit();
methods.FLUORESCENCE_METHOD_C4 = PS.KNFluorescenceModel{Float64}();
p_air = 101325.0;                   # Pa
ppar = collect(0.0:50:2000);       # μmol m⁻² s⁻¹
p_i = ones(length(ppar)) .* 20;     # Pa
t = 313.15;                         # K
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
# run N times to reach steady state solution
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(ppar, lps.auxil.a_c, label="Ac")
plot!(ppar, lps.auxil.a_j, label="Aj")
plot!(ppar, lps.auxil.a_p, label="Ap")
```
