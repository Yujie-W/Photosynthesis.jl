# C4-Vqmax model

```@example c4vqmax
import Photosynthesis as PS
using Plots

methods = PS.PhotosynthesisMethods{Float64}();
PS.use_c4_vqmax_model!(methods);
```


## CO₂ Response
```@example c4vqmax
p_air = 101325.0;                   # Pa
t = 298.15;                         # K
p_i = collect(0.5:0.5:40);          # Pa
ppar = ones(length(p_i)) * 1000;    # μmol m⁻² s⁻¹
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(p_i, lps.auxil.a_c, label="Ac")
plot!(p_i, lps.auxil.a_j, label="Aj")
plot!(p_i, lps.auxil.a_p, label="Ap")
plot!(xlabel = "CO₂ (Pa)", ylabel = "A (μmol m⁻² s⁻¹)")
```

## Light Responses
```@example c4vqmax
p_air = 101325.0;                   # Pa
t = 298.15;                         # K
ppar = collect(5.0:50:2000);        # μmol m⁻² s⁻¹
p_i = ones(length(ppar)) .* 20;     # Pa
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(ppar, lps.auxil.a_c, label="Ac")
plot!(ppar, lps.auxil.a_j, label="Aj")
plot!(ppar, lps.auxil.a_p, label="Ap")
plot!(xlabel = "PPAR (μmol m⁻² s⁻¹)", ylabel = "A (μmol m⁻² s⁻¹)")
```

## Temperature Impacts
```@example c4vqmax
p_air = 101325.0;                   # Pa
t = 298.15;                         # K
ppar = collect(5.0:50:2000);        # μmol m⁻² s⁻¹
p_i = ones(length(ppar)) .* 20;     # Pa
t = 313.15;                         # K
lps = PS.LeafPhotosystem(Float64, "C4"; dim = length(ppar));
for _ in 1:15
    PS.photosynthesis!(methods, lps, p_air, ppar, t, p_i);
end;

p = plot()
plot!(ppar, lps.auxil.a_c, label="Ac")
plot!(ppar, lps.auxil.a_j, label="Aj")
plot!(ppar, lps.auxil.a_p, label="Ap")
plot!(xlabel = "PPAR (μmol m⁻² s⁻¹)", ylabel = "A (μmol m⁻² s⁻¹)")
```
