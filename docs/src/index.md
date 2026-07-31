# Photosynthesis.jl

## Installation
```julia
using Pkg;
Pkg.add("Photosynthesis");
```

## About
Photosynthesis models for C3 and C4 photosynthesis. Photosynthesis.jl supports foour photosynthesis models and three fluorescence models. The photosynthesis models are
- C3-Jmax model based on classic C3 model, which is known as FvCB model (Farquhar et al. 1980)
- C3-Vqmaxe model based on a new C3 model developed by Johnson and Berry (2021)
- C4-Vcmax model based on classic C4 model, which is known as Collaz model (Collaz et al. 1992)
- C4-Vpmax model based on PEP activity as in Boyd et al. (2015)

We, however, made some modifications by adding a product limited photosynthetic rate to the C3 models, and a Rubisco limited photosynthetic rate to the C4 model.

Besides the tranditional photosynthesis model, we also included functions to compute fluorescence related parameters, such as fluorescence quantum yield and non-photochemical quenching. he two implemented fluorescence models are
- Van der Tol et al. (2013) fluorescence model to use with C3-Jmax and C4 models
- Han et al. (2022) fluorescence model to used with C3-Jmax and C4 models
- Johnson and Berry (2021) fluorescence model to use with C3Cytochrome model

Starting from v0.5, photosynthesis and fluorescence model selection is done by setting up the fields of a `PhotosynthesisMethods` struct (see the tutorial in the `tests` folder).

We aim to make Photosynthesis.jl a standalone package rather than just part of the Emerald model. Thus, in the documentations below, we will present examples of how to use Photosynthesis.jl at the leaf level. Same logic applies to canopy scale simulations.


## Notes
In our model, we used PPAR, rather than PAR or APAR. PPAR is the photosynthetically active radiation (PAR) that goes to photosystems; PPAR is different from APAR (PAR absorbed by leaf, some APAR does not go to photosystems).


## Photosynthesis Model Procedure
For all three photosynthesis+fluorescence combo models, photosynthetic rates are computed using the following procedure:
- Update temperature dependent variables using `photosystem_temperature_dependence!`
- Calculate electron transport rate using `photosystem_electron_transport!`
- Calculate RubisCO limited rate using `rubisco_limited_rate!`
- Calculate light limited rate using `light_limited_rate!`
- Calculate product limited rate using `product_limited_rate!`
- Calculate gross and net rates using `colimit_photosynthesis!`
- Update fluorescence related variables using `photosystem_coefficients!`

Yet, for convenience, all the listed steps are combined in one function `photosynthesis!`.

### 0. Key Structures and Methods
```@docs
Photosynthesis.C3Trait
Photosynthesis.C4Trait
Photosynthesis.C3State
Photosynthesis.C4State
Photosynthesis.LeafPhotosystemAuxil
Photosynthesis.LeafPhotosystem
Photosynthesis.PhotosynthesisMethods
```

### 1. Temperature Dependencies
```@docs
Photosynthesis.AbstractTemperatureDependency
Photosynthesis.Arrhenius
Photosynthesis.ArrheniusPeak
Photosynthesis.ArrheniusPeak2
Photosynthesis.Q10
Photosynthesis.Q10Peak
Photosynthesis.Q10PeakHT
Photosynthesis.Q10PeakLTHT
Photosynthesis.temperature_correction
Photosynthesis.temperature_corrected_value
Photosynthesis.photosystem_temperature_dependence!
```

### 2. Electron Transport
```@docs
Photosynthesis.PSRateConstants
Photosynthesis.PSIRateConstants
Photosynthesis.PSIIRateConstants
Photosynthesis.photosystem_electron_transport!
```

### 3. RubisCO Limited Rate
```@docs
Photosynthesis.AbstractAcMethod
Photosynthesis.AcMethodC3VcmaxPi
Photosynthesis.AcMethodC4Vcmax
Photosynthesis.rubisco_limited_rate!
```

### 4. Light Limited Rate
```@docs
Photosynthesis.AbstractAjMethod
Photosynthesis.AjMethodC3JmaxPi
Photosynthesis.AjMethodC3VqmaxPi
Photosynthesis.AjMethodC4JPSII
Photosynthesis.light_limited_rate!
```

### 5. Product Limited Rate
```@docs
Photosynthesis.AbstractApMethod
Photosynthesis.ApMethodC3Vcmax
Photosynthesis.ApMethodC3Inf
Photosynthesis.ApMethodC4VcmaxPi
Photosynthesis.ApMethodC4VpmaxPi
Photosynthesis.product_limited_rate!
```

### 6. Colimitation
```@docs
Photosynthesis.AbstractColimit
Photosynthesis.MinimumColimit
Photosynthesis.QuadraticColimit
Photosynthesis.SerialColimit
Photosynthesis.SquareColimit
Photosynthesis.colimited_rate
Photosynthesis.colimited_rate!
Photosynthesis.colimit_photosynthesis!
```

### 7. Fluorescence Related Parameters
```@docs
Photosynthesis.AbstractFluorescenceMethod
Photosynthesis.KNFluorescenceModel
Photosynthesis.QLFluorescenceModel
Photosynthesis.QLFluorescenceModelHan
Photosynthesis.CytochromeFluorescenceModel
Photosynthesis.photosystem_coefficients!
```

### Full Pipeline
```@docs
Photosynthesis.photosynthesis!
```
