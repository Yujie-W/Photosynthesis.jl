"""
Abstract type for photosynthesis rate calculation methods for Ac:
- AcMethodC3Vcmax
- AcMethodC4Vcmax
- AcMethodC4Vqmax
"""
abstract type AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):
```math
A_\\text{c} = V_{\\text{cmax}} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + K_{\\text{m}})}
```
"""
struct AcMethodC3Vcmax <: AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac):
```math
A_\\text{c} = V_{\\text{cmax}}
```
"""
struct AcMethodC4Vcmax <: AbstractAcMethod end;


"""
Method to compute rubisco-limited photosynthesis rate (Ac) using the Johnson et al. (2021) model for C4 photosynthesis
"""
struct AcMethodC4Vqmax <: AbstractAcMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Aj:
- AjMethodC3Jmax
- AjMethodC3Vqmax
- AjMethodC4JPSII
- AjMethodC4Vqmax
"""
abstract type AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj):
```math
A_\\text{j} = \\dfrac{J}{4} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + 2 \\Gamma)}
```
"""
struct AjMethodC3Jmax <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:
```math
J = \\dfrac{J_{\\text{PSI}}}{\\eta}
A_{\\text{j}} = \\dfrac{J}{4} \\dfrac{(P_{\\text{i}} - \\Gamma)}{(P_{\\text{i}} + 2 \\Gamma)}
```
"""
struct AjMethodC3Vqmax <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using cytochrome b6f:
```math
Aj = \\dfrac{J}{6}
```
"""
struct AjMethodC4JPSII <: AbstractAjMethod end;


"""
Method to compute electron transport-limited photosynthesis rate (Aj) using the Johnson et al. (2021) model for C4 photosynthesis
"""
struct AjMethodC4Vqmax <: AbstractAjMethod end;


"""
Abstract type for photosynthesis rate calculation methods for Ap:
- ApMethodC3Inf
- ApMethodC3Vcmax
- ApMethodC4Vcmax
- ApMethodC4Vpmax
- ApMethodC4Vqmax
"""
abstract type AbstractApMethod end;


"""
Method for models without TPU limitation:
```math
A_{\\text{p}} = \\infty
```
"""
struct ApMethodC3Inf <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:
```math
A_{\\text{p}} = \\dfrac{V_{\\text{cmax}}}{2}
```
"""
struct ApMethodC3Vcmax <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vcmax:
```math
A_{\\text{p}} = k \\cdot V_{\\text{cmax}} \\cdot P_{\\text{i}}
```
"""
struct ApMethodC4Vcmax <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) from Vpmax:
```math
A_{\\text{p}} = V_{\\text{pmax}} \\dfrac{P_{\\text{i}}}{(P_{\\text{i}} + K_{\\text{pep}})}
```
"""
struct ApMethodC4Vpmax <: AbstractApMethod end;


"""
Method to compute product-limited photosynthesis rate (Ap) using the Johnson et al. (2021) model for C4 photosynthesis
"""
struct ApMethodC4Vqmax <: AbstractApMethod end;
