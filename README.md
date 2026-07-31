# Photosynthesis.jl

<!-- Links and shortcuts -->
[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://Yujie-W.github.io/Photosynthesis.jl/dev/

[st-img]: https://github.com/Yujie-W/Photosynthesis.jl/workflows/JuliaStable/badge.svg?branch=main
[st-url]: https://github.com/Yujie-W/Photosynthesis.jl/actions?query=branch%3A"main"++workflow%3A"JuliaStable"

[cov-img]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl/branch/main/graph/badge.svg
[cov-url]: https://codecov.io/gh/Yujie-W/Photosynthesis.jl


## About
Photosynthesis models for modeling C3 and C4 photosynthesis.

| Documentation           | CI Status             | Code Coverage           |
|:-----------------------:|:---------------------:|:------------------------|
| [![][dev-img]][dev-url] | [![][st-img]][st-url] | [![][cov-img]][cov-url] |

`Photosynthesis.jl` is meant to be a module of [Emerald (formerly CliMA Land)](https://github.com/Yujie-W/EmeraldDeveloping.git). However, it should be equally convenient to use it for other purposes.

### To cite Photosynthesis.jl
- Wang et al. (2024). Beyond the visible: Accounting for ultraviolet and far-red radiation in vegetation productivity and surface energy budgets. Global Change Biology 30(5): e17346. [DOI: 10.1111/gcb.17346](10.1111/gcb.17346)

- Wang et al. (2023). Modeling global vegetation gross primary productivity, transpiration and hyperspectral canopy radiative transfer simultaneously using a next generation land surface model—CliMA Land. Journal of Advances in Modeling Earth Systems 15(3): e2021MS002964. [DOI: 10.1029/2021MS002964](10.1029/2021MS002964)

- Wang et al. (2021). Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1). Geoscientific Model Development 14(11): 6741–6763. [DOI: 10.5194/gmd-14-6741-2021](10.5194/gmd-14-6741-2021)


## Installation
```julia
julia> using Pkg;
julia> Pkg.add("Photosynthesis");
```
