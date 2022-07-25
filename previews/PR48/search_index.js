var documenterSearchIndex = {"docs":
[{"location":"#Photosynthesis.jl","page":"Home","title":"Photosynthesis.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Photosynthesis models for C3 and C4 photosynthesis.","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Pkg;\nPkg.add(\"Photosynthesis\");","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"leaf   = Leaf{Float32}(\"C3\");\nair    = AirLayer{Float32}();\np_mode = PCO₂Mode();\ng_mode = PCO₂Mode();\nleaf_photosynthesis!(leaf, air, p_mode);\nleaf_photosynthesis!(leaf, air, g_mode);","category":"page"},{"location":"API/#Photosynthesis","page":"API","title":"Photosynthesis","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CurrentModule = Photosynthesis","category":"page"},{"location":"API/#Photosynthesis-model","page":"API","title":"Photosynthesis model","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"leaf_photosynthesis!\nleaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT) where {FT<:AbstractFloat}\nleaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}\nleaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}","category":"page"},{"location":"API/#Photosynthesis.leaf_photosynthesis!","page":"API","title":"Photosynthesis.leaf_photosynthesis!","text":"Per refactored Photosynthesis module, the only things one need to know is the public function leaf_photosynthesis! and some construtors from ClimaCache. See the examples in the methods below for     details about how to use the function. The steps for computing photosynthetic rates are\n\nUpdate temperature dependent variables using photosystem_temperature_dependence!\nCalculate electron transport rate using photosystem_electron_transport!\nCalculate RubisCO limited rate using rubisco_limited_rate!\nCalculate light limited rate using light_limited_rate!\nCalculate product limited rate using product_limited_rate!\nCalculate gross and net rates using colimit_photosynthesis!\nUpdate fluorescence related variables using photosystem_coefficients!\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_photosynthesis!-Union{Tuple{FT}, Tuple{Union{ClimaCache.Leaf{FT}, ClimaCache.Leaves1D{FT}, ClimaCache.Leaves2D{FT}}, ClimaCache.AirLayer{FT}, FT, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.leaf_photosynthesis!","text":"leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT) where {FT<:AbstractFloat}\n\nUpdates leaf photosynthetic rates based on CO₂ partial pressure (for StomataModels.jl temporary use), given\n\nlf Leaf, Leaves1D, or Leaves2D type structure that stores biophysical, reaction center, and photosynthesis model structures\nair AirLayer structure for environmental conditions like O₂ partial pressure\ng_lc Leaf diffusive conductance to CO₂ in [mol m⁻² s⁻¹], default is leaf._g_CO₂\nppar APAR used for photosynthesis\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.leaf_photosynthesis!-Union{Tuple{FT}, Tuple{Union{ClimaCache.Leaf{FT}, ClimaCache.Leaves1D{FT}, ClimaCache.Leaves2D{FT}}, ClimaCache.AirLayer{FT}, Union{ClimaCache.GCO₂Mode, ClimaCache.PCO₂Mode}}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.leaf_photosynthesis!","text":"leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}\n\nUpdates leaf photosynthetic rates based on CO₂ partial pressure or CO₂ conductance, given\n\nlf Leaf, Leaves1D, or Leaves2D type structure that stores biophysical, reaction center, and photosynthesis model structures\nair AirLayer structure for environmental conditions like O₂ partial pressure\nmode GCO₂Mode or PCO₂Mode that uses CO₂ conductance or partial pressure to compute photosynthetic rates\n\n\n\nExamples\n\nleaf = Leaf{Float64}(\"C3\");\nair  = AirLayer{Float64}();\nmode = PCO₂Mode();\nleaf_photosynthesis!(leaf, air, mode);\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.leaf_photosynthesis!-Union{Tuple{FT}, Tuple{ClimaCache.MonoElementSPAC{FT}, Union{ClimaCache.GCO₂Mode, ClimaCache.PCO₂Mode}}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.leaf_photosynthesis!","text":"leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}\nleaf_photosynthesis!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}\n\nUpdates leaf photosynthetic rates for SPAC, given\n\nspac MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, or MonoMLTreeSPAC type SPAC\nmode GCO₂Mode or PCO₂Mode\n\n\n\n\n\n","category":"method"},{"location":"API/#Temperature-dependency","page":"API","title":"Temperature dependency","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"temperature_correction\ntemperature_corrected_value\nphotosystem_temperature_dependence!\n∂R∂T","category":"page"},{"location":"API/#Photosynthesis.temperature_correction","page":"API","title":"Photosynthesis.temperature_correction","text":"temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}\ntemperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}\n\nReturn the correction ratio for a temperature dependent variable, given\n\ntd Arrhenius, ArrheniusPeak, or Q10 type temperature dependency struture\nt Target temperature in K\nt_ref Reference temperature in K, default is td.T_REF (298.15 K)\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.temperature_corrected_value","page":"API","title":"Photosynthesis.temperature_corrected_value","text":"temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}\n\nReturn the temperature corrected value, given\n\ntd Q10 type temperature dependency struture\nt Target temperature in K\nt_ref Reference temperature in K, default is td.T_REF (298.15 K)\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.photosystem_temperature_dependence!","page":"API","title":"Photosynthesis.photosystem_temperature_dependence!","text":"photosystem_temperature_dependence!(psm::C3CytochromeModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}\nphotosystem_temperature_dependence!(psm::C3VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}\nphotosystem_temperature_dependence!(psm::C4VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}\n\nUpdate the temperature dependencies of C3 photosynthesis model, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for photosynthesis model\nair AirLayer structure for environmental conditions like O₂ partial pressure\nt Target temperature in K\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.∂R∂T","page":"API","title":"Photosynthesis.∂R∂T","text":"∂R∂T(leaf::Leaf{FT}) where {FT<:AbstractFloat}\n∂R∂T(leaves::Leaves1D{FT}) where {FT<:AbstractFloat}\n∂R∂T(leaves::Leaves2D{FT}) where {FT<:AbstractFloat}\n\nReturn the marginal increase in respiration rate per temperature, given\n\nleaf Leaf type leaf\nleaves Leaves1D or Leaves2D type leaf\n\n\n\n\n\n","category":"function"},{"location":"API/#Electron-transport","page":"API","title":"Electron transport","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"photosystem_electron_transport!","category":"page"},{"location":"API/#Photosynthesis.photosystem_electron_transport!","page":"API","title":"Photosynthesis.photosystem_electron_transport!","text":"photosystem_electron_transport!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nphotosystem_electron_transport!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nphotosystem_electron_transport!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, apar::FT, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the electron transport rates, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel type C3 photosynthesis model\nrc CytochromeReactionCenter or VJPReactionCenter type photosynthesis system reaction center\napar Absorbed photosynthetically active radiation in μmol m⁻² s⁻¹\np_i Internal CO₂ partial pressure in Pa, used to compute etoc\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthetic-rates","page":"API","title":"Photosynthetic rates","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"rubisco_limited_rate!\nrubisco_limited_rate!(psm::Union{C3CytochromeModel{FT},C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nrubisco_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nlight_limited_rate!\nlight_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}\nlight_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nproduct_limited_rate!\nproduct_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nproduct_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}","category":"page"},{"location":"API/#Photosynthesis.rubisco_limited_rate!","page":"API","title":"Photosynthesis.rubisco_limited_rate!","text":"This function supports two types of calculations:\n\nCalculate the rate from internal CO₂\nCalculate the rate from CO₂ conductance by solving a quadratic function\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.rubisco_limited_rate!-Union{Tuple{FT}, Tuple{Union{ClimaCache.C3CytochromeModel{FT}, ClimaCache.C3VJPModel{FT}}, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.rubisco_limited_rate!","text":"rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT},C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nrubisco_limited_rate!(psm::C4VJPModel{FT}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the RubisCO limited photosynthetic rate, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for photosynthesis model\np_i Internal CO₂ partial pressure in Pa\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.rubisco_limited_rate!-Union{Tuple{FT}, Tuple{Union{ClimaCache.C3CytochromeModel{FT}, ClimaCache.C3VJPModel{FT}}, ClimaCache.AirLayer{FT}, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.rubisco_limited_rate!","text":"rubisco_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nrubisco_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the RubisCO limited photosynthetic rate in conductance mode, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for photosynthesis model\nair AirLayer structure for environmental conditions like O₂ partial pressure\ng_lc Leaf diffusive conductance to CO₂ in [mol m⁻² s⁻¹]\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.light_limited_rate!","page":"API","title":"Photosynthesis.light_limited_rate!","text":"This function supports two types of calculations:\n\nCalculate the rate from internal CO₂\nCalculate the rate from CO₂ conductance by solving a quadratic function\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.light_limited_rate!-Union{Tuple{Union{ClimaCache.C3CytochromeModel{FT}, ClimaCache.C4VJPModel{FT}}}, Tuple{FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.light_limited_rate!","text":"light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}\nlight_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat}\n\nUpdate the electron transport limited photosynthetic rate, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for C3 photosynthesis model\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.light_limited_rate!-Union{Tuple{FT}, Tuple{ClimaCache.C3CytochromeModel{FT}, ClimaCache.CytochromeReactionCenter{FT}, ClimaCache.AirLayer{FT}, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.light_limited_rate!","text":"light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nlight_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nlight_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the electron transport limited photosynthetic rate in conductance mode, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for C3 photosynthesis model\nrc CytochromeReactionCenter or VJPReactionCenter type photosynthesis system reaction center\nair AirLayer structure for environmental conditions like O₂ partial pressure\ng_lc Leaf diffusive conductance to CO₂ in [mol m⁻² s⁻¹]\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.product_limited_rate!","page":"API","title":"Photosynthesis.product_limited_rate!","text":"This function supports two types of calculations:\n\nCalculate the rate from internal CO₂\nCalculate the rate from CO₂ conductance by solving a quadratic function\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.product_limited_rate!-Union{Tuple{FT}, Tuple{Union{ClimaCache.C3CytochromeModel{FT}, ClimaCache.C3VJPModel{FT}}, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.product_limited_rate!","text":"product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nproduct_limited_rate!(psm::C4VJPModel{FT}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the product limited photosynthetic rate, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for C3 photosynthesis model\np_i Internal CO₂ partial pressure in Pa, not used in this method\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"method"},{"location":"API/#Photosynthesis.product_limited_rate!-Union{Tuple{FT}, Tuple{Union{ClimaCache.C3CytochromeModel{FT}, ClimaCache.C3VJPModel{FT}}, ClimaCache.AirLayer{FT}, FT}} where FT<:AbstractFloat","page":"API","title":"Photosynthesis.product_limited_rate!","text":"product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nproduct_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the electron transport limited photosynthetic rate in conductance mode, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel structure for C3 photosynthesis model\nair AirLayer structure for environmental conditions like O₂ partial pressure, not used in this method\ng_lc Leaf diffusive conductance to CO₂ in [mol m⁻² s⁻¹], not used in this method\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"method"},{"location":"API/#Colimitation","page":"API","title":"Colimitation","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"colimit_photosynthesis!\ncolimited_rate","category":"page"},{"location":"API/#Photosynthesis.colimit_photosynthesis!","page":"API","title":"Photosynthesis.colimit_photosynthesis!","text":"colimit_photosynthesis!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}, C4VJPModel{FT}}; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nColimit the photosynthesis by rubisco-, light-, and product-limited photosynthetic rates, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel type photosynthesis model\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd (default is 1)\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.colimited_rate","page":"API","title":"Photosynthesis.colimited_rate","text":"colimited_rate(a_1::FT, a_2::FT, colim::MinimumColimit{FT}) where {FT<:AbstractFloat}\ncolimited_rate(a_1::FT, a_2::FT, colim::QuadraticColimit{FT}) where {FT<:AbstractFloat}\ncolimited_rate(a_1::FT, a_2::FT, colim::SerialColimit{FT}) where {FT<:AbstractFloat}\n\nReturn the minimum of two rates, given\n\na_1 Rate 1\na_2 Rate 2\ncolim MinimumColimit, QuadraticColimit, or SerialColimit type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Coefficients-and-fluorescence","page":"API","title":"Coefficients and fluorescence","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"photosystem_coefficients!","category":"page"},{"location":"API/#Photosynthesis.photosystem_coefficients!","page":"API","title":"Photosynthesis.photosystem_coefficients!","text":"photosystem_coefficients!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, apar::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\nphotosystem_coefficients!(psm::Union{C3VJPModel{FT}, C4VJPModel{FT}}, rc::VJPReactionCenter{FT}, apar::FT; β::FT = FT(1)) where {FT<:AbstractFloat}\n\nUpdate the rate constants and coefficients in reaction center, given\n\npsm C3CytochromeModel, C3VJPModel, or C4VJPModel type photosynthesis model\nrc CytochromeReactionCenter or VJPReactionCenter type photosynthesis system reaction center\napar Absorbed photosynthetically active radiation in μmol m⁻² s⁻¹\nβ Tuning factor to downregulate effective Vmax, Jmax, and Rd\n\n\n\n\n\n","category":"function"}]
}