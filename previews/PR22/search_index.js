var documenterSearchIndex = {"docs":
[{"location":"#Photosynthesis.jl","page":"Home","title":"Photosynthesis.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Photosynthesis models for C3 and C4 photosynthesis.","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Photosynthesis\n\nmod_3 = C3CLM(Float32);\nleaf  = Leaf{Float32}();\nenvir = PM.AirLayer{Float32}();\nleaf_photo_from_pi(mod_3, leaf, envir);","category":"page"},{"location":"API/#Photosynthesis","page":"API","title":"Photosynthesis","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"CurrentModule = Photosynthesis","category":"page"},{"location":"API/#Leaf-and-Environment-Structures","page":"API","title":"Leaf and Environment Structures","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"To model photosynthesis more efficiently, we use a container (Leaf     struct) to store the photosynthesis-related information. For example, many     of the physiological parameters are temperature-dependent, but these     temperature-dependent values only need to be updated when leaf temperature     changes. Therefore, use of the container significantly reduces the time     required when programing leaf gas exchange prognostically. The     Leaf struct has the following fields:","category":"page"},{"location":"API/","page":"API","title":"API","text":"Leaf","category":"page"},{"location":"API/#Photosynthesis.Leaf","page":"API","title":"Photosynthesis.Leaf","text":"mutable struct Leaf{FT}\n\nStruct to store leaf information.\n\nFields\n\nT\nTemperature [K]\nKd\nRate constant for thermal dissipation\nKf\nRate constant for fluorescence (const)\nKn\nNPQ rate constant (initially zero)\nKp\nRate constant for photochemistry (all reaction centers open)\nmaxPSII\nmax PSII yield (Kn=0, all RC open)\nPSII_frac\nFraction of absorbed light used by PSII ETR\np_i\nLeaf internal CO₂ partial pressure [Pa]\np_s\nLeaf surface CO₂ partial pressure [Pa]\np_sat\nSaturation H₂O vapor pressure [Pa]\ng_bc\nLeaf diffusive conductance to CO₂ [mol m⁻² s⁻¹]\ng_lc\nLeaf diffusive conductance to CO₂ [mol m⁻² s⁻¹]\nAc\nRubisCO limited photosynthetic rate [μmol m⁻² s⁻¹]\nAj\nLight limited photosynthetic rate [μmol m⁻² s⁻¹]\nAg\nGross photosynthetic rate [μmol m⁻² s⁻¹]\nAn\nNet photosynthetic rate [μmol m⁻² s⁻¹]\nAp\nProduct limited photosynthetic rate [μmol m⁻² s⁻¹]\nJ\nElectron transport [μmol m⁻² s⁻¹]\nJ_pot\nPotential Electron Transport Rate [μmol m⁻² s⁻¹]\nJmax\nMaximal electron transport rate [μmol m⁻² s⁻¹]\nJmax25\nMaximal electron transport rate at 298.15 K [μmol m⁻² s⁻¹]\nKc\nRubisCO coefficient Kc [Pa]\nKo\nRubisCO coefficient Ko [Pa]\nKpep\nPEP coefficient Ko [Pa]\nKm\nMichaelis-Menten's coefficient [Pa]\nRd\nRespiration rate [μmol m⁻² s⁻¹]\nRd25\nRespiration rate at 298.15 K [μmol m⁻² s⁻¹]\nVcmax\nMaximal carboxylation rate [μmol m⁻² s⁻¹]\nVcmax25\nMaximal carboxylation rate at 298.15 K [μmol m⁻² s⁻¹]\nVpmax\nMaximal PEP carboxylation rate [μmol m⁻² s⁻¹]\nVpmax25\nMaximal PEP carboxylation rate at 298.15 K [μmol m⁻² s⁻¹]\nΓ_star\nCO₂ compensation point with the absence of Rd [Pa]\nJmax25WW\nWell watered maximal electron transport rate at 298.15 K [μmol m⁻² s⁻¹]\nRd25WW\nWell watered respiration rate at 298.15 K [μmol m⁻² s⁻¹]\nVcmax25WW\nWell watered maximal carboxylation rate at 298.15 K [μmol m⁻² s⁻¹]\nVpmax25WW\nWell watered maximal PEP carboxylation rate at 298.15 K [μmol m⁻² s⁻¹]\nCO₂_per_electron\nTotal efficiency, incl. photorespiration [mol CO₂ mol⁻¹ e-]\nFm\ndark adapted yield (Kp=0)\nFm′\nlight adapted yield (Kp=0)\nFo\ndark-adapted fluorescence yield (Kp=max)\nFo′\nlight-adapted fluorescence yield in the dark (Kp=max)\nJa\nActual electron transport rate [μmol m⁻² s⁻¹]\nNPQ\nNon-Photochemical quenching\nqQ\nPhotochemical quenching\nqE\nenergy quenching\nφ\nPSII yield\nϕs\nSteady-state (light-adapted) yield (aka Fs)\nAPAR\nAbsorbed photosynthetic active radiation [μmol m⁻² s⁻¹]\n\n\n\n\n\n","category":"type"},{"location":"API/","page":"API","title":"API","text":"Also, environmental conditions are required to compute photosynthetic rate, and     these conditions are stored in AirLayer struct. An     AirLayer struct further allows for more conveniently modeling     photosynthesis the vertical CO₂ and H₂O gradients in the canopy. The     AirLayer structs has the following fields:","category":"page"},{"location":"API/","page":"API","title":"API","text":"AirLayer","category":"page"},{"location":"API/#Photosynthesis.AirLayer","page":"API","title":"Photosynthesis.AirLayer","text":"mutable struct AirLayer{FT}\n\nStruct to store environmental conditions in each air layer corresponds to one     canopy layer.\n\nFields\n\nt_air\nAir temperature [K]\np_a\nAtmospheric CO₂ partial pressure [Pa]\np_atm\nAtmospheric pressure [Pa]\np_H₂O\nAtmospheric vapor pressure [Pa]\np_O₂\nAtmospheric O₂ partial pressure [Pa]\np_sat\nSaturation vapor pressure [Pa]\nRH\nRelatiev humidity\nvpd\nVapor pressure deficit [Pa]\nwind\nWind speed [m s⁻¹]\n\n\n\n\n\n","category":"type"},{"location":"API/","page":"API","title":"API","text":"See exmaples below for how to create the structs","category":"page"},{"location":"API/","page":"API","title":"API","text":"using Photosynthesis\n\nFT = Float32;\nleaf = Leaf{FT}();\nenvir = AirLayer{FT}();","category":"page"},{"location":"API/#Temperature-Dependency-Structs","page":"API","title":"Temperature Dependency Structs","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"The temperature-dependent (TD) photosynthetic parameters include","category":"page"},{"location":"API/","page":"API","title":"API","text":"J_textmax Maximal electron transport rate\nK_textc Michaelis constant for CO₂\nK_textm Michaelis-Menten's coefficient\nK_texto Michaelis constant for O₂\nK_textpep Michaelis constant for PEP carboxylation\nR_textd Dark respiration\nV_textcmax Maximal RuBP carboxylation rate\nV_textomax Maximal RuBP oxygenation rate\nV_textpmax Maximal PEP carboxylation rate\nΓ^* CO₂ compensation point with the absence of dark respiration","category":"page"},{"location":"API/","page":"API","title":"API","text":"There are two typical types of temperature dependencies using the classic     Arrhenius equation. We define the three types as ArrheniusTD,     ArrheniusPeakTD, and Q10TD subject to     AbstractTDParameterSet type:","category":"page"},{"location":"API/","page":"API","title":"API","text":"AbstractTDParameterSet\nArrheniusTD\nArrheniusPeakTD\nQ10TD","category":"page"},{"location":"API/#Photosynthesis.AbstractTDParameterSet","page":"API","title":"Photosynthesis.AbstractTDParameterSet","text":"abstract type AbstractTDParameterSet{FT}\n\nHierarchy of the AbstractTDParameterSet:\n\nArrheniusTD\nArrheniusPeakTD\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.ArrheniusTD","page":"API","title":"Photosynthesis.ArrheniusTD","text":"struct ArrheniusTD{FT}\n\nAn AbstractTDParameterSet type struct using\n\ncorr = exp left( dfracΔHaR T_0 - dfracΔHaR T_1 right)\n\nFields\n\nVAL_25\nUncorrected value at 298.15 K\nΔHa_to_R\nRatio between ΔHa and R [K]\nΔHa_to_RT25\nRatio between ΔHa and R*K_25\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.ArrheniusPeakTD","page":"API","title":"Photosynthesis.ArrheniusPeakTD","text":"struct ArrheniusPeakTD{FT}\n\nAn AbstractTDParameterSet type struct using\n\ncorr = exp left( dfracΔHaR T_0 - dfracΔHaR T_1 right)\n       cdot\n       dfrac 1 + exp left( dfracS_v T_0 - H_dR T_0 right) \n               1 + exp left( dfracS_v T_1 - H_dR T_1 right) \n\nFields\n\nΔHa_to_RT25\nRatio between ΔHa and R*K_25\nΔHd_to_R\nRatio between ΔHd and R\nΔSv_to_R\nRatio between ΔSv and R\nC\nCorrection factor C = 1 + exp( Sv/R + Hd/(RT0) )\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.Q10TD","page":"API","title":"Photosynthesis.Q10TD","text":"struct Q10TD{FT}\n\nAn AbstractTDParameterSet type struct using\n\nVAL = VAL_REF left( dfracT_1 - T_REF10 right)^Q_10\n\nFields\n\nVAL_REF\nUncorrected value at reference temperature\nT_REF\nReference temperature [K]\nQ_10\nPower of Q10 correction\n\n\n\n\n\n","category":"type"},{"location":"API/","page":"API","title":"API","text":"There are many published parameter sets for the various temperature     dependencies, and to ease the modeling we predefined most of the structs:","category":"page"},{"location":"API/","page":"API","title":"API","text":"JmaxTDBernacchi\nJmaxTDCLM\nJmaxTDLeuning\nKcTDBernacchi\nKcTDCLM\nKoTDBernacchi\nKoTDCLM\nKpepTDBoyd\nKpepTDCLM\nRespirationTDBernacchi\nRespirationTDCLM\nVcmaxTDBernacchi\nVcmaxTDCLM\nVcmaxTDLeuning\nVomaxTDBernacchi\nVpmaxTDBoyd\nΓStarTDBernacchi\nΓStarTDCLM","category":"page"},{"location":"API/#Photosynthesis.JmaxTDBernacchi","page":"API","title":"Photosynthesis.JmaxTDBernacchi","text":"ArrheniusPeakTD type Jmax TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.JmaxTDCLM","page":"API","title":"Photosynthesis.JmaxTDCLM","text":"ArrheniusPeakTD type Jmax TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.JmaxTDLeuning","page":"API","title":"Photosynthesis.JmaxTDLeuning","text":"ArrheniusPeakTD type Jmax TD from Leuning's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KcTDBernacchi","page":"API","title":"Photosynthesis.KcTDBernacchi","text":"ArrheniusTD type Kc TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KcTDCLM","page":"API","title":"Photosynthesis.KcTDCLM","text":"ArrheniusTD type Kc TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KoTDBernacchi","page":"API","title":"Photosynthesis.KoTDBernacchi","text":"ArrheniusTD type Ko TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KoTDCLM","page":"API","title":"Photosynthesis.KoTDCLM","text":"ArrheniusTD type Ko TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KpepTDBoyd","page":"API","title":"Photosynthesis.KpepTDBoyd","text":"ArrheniusTD type Kpep TD from Boyd's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.KpepTDCLM","page":"API","title":"Photosynthesis.KpepTDCLM","text":"ArrheniusTD type Kpep TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.RespirationTDBernacchi","page":"API","title":"Photosynthesis.RespirationTDBernacchi","text":"ArrheniusTD type Respiration TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.RespirationTDCLM","page":"API","title":"Photosynthesis.RespirationTDCLM","text":"ArrheniusPeakTD type Respiration TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VcmaxTDBernacchi","page":"API","title":"Photosynthesis.VcmaxTDBernacchi","text":"ArrheniusTD type Vcmax TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VcmaxTDCLM","page":"API","title":"Photosynthesis.VcmaxTDCLM","text":"ArrheniusPeakTD type Vcmax TD \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VcmaxTDLeuning","page":"API","title":"Photosynthesis.VcmaxTDLeuning","text":"ArrheniusPeakTD type Vcmax TD from Leuning's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VomaxTDBernacchi","page":"API","title":"Photosynthesis.VomaxTDBernacchi","text":"ArrheniusTD type Vomax TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VpmaxTDBoyd","page":"API","title":"Photosynthesis.VpmaxTDBoyd","text":"ArrheniusPeakTD type Vpmax TD from Boyd's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.ΓStarTDBernacchi","page":"API","title":"Photosynthesis.ΓStarTDBernacchi","text":"ArrheniusTD type Γ^* TD from Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.ΓStarTDCLM","page":"API","title":"Photosynthesis.ΓStarTDCLM","text":"ArrheniusTD type Γ* TD \n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"The TDs can be easily created using commands like","category":"page"},{"location":"API/","page":"API","title":"API","text":"using Photosynthesis\n\nFT = Float32;\n_td_1 = JmaxTDBernacchi(FT);\n_td_2 = VcmaxTDCLM(FT);","category":"page"},{"location":"API/","page":"API","title":"API","text":"However, be aware that these pre-defined TD structs are not mutable, to create     customized TD struct, code like this will be useful","category":"page"},{"location":"API/","page":"API","title":"API","text":"using Photosynthesis\n\nFT = Float32;\n_td_1 = ArrheniusTD{FT}(1, 10000, 30);\n_td_1 = ArrheniusPeakTD{FT}(1, 10000, 30, 1);","category":"page"},{"location":"API/","page":"API","title":"API","text":"To further simplify the use of Photosynthesis module, we provide a few     collections/structs of temperature dependencies as well as other parameter     sets like FluoParaSet. The structs are catergorized to     C3ParaSet and C4ParaSet subject to an     AbstractPhotoModelParaSet type, and the structs are meant for     modeling C3 photosynthesis and C4 photosynthesis, respectively.","category":"page"},{"location":"API/","page":"API","title":"API","text":"AbstractPhotoModelParaSet\nC3ParaSet\nC4ParaSet","category":"page"},{"location":"API/#Photosynthesis.AbstractPhotoModelParaSet","page":"API","title":"Photosynthesis.AbstractPhotoModelParaSet","text":"abstract type AbstractPhotoModelParaSet{FT}\n\nHierarchy of the AbstractPhotoModelParaSet:\n\nC3ParaSet\nC4ParaSet\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.C3ParaSet","page":"API","title":"Photosynthesis.C3ParaSet","text":"mutable struct C3Paraset{FT}\n\nParameter sets for C3 photosynthesis.\n\nFields\n\nJT\nJmax temperature dependency\nKcT\nKc temperature dependency\nKoT\nKo temperature dependency\nReT\nRespiration temperature dependency\nVcT\nVcmax temperature dependency\nΓsT\nΓ_star temperature dependency\nFlu\nFluorescence model\nVR\nVcmax25 and respiration correlation\nEff_1\nCoefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively\nEff_2\nCoefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.C4ParaSet","page":"API","title":"Photosynthesis.C4ParaSet","text":"mutable struct C4ParaSet{FT}\n\nParameter sets for C3 photosynthesis.\n\nFields\n\nKpT\nKpep temperature dependency\nReT\nRespiration temperature dependency\nVcT\nVcmax temperature dependency\nVpT\nVpmax temperature dependency\nFlu\nFluorescence model\nVR\nVcmax25 and respiration correlation\n\n\n\n\n\n","category":"type"},{"location":"API/","page":"API","title":"API","text":"Again, to guarantee a quick start, we provided a few pre-defined parameter     sets:","category":"page"},{"location":"API/","page":"API","title":"API","text":"C3Bernacchi\nC3CLM\nC4CLM","category":"page"},{"location":"API/#Photosynthesis.C3Bernacchi","page":"API","title":"Photosynthesis.C3Bernacchi","text":"C3ParaSet type C3 photosynthesis using Bernacchi's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.C3CLM","page":"API","title":"Photosynthesis.C3CLM","text":"C3ParaSet type C3 photosynthesis using CLM5's data \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.C4CLM","page":"API","title":"Photosynthesis.C4CLM","text":"C4ParaSet type C4 photosynthesis using CLM5's data \n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Examples:","category":"page"},{"location":"API/","page":"API","title":"API","text":"using Photosynthesis\n\nFT = Float32;\nset_b = C3Bernacchi(FT);\nset_3 = C3CLM(FT);\nset_4 = C4CLM(FT);","category":"page"},{"location":"API/","page":"API","title":"API","text":"Note it here that the C3ParaSet and C4ParaSet structs are     mutable, and the fields can be changed to another non-mutable TD struct.     We'd like to mention that in some cases, leaf respiration rate is not     measured, and in this case, the dark respiration rate will be computed from     V_textcmax using a multiplier","category":"page"},{"location":"API/","page":"API","title":"API","text":"VtoRCollatz\nVtoRDefault","category":"page"},{"location":"API/#Photosynthesis.VtoRCollatz","page":"API","title":"Photosynthesis.VtoRCollatz","text":"A constant of 0.01 \n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.VtoRDefault","page":"API","title":"Photosynthesis.VtoRDefault","text":"A constant of 0.015 \n\n\n\n\n\n","category":"function"},{"location":"API/#Temperature-Dependency","page":"API","title":"Temperature Dependency","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"As mentioned above, temperature corrections only need to be done once per     temperature change, and storing the temperature corrected values will     significantly boost the code speed. Here we provide a few functions to     change the stored values. First of all, all the temperature corrections are     made with temperature_correction:","category":"page"},{"location":"API/","page":"API","title":"API","text":"temperature_correction","category":"page"},{"location":"API/#Photosynthesis.temperature_correction","page":"API","title":"Photosynthesis.temperature_correction","text":"temperature_correction(td_set::AbstractTDParameterSet{FT},\n                       T::FT) where {FT<:AbstractFloat}\n\nA correction factor based on arrhenius's fitting procedure, given\n\ntd_set ArrheniusTD or ArrheniusPeakTD type struct\nT Leaf temperature in [K]\n\nThe equation used for ArrheniusTD is\n\ncorr = exp left( dfracΔHaR T_0 - dfracΔHaR T_1 right)\n\nThe equations used for ArrheniusPeakTD are\n\ncorr = exp left( dfracΔHaR T_0 - dfracΔHaR T_1 right)\n       cdot\n       dfrac 1 + exp left( dfracS_v T_0 - H_dR T_0 right) \n               1 + exp left( dfracS_v T_1 - H_dR T_1 right) \n\nThe equation used for Q10TD is\n\ncorr = left( dfracT_1 - T_REF10 right)^Q_10\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Second, depending on which physiological parameter to correct, some corrections     use the VAL_25 field in the ArrheniusTD, like K_textc,     K_texto, and K_textpep:","category":"page"},{"location":"API/","page":"API","title":"API","text":"photo_TD_from_set","category":"page"},{"location":"API/#Photosynthesis.photo_TD_from_set","page":"API","title":"Photosynthesis.photo_TD_from_set","text":"photo_TD_from_set(td_set::ArrheniusTD{FT}, T::FT) where {FT<:AbstractFloat}\n\nMake temperature correction from parameter set, given\n\ntd_set ArrheniusTD type parameter set, which has a VAL_25 field\nT Leaf temperature\n\nUseful for Kc, Ko, Kpep, and Γ^*.\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Some corrections use the reference values from the Leaf struct, like     V_textcmax and J_textmax:","category":"page"},{"location":"API/","page":"API","title":"API","text":"photo_TD_from_val","category":"page"},{"location":"API/#Photosynthesis.photo_TD_from_val","page":"API","title":"Photosynthesis.photo_TD_from_val","text":"photo_TD_from_val(td_set::AbstractTDParameterSet{FT}, val::FT, T::FT) where {FT<:AbstractFloat}\n\nMake temperature correction from a given value, given\n\ntd_set ArrheniusTD or ArrheniusPeakTD type struct\nval Uncorrected value at 298.15 K\nT Leaf temperature\n\nUseful for Vcmax, Vomac, Vpmax, Jmax, and Respiration.\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"The functions to make temperature corrections to each individual variables are","category":"page"},{"location":"API/","page":"API","title":"API","text":"leaf_jmax!\nleaf_kc!\nleaf_km!\nleaf_ko!\nleaf_kpep!\nleaf_rd!\nleaf_vcmax!\nleaf_vpmax!\nleaf_Γstar!","category":"page"},{"location":"API/#Photosynthesis.leaf_jmax!","page":"API","title":"Photosynthesis.leaf_jmax!","text":"leaf_jmax!(td_set::AbstractTDParameterSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate maximal electron transport rate at leaf temperature, given\n\ntd_set AbstractTDParameterSet type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_kc!","page":"API","title":"Photosynthesis.leaf_kc!","text":"leaf_kc!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate Kc at leaf temperature, given\n\ntd_set ArrheniusTD type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_km!","page":"API","title":"Photosynthesis.leaf_km!","text":"leaf_km!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\n\nUpdate Ko at leaf temperature, given\n\nphoto_set C3ParaSet type photosynthesis parameter set\nleaf Leaf type struct\nenvir AirLayer type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_ko!","page":"API","title":"Photosynthesis.leaf_ko!","text":"leaf_ko!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate Ko at leaf temperature, given\n\ntd_set ArrheniusTD type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_kpep!","page":"API","title":"Photosynthesis.leaf_kpep!","text":"leaf_kpep!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate Kpep at leaf temperature, given\n\ntd_set ArrheniusTD type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_rd!","page":"API","title":"Photosynthesis.leaf_rd!","text":"leaf_rd!(td_set::AbstractTDParameterSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate leaf dark respiration rate at leaf temperature, given\n\ntd_set AbstractTDParameterSet type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_vcmax!","page":"API","title":"Photosynthesis.leaf_vcmax!","text":"leaf_vcmax!(td_set::AbstractTDParameterSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate leaf maximal carboxylation rate at leaf temperature, given\n\ntd_set AbstractTDParameterSet type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_vpmax!","page":"API","title":"Photosynthesis.leaf_vpmax!","text":"leaf_vpmax!(td_set::AbstractTDParameterSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate leaf maximal PEP carboxylation rate at leaf temperature, given\n\ntd_set AbstractTDParameterSet type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_Γstar!","page":"API","title":"Photosynthesis.leaf_Γstar!","text":"leaf_Γstar!(td_set::ArrheniusTD{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate Γ^* at leaf temperature, given\n\ntd_set ArrheniusTD type TD parameter set\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Again to ease the coding, we provide a function to run all the temperature     dependencies:","category":"page"},{"location":"API/","page":"API","title":"API","text":"leaf_temperature_dependence!","category":"page"},{"location":"API/#Photosynthesis.leaf_temperature_dependence!","page":"API","title":"Photosynthesis.leaf_temperature_dependence!","text":"leaf_temperature_dependence!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\nleaf_temperature_dependence!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\nleaf_temperature_dependence!(photo_set::AbstractPhotoModelParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}, T::FT) where {FT<:AbstractFloat}\n\nUpdate the temperature dependent photosynthesis only, given\n\nphoto_set AbstractPhotoModelParaSet type parameter set\nleaf Leaf type struct\nenvir AirLayer type struct\nT Given leaf temperature\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Note it here that function leaf_temperature_dependence! updates     saturated vapor pressure from leaf temperature as well.","category":"page"},{"location":"API/","page":"API","title":"API","text":"Example:","category":"page"},{"location":"API/","page":"API","title":"API","text":"using Photosynthesis\n\nFT = Float32;\nleaf = Leaf{FT}();\nenvir = AirLayer{FT}();\nset_3 = C3CLM(FT);\n\nleaf_temperature_dependence!(c3_set, leaf, envir);\nleaf_temperature_dependence!(c3_set, leaf, envir, FT(300));","category":"page"},{"location":"API/#RubisCO-limited-Photosynthesis","page":"API","title":"RubisCO-limited Photosynthesis","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"By default, Photosynthesis module computes gross photosynthetic rate as the     minimal of the three:","category":"page"},{"location":"API/","page":"API","title":"API","text":"A_textc RubisCO-limited photosynthetic rate\nA_textj Light-limited photosynthetic rate\nA_textp Product-limited photosynthetic rate","category":"page"},{"location":"API/","page":"API","title":"API","text":"If leaf internal CO₂ is known, A_textc (gross rate) can be computed using","category":"page"},{"location":"API/","page":"API","title":"API","text":"rubisco_limited_rate!","category":"page"},{"location":"API/#Photosynthesis.rubisco_limited_rate!","page":"API","title":"Photosynthesis.rubisco_limited_rate!","text":"rubisco_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\nrubisco_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nCalculate the RubisCO limited photosynthetic rate, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"If total leaf diffusive conductance to CO₂ is known, A_textc can be     computed analytically by solving the quadratic function using","category":"page"},{"location":"API/","page":"API","title":"API","text":"lower_quadratic","category":"page"},{"location":"API/#Photosynthesis.lower_quadratic","page":"API","title":"Photosynthesis.lower_quadratic","text":"lower_quadratic(a::FT, b::FT, c::FT) where {FT<:AbstractFloat}\n\nReturn the lower quadratic solution or NaN, given\n\na Parameter in a*x^2 + b*x + c = 0\nb Parameter in a*x^2 + b*x + c = 0\nc Parameter in a*x^2 + b*x + c = 0\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"The function to analytically compute A_textc is","category":"page"},{"location":"API/","page":"API","title":"API","text":"rubisco_limited_rate_glc!","category":"page"},{"location":"API/#Photosynthesis.rubisco_limited_rate_glc!","page":"API","title":"Photosynthesis.rubisco_limited_rate_glc!","text":"rubisco_limited_rate_glc!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\n\nCalculate the RubisCO limited photosynthetic rate from glc, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\nenvir AirLayer type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Note it here that rubisco_limited_rate_glc! only applies to C3     photosynthesis as the RubisCO-limited rate for C4 plants is     V_textcmax.","category":"page"},{"location":"API/#Light-limited-Photosynthesis","page":"API","title":"Light-limited Photosynthesis","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"If leaf internal CO₂ is known, A_textj (gross rate) can be computed using","category":"page"},{"location":"API/","page":"API","title":"API","text":"light_limited_rate!","category":"page"},{"location":"API/#Photosynthesis.light_limited_rate!","page":"API","title":"Photosynthesis.light_limited_rate!","text":"light_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\nlight_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nCalculate the Light limited photosynthetic rate, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"If total leaf diffusive conductance to CO₂ is known, A_textj can be     computed analytically using","category":"page"},{"location":"API/","page":"API","title":"API","text":"light_limited_rate_glc!","category":"page"},{"location":"API/#Photosynthesis.light_limited_rate_glc!","page":"API","title":"Photosynthesis.light_limited_rate_glc!","text":"light_limited_rate_glc!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\n\nCalculate the Light limited photosynthetic rate from glc, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\nenvir AirLayer type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Note it here that light_limited_rate_glc! only applies to C3     photosynthesis as the RubisCO-limited rate for C4 plants is the electron     transport rate.","category":"page"},{"location":"API/","page":"API","title":"API","text":"Be aware that leaf electron transport rate needs to be calculated before the     light-limited rate:","category":"page"},{"location":"API/","page":"API","title":"API","text":"leaf_ETR!","category":"page"},{"location":"API/#Photosynthesis.leaf_ETR!","page":"API","title":"Photosynthesis.leaf_ETR!","text":"leaf_ETR!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\nleaf_ETR!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nUpdate the electron transport variables in the leaf struct, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/#Product-limited-Photosynthesis","page":"API","title":"Product-limited Photosynthesis","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"If leaf internal CO₂ is known, A_textp (gross rate) can be computed using","category":"page"},{"location":"API/","page":"API","title":"API","text":"product_limited_rate!","category":"page"},{"location":"API/#Photosynthesis.product_limited_rate!","page":"API","title":"Photosynthesis.product_limited_rate!","text":"product_limited_rate!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\nproduct_limited_rate!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nCalculate the Product limited photosynthetic rate, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"If total leaf diffusive conductance to CO₂ is known, A_textp can be     computed analytically using","category":"page"},{"location":"API/","page":"API","title":"API","text":"product_limited_rate_glc!","category":"page"},{"location":"API/#Photosynthesis.product_limited_rate_glc!","page":"API","title":"Photosynthesis.product_limited_rate_glc!","text":"product_limited_rate_glc!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\n\nCalculate the Product limited photosynthetic rate from glc, given\n\nphoto_set C3ParaSet or C4ParaSet type struct\nleaf Leaf type struct\nenvir AirLayer type struct\n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"Note it here that product_limited_rate_glc! only applies to C4     photosynthesis as the RubisCO-limited rate for C4 plants is     V_textcmax/2.","category":"page"},{"location":"API/#Photosynthetic-Rates","page":"API","title":"Photosynthetic Rates","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"For empirical and optimization stomatal models, iterations are required to get     the final solution as in StomataModels module. In this case, more     conveniently computing photosynthetic rates for each leaf is preferable. In     this case, leaf_photo_from_pi! and leaf_photo_from_glc!     are better options:","category":"page"},{"location":"API/","page":"API","title":"API","text":"leaf_photo_from_pi!\nleaf_photo_from_glc!","category":"page"},{"location":"API/#Photosynthesis.leaf_photo_from_pi!","page":"API","title":"Photosynthesis.leaf_photo_from_pi!","text":"leaf_photo_from_pi!(photo_set::AbstractPhotoModelParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\nleaf_photo_from_pi!(photo_set::AbstractPhotoModelParaSet{FT}, leaf::Leaf{FT}, p_i::FT) where {FT<:AbstractFloat}\n\nCompute leaf photosynthetic rates, given\n\nphoto_set AbstractPhotoModelParaSet type parameter set\nleaf Leaf type struct\np_i Given leaf internal CO₂\n\nThe C3 photosynthesis model is from Farquhar et al. (1980) \"A biochemical model     of photosynthetic CO₂ assimilation in leaves of C3 species.\"\n\nThe C4 photosynthesis model is adapted from Collatz et al. (1992) \"Coupled     photosynthesis-stomatal conductance model for leaves of C4 plants.\"\n\n\n\n\n\n","category":"function"},{"location":"API/#Photosynthesis.leaf_photo_from_glc!","page":"API","title":"Photosynthesis.leaf_photo_from_glc!","text":"leaf_photo_from_glc!(photo_set::C3ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\nleaf_photo_from_glc!(photo_set::C4ParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT}) where {FT<:AbstractFloat}\nleaf_photo_from_glc!(photo_set::AbstractPhotoModelParaSet{FT}, leaf::Leaf{FT}, envir::AirLayer{FT},g_lc::FT) where {FT<:AbstractFloat}\n\nUpdate leaf photosynthetic rates from a known leaf diffusive conductance, given\n\nphoto_set AbstractPhotoModelParaSet type parameter set\nleaf Leaf type struct\nenvir AirLayer type struct\ng_lc Given leaf diffusive conductance to CO₂\n\n\n\n\n\n","category":"function"},{"location":"API/#Fluorescence","page":"API","title":"Fluorescence","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Photosynthesis module also provide ways to compute leaf fluorescence. By     default, the modules uses fluorescence parameters from Flexas et al. with     struct FluorescenceFlexas:","category":"page"},{"location":"API/","page":"API","title":"API","text":"AbstractFluoModelParaSet\nFluoParaSet\nFluorescenceFlexas","category":"page"},{"location":"API/#Photosynthesis.AbstractFluoModelParaSet","page":"API","title":"Photosynthesis.AbstractFluoModelParaSet","text":"abstract type AbstractFluoModelParaSet{FT}\n\nHierarchy of the AbstractFluoModelParaSet:\n\nFluoParaSet\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.FluoParaSet","page":"API","title":"Photosynthesis.FluoParaSet","text":"mutable struct FluoParaSet{FT}\n\nA AbstractFluoModelParaSet type paramter set.\n\nFields\n\nKn1\nFluorescence model coefficient\nKn2\nFluorescence model coefficient\nKn3\nFluorescence model coefficient\n\n\n\n\n\n","category":"type"},{"location":"API/#Photosynthesis.FluorescenceFlexas","page":"API","title":"Photosynthesis.FluorescenceFlexas","text":"FluoParaSet type parameter set using Flexas's data \n\n\n\n\n\n","category":"function"},{"location":"API/","page":"API","title":"API","text":"The function that is used to compute fluorescene is","category":"page"},{"location":"API/","page":"API","title":"API","text":"leaf_fluorescence!","category":"page"},{"location":"API/#Photosynthesis.leaf_fluorescence!","page":"API","title":"Photosynthesis.leaf_fluorescence!","text":"leaf_fluorescence!(fluo_set::FluoParaSet{FT}, leaf::Leaf{FT}) where {FT<:AbstractFloat}\n\nCompute fluorescence yield, Kn, and Kp for leaf, given\n\nfluo_set FluoParaSet type parameter set\nleaf Leaf struct\n\n\n\n\n\n","category":"function"}]
}
