# Settings

## Default PhotosynthesisMethods
| Field       | Default                |
|:------------|:-----------------------|
| TD_JMAX     | JmaxTDCLM              |
| TD_KC       | KcTDCLM                |
| TD_KO       | KoTDCLM                |
| TD_KPEP     | KpepTDBoyd             |
| TD_KPEP_CLM | Q10TDKpepCLM           |
| TD_KQ       | KqTDJohnson            |
| TD_R_C3     | RespirationTDCLMC3     |
| TD_R_C4     | RespirationTDCLMC4     |
| TD_VCMAX_C3 | VcmaxTDCLMC3           |
| TD_VCMAX_C4 | VcmaxTDCLMC4           |
| TD_VPMAX    | VpmaxTDBoyd            |
| TD_Γ        | ΓStarTDCLM             |
| TD_ηC       | ηCTDWang               |
| TD_ηL       | ηLTDWang               |
|||

## Rate constant
| Model                     | Prameters                                |
|:--------------------------|:-----------------------------------------|
| C3-Vqmax                  | KD1 = 0.55                               |
|                           | KF1 = 0.05                               |
|                           | KP1max = 14.5                            |
|                           | KX1 = 14.5                               |
|                           | KD2 = 0.55                               |
|                           | KF2 = 0.05                               |
|                           | KP2max = 4.5                             |
|                           | KU = 2                                   |
| C3-Jmax/C4-Vcmax/C4-Vpmax | KD = 0.95                                |
|                           | KF = 0.05                                |
|                           | KPmax = 4                                |
| For KN model              | KD = max[0.8733, 0.0301*(t-T₀())+0.0773] |
|||

## Temperature dependence
| Function               | Type             | Parameter                      |
|:---------------------- |:-----------------|:-------------------------------|
| KcTDBernacchi          | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 41.0264925           |
|                        |                  | ΔHA = 79430.0                  |
| KcTDCLM                | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 40.49                |
|                        |                  | ΔHA = 79430.0                  |
| KoTDBernacchi          | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 28208.88             |
|                        |                  | ΔHA = 36380.0                  |
| KoTDCLM                | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 27840.0              |
|                        |                  | ΔHA = 36380.0                  |
| KpepTDBoyd             | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 16.0                 |
|                        |                  | ΔHA = 36300.0                  |
| KqTDJohnson            | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 300                  |
|                        |                  | ΔHA = 37000.0                  |
| RespirationTDBernacchi | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 46390.0                  |
| VcmaxTDBernacchi       | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 65330.0                  |
| VomaxTDBernacchi       | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 60110.0                  |
| ΓStarTDBernacchi       | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 4.33164375           |
|                        |                  | ΔHA = 37830.0                  |
| ΓStarTDCLM             | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 4.275                |
|                        |                  | ΔHA = 37830.0                  |
| JmaxTDBernacchi        | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 57500.0                  |
|                        |                  | ΔHD = 439000.0                 |
|                        |                  | ΔSV = 1400.0                   |
| JmaxTDCLM              | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 50000.0                  |
|                        |                  | ΔHD = 200000.0                 |
|                        |                  | ΔSV = 659.70 - 0.75 * (t-T₀()) |
| JmaxTDLeuning          | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 50300.0                  |
|                        |                  | ΔHD = 152044.0                 |
|                        |                  | ΔSV = 495.0                    |
| RespirationTDCLMC3     | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 46390.0                  |
|                        |                  | ΔHD = 150650.0                 |
|                        |                  | ΔSV = 490.0                    |
| VcmaxTDCLMC3           | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 72000.0                  |
|                        |                  | ΔHD = 200000.0                 |
|                        |                  | ΔSV = 668.39 - 1.07 * (t-T₀()) |
| VcmaxTDLeuning         | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 73637.0                  |
|                        |                  | ΔHD = 149252.0                 |
|                        |                  | ΔSV = 486.0                    |
| VpmaxTDBoyd            | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 94800.0                  |
|                        |                  | ΔHD = 73300.0                  |
|                        |                  | ΔSV = 250.0                    |
| ηCTDJohnson            | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 1.0                  |
|                        |                  | ΔHA = 0.0                      |
|                        |                  | ΔHD = 220000.0                 |
|                        |                  | ΔSV = 710.0                    |
| ηLTDJohnson            | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 0.75                 |
|                        |                  | ΔHA = 0.0                      |
|                        |                  | ΔHD = 220000.0                 |
|                        |                  | ΔSV = 710.0                    |
| ηCTDWang               | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 1.0                  |
|                        |                  | ΔHA = 0.0                      |
|                        |                  | ΔHD = 225100.0                 |
|                        |                  | ΔSV = 710.0                    |
| ηLTDWang               | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 0.75                 |
|                        |                  | ΔHA = 0.0                      |
|                        |                  | ΔHD = 225100.0                 |
|                        |                  | ΔSV = 710.0                    |
| Q10TDAngiosperm        | Q10              | Q_10 = 1.4                     |
|                        |                  | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 2 * 0.0140 / 8.760   |
| Q10TDGymnosperm        | Q10              | Q_10 = 1.7                     |
|                        |                  | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 2 * 0.0425 / 8.760   |
| Q10TDKpepCLM           | Q10              | Q_10 = 2.0                     |
|                        |                  | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 0.2                  |
| RespirationTDCLMC4     | Q10PeakHT        | Q_10 = 2.0                     |
|                        |                  | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔT_REF = 328.15                |
|                        |                  | ΔT_SLOPE = 1.3                 |
| VcmaxTDCLMC4           | Q10PeakLTHT      | Q_10 = 2.0                     |
|                        |                  | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHT_REF = 313.15               |
|                        |                  | ΔHT_SLOPE = 0.3                |
|                        |                  | ΔLT_REF = 288.15               |
|                        |                  | ΔLT_SLOPE = 0.2                |
| **in prep**            |                  |                                |
| ΓStarTDWang2024        | Arrhenius        | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 4.56,                |
|                        |                  | ΔHA = 11800.0                  |
| JmaxTDWang2024         | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 50000                    |
|                        |                  | ΔHD = 201000                   |
|                        |                  | ΔSV = 659.70 - 0.75 * (t-T₀()) |
| KqTDWang2024           | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 300                  |
|                        |                  | ΔHA = 28500                    |
|                        |                  | ΔHD = 223500                   |
|                        |                  | ΔSV = 700                      |
| VcmaxTDWang2024        | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = NaN                  |
|                        |                  | ΔHA = 63000                    |
|                        |                  | ΔHD = 204000                   |
|                        |                  | ΔSV = 668.39 - 1.07 * (t-T₀()) |
| ηCTDWang2024           | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 4*3/14               |
|                        |                  | ΔHA = 28500                    |
|                        |                  | ΔHD = 223500                   |
|                        |                  | ΔSV = 700                      |
| ηLTDWang2024           | ArrheniusPeak    | T_REF = T₂₅()                  |
|                        |                  | VAL_REF = 3*3/14               |
|                        |                  | ΔHA = 28500                    |
|                        |                  | ΔHD = 223500                   |
|                        |                  | ΔSV = 700                      |
||||

## Smoothing function
| Function       | Type             | Parameter |
|:---------------|:-----------------|:----------|
| ColimitCJCLMC3 | QuadraticColimit | 0.98      |
| ColimitCJCLMC4 | QuadraticColimit | 0.8       |
| ColimitIPCLM   | QuadraticColimit | 0.95      |
| ColimitJCLM    | QuadraticColimit | 0.7       |
||||

## Fluorescence model
| Function                   | Type                   | Parameter            |
|:---------------------------|:-----------------------|:---------------------|
| KNFluorescenceModelAll     | KNFluorescenceModel    | K_0 = 2.48           |
|                            |                        | K_A = 2.83           |
|                            |                        | K_B = 0.114          |
| KNFluorescenceModelDrought | KNFluorescenceModel    | K_0 = 5.01           |
|                            |                        | K_A = 1.93           |
|                            |                        | K_B = 10             |
| QLFluorescenceModelC3      | QLFluorescenceModel    | K_B = 0.95e-3 / 0.85 |
| QLFluorescenceModelC4      | QLFluorescenceModel    | K_B = 0.63e-3 / 0.85 |
| QLFluorescenceModelHanC3   | QLFluorescenceModelHan | K_A = 0.8            |
|                            |                        | K_B = 0.95e-3 / 0.85 |
| QLFluorescenceModelHanC4   | QLFluorescenceModelHan | K_A = 0.83           |
|                            |                        | K_B = 0.63e-3 / 0.85 |
||||
