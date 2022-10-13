"""

    update_clm_vcmax_td!(td::ArrheniusPeakTD{FT}, t_opt::FT) where {FT<:AbstractFloat}

Update the TD for Vcmax, given
- `td` `ArrheniusPeakTD` type temperature dependency
- `t_opt` 10 day average temperature in `[K]`

"""
function update_clm_vcmax_td!(td::ArrheniusPeakTD{FT}, t_opt::FT) where {FT<:AbstractFloat}
    td.ΔHa_to_RT25 = 72000.0 / GAS_R() / T₂₅();
    td.ΔHd_to_R = 200000.0 / GAS_R();
    td.ΔSv_to_R = (668.39 - 1.07 * (t_opt - T₀())) / GAS_R();
    td.C = 1 + exp( td.ΔSv_to_R - td.ΔHd_to_R/T₂₅(FT) );

    return nothing
end


"""

    update_clm_vcmax_td!(td::ArrheniusPeakTD{FT}, t_opt::FT) where {FT<:AbstractFloat}

Update the TD for Jmax, given
- `td` `ArrheniusPeakTD` type temperature dependency
- `t_opt` 10 day average temperature in `[K]`

"""
function update_clm_jmax_td!(td::ArrheniusPeakTD{FT}, t_opt::FT) where {FT<:AbstractFloat}
    td.ΔHa_to_RT25 = 50000.0 / GAS_R() / T₂₅();
    td.ΔHd_to_R = 200000.0 / GAS_R();
    td.ΔSv_to_R = (659.70 - 0.75 * (t_opt - T₀())) / GAS_R();
    td.C = 1 + exp( td.ΔSv_to_R - td.ΔHd_to_R/T₂₅(FT) );

    return nothing
end
