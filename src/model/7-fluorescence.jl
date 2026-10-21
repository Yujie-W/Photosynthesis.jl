# To do
#     TODO: add more calculations such as NPQ when the model is ready
"""

    photosystem_coefficients!(
                methods::PhotosynthesisMethods{FT},
                ps::LeafPhotosystem{FT},
                ppar::Vector{FT};
                β::FT = FT(1)
        ) where {FT}

Update the rate constants and coefficients in reaction center, given
- `methods` `PhotosynthesisMethods` structure that contains rate constants, temperature dependency methods, colimitation methods, and model features
- `ps` `LeafPhotosystem` type struct
- `ppar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
function photosystem_coefficients! end;

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            ps::LeafPhotosystem{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(methods, ps.trait, ps.state, ps.auxil, ppar; β = β);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(methods, pst, pss, psa, methods.FLUORESCENCE_METHOD_C3, ppar; β = β);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = photosystem_coefficients!(methods, pst, pss, psa, methods.FLUORESCENCE_METHOD_C4, ppar; β = β);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::C3Trait{FT},
            pss::C3State{FT},
            psa::LeafPhotosystemAuxil{FT},
            flm::CytochromeFluorescenceModel,
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PSI_RATE_CONSTANTS, PSII_RATE_CONSTANTS) = methods;

    # unpack the vars from the cache
    ϕ_P1_a   = psa.cache1;
    ϕ_P2_a   = psa.cache2;
    q1       = psa.cache3;
    q2       = psa.cache4;
    k_sum_nb = psa.cache5;
    k_sum_nc = psa.cache6;
    k_sum    = psa.cache7;
    k_sum_1  = psa.cache8;
    k_sum_2  = psa.cache9;
    ϕ_U2_a   = psa.cache0;

    # adapted from https://github.com/jenjohnson/johnson-berry-2021-pres/blob/main/scripts/model_fun.m
    @. ϕ_P1_a = psa.a_g * psa.η / (psa.e2c * ppar * (1 - psa.f_psii));
    @. ϕ_P2_a = psa.a_g / (psa.e2c * ppar * psa.f_psii);
    @. q1 = ϕ_P1_a / psa.ϕ_psi_max;
    @. q2 = 1 - psa.j_psi / (β * psa.v_qmax);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    @. k_sum_nb = -1 * (PSII_RATE_CONSTANTS.K_U * ϕ_P2_a + PSII_RATE_CONSTANTS.K_P * (q2 - ϕ_P2_a));
    @. k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * PSII_RATE_CONSTANTS.K_U * PSII_RATE_CONSTANTS.K_P);
    @. k_sum = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    @. psa.k_n = k_sum - PSII_RATE_CONSTANTS.K_F - PSII_RATE_CONSTANTS.K_U - PSII_RATE_CONSTANTS.K_D;

    # compute PSII and PSI yeilds
    @. k_sum_1 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n;
    @. k_sum_2 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n + PSII_RATE_CONSTANTS.K_P;
    k_sum_3 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P;
    k_sum_4 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_X;
    @. ϕ_U2_a = (q2 * PSII_RATE_CONSTANTS.K_U / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_U / k_sum_1);
    @. psa.ϕ_f2 = (q2 * PSII_RATE_CONSTANTS.K_F / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_F / k_sum_1) / (1 - ϕ_U2_a);
    @. psa.ϕ_f1 = PSI_RATE_CONSTANTS.K_F / k_sum_3 * q1 + PSI_RATE_CONSTANTS.K_F / k_sum_4 * (1 - q1);
    _ϕ_d1 = @. PSI_RATE_CONSTANTS.K_D / PSI_RATE_CONSTANTS.K_F * psa.ϕ_f1;
    _ϕ_d2 = @. PSII_RATE_CONSTANTS.K_D / PSII_RATE_CONSTANTS.K_F * psa.ϕ_f2;
    _ϕ_n1 = @. 1 - psa.ϕ_f1 - _ϕ_d1 - ϕ_P1_a;
    _ϕ_n2 = @. 1 - psa.ϕ_f2 - _ϕ_d2 - ϕ_P2_a;

    # save the weighted fluorescence and photosynthesis yields in reaction center
    @. psa.ϕ_d = _ϕ_d1 * (1 - psa.f_psii) + _ϕ_d2 * psa.f_psii;
    @. psa.ϕ_f = psa.ϕ_f1 * (1 - psa.f_psii) + psa.ϕ_f2 * psa.f_psii;
    @. psa.ϕ_n = _ϕ_n1 * (1 - psa.f_psii) + _ϕ_n2 * psa.f_psii;
    @. psa.ϕ_p = ϕ_P1_a * (1 - psa.f_psii) + ϕ_P2_a * psa.f_psii;

    if any(isnan.(psa.ϕ_f)) || any(isnan.(psa.ϕ_p))
        println();
        @show psa.a_g psa.η psa.e2c ppar psa.f_psii;
        error("NaN detected in photosystem coefficients calculation");
    end;

    return nothing
);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::C4Trait{FT},
            pss::C4State{FT},
            psa::LeafPhotosystemAuxil{FT},
            flm::CytochromeFluorescenceModel,
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PSI_RATE_CONSTANTS, PSII_RATE_CONSTANTS) = methods;

    # unpack the vars from the cache
    ϕ_P1_a   = psa.cache1;
    ϕ_P2_a   = psa.cache2;
    q1       = psa.cache3;
    q2       = psa.cache4;
    k_sum_nb = psa.cache5;
    k_sum_nc = psa.cache6;
    k_sum    = psa.cache7;
    k_sum_1  = psa.cache8;
    k_sum_2  = psa.cache9;
    ϕ_U2_a   = psa.cache0;

    # adapted from the C3 model

    #
    #
    # for the mesophyll
    #
    #
    @. ϕ_P1_a = psa.a_p * psa.η_m / (FT(0.5) * ppar * (1 - pst.f_b6f_bs) * (1 - psa.f_psii_m));
    @. ϕ_P2_a = psa.a_p / (FT(0.5) * ppar * (1 - pst.f_b6f_bs) * psa.f_psii_m);
    @. q1 = ϕ_P1_a / psa.ϕ_psi_max;
    @. q2 = 1 - psa.j_psi_m / (β * psa.v_qmax_m);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    @. k_sum_nb = -1 * (PSII_RATE_CONSTANTS.K_U * ϕ_P2_a + PSII_RATE_CONSTANTS.K_P * (q2 - ϕ_P2_a));
    @. k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * PSII_RATE_CONSTANTS.K_U * PSII_RATE_CONSTANTS.K_P);
    @. k_sum = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    @. psa.k_n_m = k_sum - PSII_RATE_CONSTANTS.K_F - PSII_RATE_CONSTANTS.K_U - PSII_RATE_CONSTANTS.K_D;

    # compute PSII and PSI yeilds
    @. k_sum_1 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n_m;
    @. k_sum_2 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n_m + PSII_RATE_CONSTANTS.K_P;
    k_sum_3 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P;
    k_sum_4 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_X;
    @. ϕ_U2_a = (q2 * PSII_RATE_CONSTANTS.K_U / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_U / k_sum_1);
    @. psa.ϕ_f2_m = (q2 * PSII_RATE_CONSTANTS.K_F / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_F / k_sum_1) / (1 - ϕ_U2_a);
    @. psa.ϕ_f1_m = PSI_RATE_CONSTANTS.K_F / k_sum_3 * q1 + PSI_RATE_CONSTANTS.K_F / k_sum_4 * (1 - q1);
    @. psa.ϕ_p1_m = ϕ_P1_a;
    @. psa.ϕ_p2_m = ϕ_P2_a;
    @. psa.ϕ_d1_m = PSI_RATE_CONSTANTS.K_D / PSI_RATE_CONSTANTS.K_F * psa.ϕ_f1_m;
    @. psa.ϕ_d2_m = PSII_RATE_CONSTANTS.K_D / PSII_RATE_CONSTANTS.K_F * psa.ϕ_f2_m;
    @. psa.ϕ_n1_m = 1 - psa.ϕ_f1_m - psa.ϕ_d1_m - psa.ϕ_p1_m;
    @. psa.ϕ_n2_m = 1 - psa.ϕ_f2_m - psa.ϕ_d2_m - psa.ϕ_p2_m;

    #
    #
    # for bundle sheath, it is the same as the C3 model
    #
    #
    @. ϕ_P1_a = psa.a_g * psa.η_bs / (psa.e2c * ppar * pst.f_b6f_bs * (1 - psa.f_psii_bs));
    @. ϕ_P2_a = psa.a_g / (psa.e2c * ppar * pst.f_b6f_bs * psa.f_psii_bs);
    @. q1 = ϕ_P1_a / psa.ϕ_psi_max;
    @. q2 = 1 - psa.j_psi_bs / (β * psa.v_qmax_bs);

    # solve PSII K_N
    k_sum_na = ϕ_P2_a;
    @. k_sum_nb = -1 * (PSII_RATE_CONSTANTS.K_U * ϕ_P2_a + PSII_RATE_CONSTANTS.K_P * (q2 - ϕ_P2_a));
    @. k_sum_nc = -1 * (ϕ_P2_a * (1 - q2) * PSII_RATE_CONSTANTS.K_U * PSII_RATE_CONSTANTS.K_P);
    @. k_sum = upper_quadratic(k_sum_na, k_sum_nb, k_sum_nc);
    @. psa.k_n_bs = k_sum - PSII_RATE_CONSTANTS.K_F - PSII_RATE_CONSTANTS.K_U - PSII_RATE_CONSTANTS.K_D;

    # compute PSII and PSI yeilds
    @. k_sum_1 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n_bs;
    @. k_sum_2 = PSII_RATE_CONSTANTS.K_D + PSII_RATE_CONSTANTS.K_F + PSII_RATE_CONSTANTS.K_U + psa.k_n_bs + PSII_RATE_CONSTANTS.K_P;
    k_sum_3 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_P;
    k_sum_4 = PSI_RATE_CONSTANTS.K_D + PSI_RATE_CONSTANTS.K_F + PSI_RATE_CONSTANTS.K_X;
    @. ϕ_U2_a = (q2 * PSII_RATE_CONSTANTS.K_U / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_U / k_sum_1);
    @. psa.ϕ_f2_bs = (q2 * PSII_RATE_CONSTANTS.K_F / k_sum_2 + (1 - q2) * PSII_RATE_CONSTANTS.K_F / k_sum_1) / (1 - ϕ_U2_a);
    @. psa.ϕ_f1_bs = PSI_RATE_CONSTANTS.K_F / k_sum_3 * q1 + PSI_RATE_CONSTANTS.K_F / k_sum_4 * (1 - q1);
    @. psa.ϕ_p1_bs = ϕ_P1_a;
    @. psa.ϕ_p2_bs = ϕ_P2_a;
    @. psa.ϕ_d1_bs = PSI_RATE_CONSTANTS.K_D / PSI_RATE_CONSTANTS.K_F * psa.ϕ_f1_bs;
    @. psa.ϕ_d2_bs = PSII_RATE_CONSTANTS.K_D / PSII_RATE_CONSTANTS.K_F * psa.ϕ_f2_bs;
    @. psa.ϕ_n1_bs = 1 - psa.ϕ_f1_bs - psa.ϕ_d1_bs - psa.ϕ_p1_bs;
    @. psa.ϕ_n2_bs = 1 - psa.ϕ_f2_bs - psa.ϕ_d2_bs - psa.ϕ_p2_bs;

    #
    #
    # now for the leaf level weighted value
    #
    #
    @. psa.ϕ_d = (psa.ϕ_d1_m  * (1 - psa.f_psii_m ) + psa.ϕ_d2_m  * psa.f_psii_m ) * (1 - pst.f_b6f_bs) +
                 (psa.ϕ_d1_bs * (1 - psa.f_psii_bs) + psa.ϕ_d2_bs * psa.f_psii_bs) * pst.f_b6f_bs;
    @. psa.ϕ_f = (psa.ϕ_f1_m  * (1 - psa.f_psii_m ) + psa.ϕ_f2_m  * psa.f_psii_m ) * (1 - pst.f_b6f_bs) +
                 (psa.ϕ_f1_bs * (1 - psa.f_psii_bs) + psa.ϕ_f2_bs * psa.f_psii_bs) * pst.f_b6f_bs;
    @. psa.ϕ_n = (psa.ϕ_n1_m  * (1 - psa.f_psii_m ) + psa.ϕ_n2_m  * psa.f_psii_m ) * (1 - pst.f_b6f_bs) +
                 (psa.ϕ_n1_bs * (1 - psa.f_psii_bs) + psa.ϕ_n2_bs * psa.f_psii_bs) * pst.f_b6f_bs;
    @. psa.ϕ_p = (psa.ϕ_p1_m  * (1 - psa.f_psii_m ) + psa.ϕ_p2_m  * psa.f_psii_m ) * (1 - pst.f_b6f_bs) +
                 (psa.ϕ_p1_bs * (1 - psa.f_psii_bs) + psa.ϕ_p2_bs * psa.f_psii_bs) * pst.f_b6f_bs;

    if any(isnan.(psa.ϕ_f)) || any(isnan.(psa.ϕ_p))
        println();
        @show psa.ϕ_f psa.ϕ_p;
        @show psa.a_g psa.η psa.e2c ppar psa.f_psii;
        error("NaN detected in photosystem coefficients calculation");
    end;

    return nothing
);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::Union{C3Trait{FT}, C4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::KNFluorescenceModel{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = methods;

    # calculate photochemical yield
    psa.ϕ_p .= psa.a_g ./ (psa.e2c .* psa.f_psii .* ppar);

    # unpack the cache
    x  = psa.cache1;
    xᵅ = psa.cache2;

    # calculate rate constants
    @. x  = max(0, 1 - psa.ϕ_p / psa.ϕ_psii_max);
    @. xᵅ = x ^ flm.K_A;
    @. psa.k_n = flm.K_0 .* (1 + flm.K_B) .* xᵅ ./ (flm.K_B .+ xᵅ);
    @. psa.k_p = max.(0, psa.ϕ_p .* (PS_RATE_CONSTANTS.K_F .+ psa.k_d .+ psa.k_n .+ pss.k_npq_sus) ./ (1 .- psa.ϕ_p) );

    # TODO: whether to consider sustained K_N in the calculations of f_o and f_m
    # rc._f_o  = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus);
    # rc._f_o′ = K_F / (K_F + K_PSII + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);
    # rc._f_m  = K_F / (K_F + rc._k_d + rc.k_npq_sus);
    # rc._f_m′ = K_F / (K_F + rc._k_d + rc.k_npq_sus + rc._k_npq_rev);

    # calculate fluorescence quantum yield
    psa.f_o = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;
    @. psa.ϕ_d  = psa.ϕ_f .* psa.k_d ./ PS_RATE_CONSTANTS.K_F;
    @. psa.ϕ_n  = 1 .- psa.ϕ_f .- psa.ϕ_d .- psa.ϕ_p;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::Union{C3Trait{FT}, C4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::QLFluorescenceModel{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = methods;

    # calculate photochemical yield
    @. psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = psa.cache1;
    @. q_l = exp(-flm.K_B * ppar);
    @. psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    @. psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;
    @. psa.ϕ_d  = psa.ϕ_f * psa.k_d / PS_RATE_CONSTANTS.K_F;
    @. psa.ϕ_n  = 1 - psa.ϕ_f - psa.ϕ_d - psa.ϕ_p;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);

photosystem_coefficients!(
            methods::PhotosynthesisMethods{FT},
            pst::Union{C3Trait{FT}, C4Trait{FT}},
            pss::Union{C3State{FT}, C4State{FT}},
            psa::LeafPhotosystemAuxil{FT},
            flm::QLFluorescenceModelHan{FT},
            ppar::Vector{FT};
            β::FT = FT(1)) where {FT} = (
    if ppar[1] == 0 && ppar[end] == 0
        psa.ϕ_f .= 0;
        psa.ϕ_p .= 0;

        return nothing
    end;

    (; PS_RATE_CONSTANTS) = methods;

    # calculate photochemical yield
    @. psa.ϕ_p = psa.a_g / (psa.e2c * psa.f_psii * ppar);

    # calculate the qL
    q_l = psa.cache1;
    @. q_l = flm.K_A * exp(-flm.K_B * ppar);
    @. psa.k_p = PS_RATE_CONSTANTS.K_P * q_l;
    @. psa.k_n = (psa.k_p - psa.ϕ_p * (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_p)) / psa.ϕ_p;

    # calculate fluorescence quantum yield
    psa.f_o  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + PS_RATE_CONSTANTS.K_P);
    psa.f_m  = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d);
    @. psa.f_o′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus + PS_RATE_CONSTANTS.K_P);
    @. psa.f_m′ = PS_RATE_CONSTANTS.K_F / (PS_RATE_CONSTANTS.K_F + psa.k_d + psa.k_n + pss.k_npq_sus);
    @. psa.ϕ_f  = psa.f_m′ * (1 - psa.ϕ_p);
    @. psa.ϕ_f1 = psa.ϕ_f;
    @. psa.ϕ_f2 = psa.ϕ_f;
    @. psa.ϕ_d  = psa.ϕ_f * psa.k_d / PS_RATE_CONSTANTS.K_F;
    @. psa.ϕ_n  = 1 - psa.ϕ_f - psa.ϕ_d - psa.ϕ_p;

    # TODO: if K_N is used above, do we need to recalculate _npq
    # rc._npq = (rc._k_npq_rev + rc.k_npq_sus) / (K_F + rc._k_d + rc.k_npq_sus);

    # calculate quenching rates
    @. psa.q_e = 1 - (psa.f_m - psa.f_o′) / (psa.f_m′ - psa.f_o);
    @. psa.q_p = 1 - (psa.ϕ_f - psa.f_o′) / (psa.f_m - psa.f_o′);
    @. psa.npq = (psa.k_n + pss.k_npq_sus) / (PS_RATE_CONSTANTS.K_F + psa.k_d);

    return nothing
);
