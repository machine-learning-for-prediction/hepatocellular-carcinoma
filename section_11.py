coverage_rows = []

for cohort_name, ids in [
    ("All selected tumours", tids),
    ("Patient-disjoint tumours", dids),
    ("Clinical baseline complete cases", BASE_CLIN.index.tolist()),
]:
    values = purity.reindex(ids)
    nonmissing = values.dropna()

    coverage_rows.append({
        "Cohort": cohort_name,
        "N_total": len(ids),
        "N_ABSOLUTE_available": int(values.notna().sum()),
        "N_ABSOLUTE_missing": int(values.isna().sum()),
        "Fraction_ABSOLUTE_available": (
            float(values.notna().mean()) if len(values) else np.nan
        ),
        "Median_ABSOLUTE_purity": (
            float(nonmissing.median()) if len(nonmissing) else np.nan
        ),
        "Q1_ABSOLUTE_purity": (
            float(nonmissing.quantile(.25)) if len(nonmissing) else np.nan
        ),
        "Q3_ABSOLUTE_purity": (
            float(nonmissing.quantile(.75)) if len(nonmissing) else np.nan
        ),
        "N_purity_ge_070": int((values >= .70).sum()),
    })

tab(pd.DataFrame(coverage_rows), "S39_ABSOLUTE_purity_coverage.csv")

combined_covariates = pd.DataFrame(index=tids)
combined_covariates["ABSOLUTE Purity"] = purity

if PU_ESTIMATE_AVAILABLE:
    combined_covariates = combined_covariates.join(
        PU_ESTIMATE.reindex(tids)
    )

tab(
    combined_covariates.reset_index(names="sample"),
    "S39_tumour_purity_ESTIMATE_covariates.csv",
)

purity_association_rows = []

for feature in FOLLOWUP_FEATURES:
    values = feature_values(feature)
    for covariate in combined_covariates:
        rho, pvalue, count = sp(values, combined_covariates[covariate])
        purity_association_rows.append({
            "Feature": feature,
            "Covariate": covariate,
            "rho": rho,
            "p": pvalue,
            "N": count,
        })

purity_association_table = pd.DataFrame(purity_association_rows)
purity_association_table["q_exploratory_family"] = bh(
    purity_association_table.p
)
tab(
    purity_association_table,
    "S36_score_purity_ESTIMATE_associations.csv",
)

if PU_ESTIMATE_AVAILABLE:
    crosscheck_rows = []

    for covariate in [
        "ESTIMATE Stromal", "ESTIMATE Immune", "ESTIMATE Score"
    ]:
        rho, pvalue, count = sp(
            purity, PU_ESTIMATE.reindex(tids)[covariate]
        )
        crosscheck_rows.append({
            "Comparison": "ABSOLUTE Purity vs " + covariate,
            "rho": rho,
            "p": pvalue,
            "N": count,
        })

    crosscheck = pd.DataFrame(crosscheck_rows)
    crosscheck["q_crosscheck_family"] = bh(crosscheck.p)
    tab(crosscheck, "S39_ABSOLUTE_ESTIMATE_crosscheck.csv")

pu_high_count = int((purity.reindex(dids) >= .70).sum())

PU_STATUS.append({
    "Analysis": "High-purity edge sensitivity",
    "Status": (
        "Completed"
        if (TABLES / "S23_high_purity_edges.csv").exists()
        else "Not performed"
    ),
    "N_high_purity_disjoint_ge_070": pu_high_count,
    "N_target_pairs": len(FOLLOWUP_EDGE_PAIRS),
    "Target_policy": EDGE_FOLLOWUP_POLICY,
    "Note": (
        "Requires at least 30 disjoint tumours with purity >= 0.70 "
        "and at least one primary rewired pair; adjacent reference unchanged"
    ),
})

estimate_os_rows = []

if PU_ESTIMATE_AVAILABLE:
    for extra in [
        ["ESTIMATE Score"],
        ["ESTIMATE Immune", "ESTIMATE Stromal"],
    ]:
        design = BASE_CLIN.join(
            PU_ESTIMATE.reindex(tids)[extra]
        ).replace([np.inf, -np.inf], np.nan).dropna()

        if (
            len(design) < 60
            or Y.loc[design.index, "E"].sum() < 15
            or not MODEL_FEATURES
        ):
            PU_STATUS.append({
                "Analysis": "Clinical + " + " + ".join(extra),
                "Status": "Skipped",
                "Reason": "Insufficient complete cases/events or no targets",
                "N": len(design),
            })
            continue

        for feature in MODEL_FEATURES:
            row, _ = association(
                feature,
                feature_values(feature),
                design,
                "Clinical_plus_" + "_".join(extra),
            )
            row["Adjustment"] = (
                "Age, sex, categorical stage, " + ", ".join(extra)
            )
            estimate_os_rows.append(row)

    estimate_os_table = pd.DataFrame(estimate_os_rows)
    if len(estimate_os_table):
        estimate_os_table[
            "q_ESTIMATE_OS_postselection_exploratory"
        ] = bh(estimate_os_table.p)
    tab(estimate_os_table, "S38_ESTIMATE_adjusted_OS_sensitivity.csv")


def pu_partial_rank_corr(expression, covariates):
    expression = np.asarray(expression, float)
    covariates = np.asarray(covariates, float)

    if expression.shape[1] != covariates.shape[0]:
        raise ValueError("Expression and covariate sample counts differ.")
    if not np.isfinite(expression).all():
        raise ValueError("Non-finite residual-analysis expression.")
    if not np.isfinite(covariates).all():
        raise ValueError("Non-finite residual-analysis covariates.")

    return residual_corr(expression, covariates)


estimate_adjusted_edge_rows = []

if PU_ESTIMATE_AVAILABLE:
    differentiation_name = "Hepatocyte differentiation proxy"
    edge_covariates = PU_ESTIMATE[
        ["ESTIMATE Immune", "ESTIMATE Stromal"]
    ].join(PROXIES[[differentiation_name]], how="left")

    tumour_edge_ids = [
        sample for sample in dids
        if sample in edge_covariates.index
        and np.isfinite(edge_covariates.loc[sample].to_numpy(float)).all()
    ]
    adjacent_edge_ids = [
        sample for sample in nids
        if sample in edge_covariates.index
        and np.isfinite(edge_covariates.loc[sample].to_numpy(float)).all()
    ]

    if (
        len(tumour_edge_ids) >= 30
        and len(adjacent_edge_ids) >= 20
        and FOLLOWUP_EDGE_PAIRS
    ):
        tvalues = T.loc[genes, tumour_edge_ids].to_numpy()
        nvalues = N.loc[genes, adjacent_edge_ids].to_numpy()
        tcov = edge_covariates.loc[tumour_edge_ids].to_numpy(float)
        ncov = edge_covariates.loc[adjacent_edge_ids].to_numpy(float)

        raw_tcorr, raw_ncorr = corr(tvalues), corr(nvalues)
        adjusted_tcorr = pu_partial_rank_corr(tvalues, tcov)
        adjusted_ncorr = pu_partial_rank_corr(nvalues, ncov)
        adjusted_delta = adjusted_tcorr - adjusted_ncorr

        target_pairs = FOLLOWUP_EDGE_PAIRS
        pair_indices = [
            (genes.index(a), genes.index(b))
            for a, b in target_pairs
        ]
        rng_estimate = np.random.default_rng(260101)
        boot_estimate = np.full(
            (PU_N_ESTIMATE_BOOT, len(target_pairs)), np.nan
        )

        for repeat in tqdm(
            range(PU_N_ESTIMATE_BOOT),
            desc="ESTIMATE-adjusted bootstrap: discovery-rewired pairs",
        ):
            it = rng_estimate.integers(
                0, len(tumour_edge_ids), len(tumour_edge_ids)
            )
            inn = rng_estimate.integers(
                0, len(adjacent_edge_ids), len(adjacent_edge_ids)
            )
            delta_boot = (
                pu_partial_rank_corr(tvalues[:, it], tcov[it])
                - pu_partial_rank_corr(nvalues[:, inn], ncov[inn])
            )

            for k, (i, j) in enumerate(pair_indices):
                boot_estimate[repeat, k] = delta_boot[i, j]

        required_valid = max(
            20, int(np.ceil(.90 * PU_N_ESTIMATE_BOOT))
        )

        for k, ((a, b), (i, j)) in enumerate(
            zip(target_pairs, pair_indices)
        ):
            valid_boot = boot_estimate[:, k]
            valid_boot = valid_boot[np.isfinite(valid_boot)]
            low, high_ci = ci(valid_boot, required_valid)
            raw_delta = raw_tcorr[i, j] - raw_ncorr[i, j]
            adjusted_value = adjusted_delta[i, j]

            estimate_adjusted_edge_rows.append({
                "Gene1": a,
                "Gene2": b,
                "N_tumour_disjoint": len(tumour_edge_ids),
                "N_adjacent": len(adjacent_edge_ids),
                "Primary_TCGA_delta": D[i, j],
                "Same_subset_unadjusted_delta": raw_delta,
                "Adjusted_tumour_partial_rank_rho": adjusted_tcorr[i, j],
                "Adjusted_adjacent_partial_rank_rho": adjusted_ncorr[i, j],
                "ESTIMATE_adjusted_delta": adjusted_value,
                "Adjusted_boot_CI_low": low,
                "Adjusted_boot_CI_high": high_ci,
                "Bootstrap_valid": len(valid_boot),
                "Bootstrap_requested": PU_N_ESTIMATE_BOOT,
                "Same_direction_as_same_subset_unadjusted": (
                    np.sign(adjusted_value) == np.sign(raw_delta)
                ),
                "Adjustment": (
                    "ESTIMATE immune + ESTIMATE stromal + "
                    "hepatocyte differentiation proxy"
                ),
                "Target_policy": EDGE_FOLLOWUP_POLICY,
                "Interpretation": (
                    "Same-cohort post-selection bulk-expression sensitivity; "
                    "pointwise CI; no primary BH rewiring classification"
                ),
            })

        tab(
            pd.DataFrame(estimate_adjusted_edge_rows),
            "S37_ESTIMATE_adjusted_edge_sensitivity.csv",
        )

        PU_STATUS.append({
            "Analysis": "ESTIMATE-adjusted edge sensitivity",
            "Status": "Completed",
            "N_tumour_disjoint": len(tumour_edge_ids),
            "N_adjacent": len(adjacent_edge_ids),
            "N_pairs": len(estimate_adjusted_edge_rows),
            "Bootstrap_resamples": PU_N_ESTIMATE_BOOT,
            "Target_policy": EDGE_FOLLOWUP_POLICY,
        })
    else:
        PU_STATUS.append({
            "Analysis": "ESTIMATE-adjusted edge sensitivity",
            "Status": "Skipped",
            "N_tumour_disjoint": len(tumour_edge_ids),
            "N_adjacent": len(adjacent_edge_ids),
            "N_target_pairs": len(FOLLOWUP_EDGE_PAIRS),
            "Reason": (
                "Insufficient complete-case samples or no primary rewired pairs"
            ),
        })

pu_fig, pu_axes = plt.subplots(1, 2, figsize=(12, 4.5))
pu_axes[0].hist(
    purity.dropna(), bins=20, color="#20639B", edgecolor="white"
)
pu_axes[0].axvline(
    .70, color="black", linestyle="--",
    label="High-purity sensitivity threshold",
)
pu_axes[0].set_xlabel("Published ABSOLUTE purity")
pu_axes[0].set_ylabel("Number of selected tumours")
pu_axes[0].set_title("A: available tumour-purity estimates", loc="left")
pu_axes[0].legend(fontsize=8)

if PU_ESTIMATE_AVAILABLE:
    plot_data = pd.concat([
        purity.rename("ABSOLUTE Purity"),
        PU_ESTIMATE.reindex(tids)["ESTIMATE Score"],
    ], axis=1).dropna()

    pu_axes[1].scatter(
        plot_data["ABSOLUTE Purity"],
        plot_data["ESTIMATE Score"],
        s=18, alpha=.65, color="#20639B",
    )
    pu_axes[1].set_xlabel("Published ABSOLUTE purity")
    pu_axes[1].set_ylabel("RNA-seq ESTIMATE score")
    pu_axes[1].set_title(
        "B: DNA-based purity versus expression score", loc="left"
    )
else:
    pu_axes[1].axis("off")
    pu_axes[1].text(
        .5, .5, "ESTIMATE unavailable or disabled",
        ha="center", va="center",
    )

pu_fig.tight_layout()
figsave(pu_fig, "FigS_ABSOLUTE_and_ESTIMATE")
tab(pd.DataFrame(PU_STATUS), "S39_purity_ESTIMATE_analysis_status.csv")