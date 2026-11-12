rng = np.random.default_rng(48)
balance_rows = []
m = min(TD.shape[1], N.shape[1])

if FOLLOWUP_EDGE_PAIRS:
    for repeat in tqdm(
        range(N_BALANCED_REPEATS),
        desc="Balanced-size sensitivity: discovery-rewired pairs",
    ):
        index = rng.choice(TD.shape[1], m, replace=False)
        small = TD.to_numpy()[:, index]
        correlations = corr(small)
        delta = correlations - RN
        pvalues = bw(correlations[IU], RN[IU], m, N.shape[1])
        qvalues = bh(pvalues)

        bl, bhigh, bv = edge_intervals(
            small,
            N.to_numpy(),
            N_BALANCED_BOOT,
            1000 + repeat,
            pair_indices=selected_iu,
        )

        for selected_k, (a, b) in enumerate(FOLLOWUP_EDGE_PAIRS):
            k = edge_key[frozenset((a, b))]
            i, j = genes.index(a), genes.index(b)

            balance_rows.append({
                "Repeat": repeat,
                "Gene1": a,
                "Gene2": b,
                "Delta": delta[i, j],
                "Same_direction": (
                    np.sign(delta[i, j]) == np.sign(D[i, j])
                ),
                "BW_q_lt_05": qvalues[k] < .05,
                "BW_and_CI_supported": bool(
                    qvalues[k] < .05
                    and bv[selected_k] >= .9 * N_BALANCED_BOOT
                    and (bl[selected_k] > 0 or bhigh[selected_k] < 0)
                ),
                "Target_policy": EDGE_FOLLOWUP_POLICY,
                "BH_family": "All original retained gene pairs",
                "Interpretation": "Same-cohort post-selection sensitivity",
            })

BAL = pd.DataFrame(balance_rows, columns=[
    "Repeat", "Gene1", "Gene2", "Delta", "Same_direction",
    "BW_q_lt_05", "BW_and_CI_supported", "Target_policy",
    "BH_family", "Interpretation",
])
tab(BAL, "S25_balanced_edge_iterations.csv")

if not BAL.empty:
    bs = BAL.groupby(["Gene1", "Gene2"]).agg(
        Median_delta=("Delta", "median"),
        Direction_frequency=("Same_direction", "mean"),
        BW_q_frequency=("BW_q_lt_05", "mean"),
        BW_and_CI_frequency=("BW_and_CI_supported", "mean"),
    ).reset_index()
else:
    bs = pd.DataFrame(columns=[
        "Gene1", "Gene2", "Median_delta", "Direction_frequency",
        "BW_q_frequency", "BW_and_CI_frequency",
    ])

bs["Target_policy"] = EDGE_FOLLOWUP_POLICY
bs["BH_family"] = "All original retained gene pairs"
bs["Interpretation"] = "Same-cohort post-selection sensitivity"
tab(bs, "S25_balanced_edge_summary.csv")

rng = np.random.default_rng(49)
refrows = []

for repeat in tqdm(range(N_REFERENCE_BOOT), desc="Reference stability"):
    index = rng.integers(0, N.shape[1], N.shape[1])
    nn = N.iloc[:, index].copy()
    nn.columns = [f"b{i}" for i in range(nn.shape[1])]

    try:
        ref = fit_ref(T, nn)
        scores = make_scores(T, ref)
        row = {"Repeat": repeat, "Success": True}

        for feature in [
            "Coordination Deviation",
            "Dysregulation Burden",
            "Activation Burden",
        ]:
            row[feature + "_rho"] = sp(
                scores[feature], SCORES[feature]
            )[0]

        refrows.append(row)

    except Exception as error:
        refrows.append({
            "Repeat": repeat,
            "Success": False,
            "Error": str(error),
        })

REFBOOT = pd.DataFrame(refrows)
tab(REFBOOT, "S10_reference_bootstrap.csv")

ZN = T.sub(REF["nm"], axis=0).div(REF["ns"], axis=0).to_numpy()
threshold_rows = []

for zthreshold in [.5, 1., 1.5]:
    for rthreshold in [.1, .2, .3]:
        for minimum_pairs in [3, 5, 10]:
            default, weighted, continuous, count, k = coordination(
                ZN, REF["R"], zthreshold, rthreshold, minimum_pairs
            )
            values = pd.Series(default, index=tids)
            number_imputed = int(values.isna().sum())
            median = values.median()

            if pd.isna(median):
                threshold_rows.append({
                    "z_threshold": zthreshold,
                    "rho_threshold": rthreshold,
                    "minimum_pairs": minimum_pairs,
                    "N_imputed": number_imputed,
                    "N_reference_pairs": k,
                    "Note": "Entire score undefined; not imputed",
                })
                continue

            values = values.fillna(median)
            row, _ = association(
                "Threshold variant", values, name="Threshold sensitivity"
            )
            row.update({
                "z_threshold": zthreshold,
                "rho_threshold": rthreshold,
                "minimum_pairs": minimum_pairs,
                "N_imputed": number_imputed,
                "Fraction_imputed": number_imputed / len(tids),
                "N_reference_pairs": k,
                "rho_vs_default": sp(
                    values, SCORES["Coordination Deviation"]
                )[0],
            })
            threshold_rows.append(row)

TH = pd.DataFrame(threshold_rows)
if "p" in TH:
    TH["q_threshold_family"] = bh(TH.p)
tab(TH, "S26_coordination_thresholds.csv")

rng = np.random.default_rng(50)
sg = 12
sigma = .6 ** np.abs(
    np.subtract.outer(np.arange(sg), np.arange(sg))
)
flip = np.diag(np.r_[
    np.ones(sg // 2), -np.ones(sg - sg // 2)
])

cov_scenarios = {
    "Unchanged": sigma,
    "Weakened": .3 * sigma + .7 * np.eye(sg),
    "Block_sign_changed": flip @ sigma @ flip,
}

mean_scenarios = {
    "No_shift": np.zeros(sg),
    "Concordant_shift": np.ones(sg),
    "Opposing_shift": np.r_[
        np.ones(sg // 2), -np.ones(sg - sg // 2)
    ],
}

simrows = []

for repeat in range(N_SIM_REPEATS):
    reference = rng.multivariate_normal(np.zeros(sg), sigma, 50).T
    mean = reference.mean(axis=1, keepdims=True)
    sd = reference.std(axis=1, ddof=1, keepdims=True)
    reference_corr = corr(reference)

    for covariance_name, covariance in cov_scenarios.items():
        for mean_name, mean_vector in mean_scenarios.items():
            tumour = rng.multivariate_normal(
                mean_vector, covariance, 300
            ).T
            standardized = (tumour - mean) / sd
            default, weighted, continuous, count, k = coordination(
                standardized, reference_corr
            )

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                median_default = np.nanmedian(default)

            simrows.append({
                "Repeat": repeat,
                "Covariance": covariance_name,
                "Mean_shift": mean_name,
                "Median_default_evaluable": median_default,
                "Mean_continuous": np.nanmean(continuous),
                "Fraction_default_undefined": np.mean(
                    ~np.isfinite(default)
                ),
                "Observed_mean_abs_delta_rho": np.abs(
                    (corr(tumour) - reference_corr)[
                        np.triu_indices(sg, 1)
                    ]
                ).mean(),
            })

SIM = pd.DataFrame(simrows)
tab(SIM, "S27_construct_simulation_iterations.csv")

tab(
    SIM.groupby(["Covariance", "Mean_shift"]).agg(
        Median_default=("Median_default_evaluable", "median"),
        Mean_continuous=("Mean_continuous", "mean"),
        Mean_undefined=("Fraction_default_undefined", "mean"),
        Mean_abs_delta_rho=("Observed_mean_abs_delta_rho", "mean"),
    ).reset_index(),
    "S27_construct_simulation_summary.csv",
)

nsd = N.std(axis=1)
nmed = N.median(axis=1)
mad = N.sub(nmed, axis=0).abs().median(axis=1) * 1.4826
iqr = (N.quantile(.75, axis=1) - N.quantile(.25, axis=1)) / 1.349

robust_sd = mad.copy()
fallback = robust_sd <= EPS
robust_sd[fallback] = iqr[fallback]
unusable = robust_sd <= EPS

reference_gene = pd.DataFrame({
    "Gene": genes,
    "Adjacent_SD": nsd.values,
    "Scaled_MAD": mad.values,
    "IQR_scale": iqr.values,
    "Robust_scale_used": robust_sd.values,
    "MAD_to_IQR_fallback": fallback.values,
    "Robust_unusable": unusable.values,
    "Mean_abs_z_contribution": np.abs(ZN).mean(axis=1),
    "Mean_positive_z_contribution": np.maximum(ZN, 0).mean(axis=1),
})
reference_gene["Abs_contribution_share"] = (
    reference_gene.Mean_abs_z_contribution
    / reference_gene.Mean_abs_z_contribution.sum()
)
tab(reference_gene, "S28_reference_variance_gene_contributions.csv")

robust_rows = []
usable = robust_sd.index[~unusable]

if len(usable) >= 6:
    rz = T.loc[usable].sub(nmed.loc[usable], axis=0).div(
        robust_sd.loc[usable], axis=0
    )

    for feature, values in [
        ("Robust Dysregulation Burden", rz.abs().mean(axis=0)),
        ("Robust Activation Burden", rz.clip(lower=0).mean(axis=0)),
    ]:
        row, _ = association(feature, values, name="Robust scaling")
        standard = (
            "Dysregulation Burden"
            if "Dysregulation" in feature else "Activation Burden"
        )
        row["rho_vs_standard"] = sp(values, SCORES[standard])[0]
        row["N_genes"] = len(usable)
        robust_rows.append(row)

floor = float(nsd.quantile(.10))
fz = T.sub(REF["nm"], axis=0).div(nsd.clip(lower=floor), axis=0)

for feature, values in [
    ("SDfloor Dysregulation Burden", fz.abs().mean(axis=0)),
    ("SDfloor Activation Burden", fz.clip(lower=0).mean(axis=0)),
]:
    row, _ = association(feature, values, name="SD floor")
    row["SD_floor"] = floor
    robust_rows.append(row)

rob = pd.DataFrame(robust_rows)
if len(rob):
    rob["q_scaling_sensitivity_family"] = bh(rob.p)
tab(rob, "S28_robust_scaling_OS.csv")

loo = []
znframe = pd.DataFrame(ZN, index=genes, columns=tids)

for omitted in genes:
    reduced = znframe.drop(index=omitted)

    for feature, score in [
        ("Dysregulation Burden", reduced.abs().mean(axis=0)),
        ("Activation Burden", reduced.clip(lower=0).mean(axis=0)),
    ]:
        row, _ = association(feature, score, name="Leave_one_gene_out")
        row["Omitted_gene"] = omitted
        row["rho_vs_full"] = sp(score, SCORES[feature])[0]
        loo.append(row)

LOO = pd.DataFrame(loo)
LOO["q_LOO_sensitivity_family"] = bh(LOO.p)
tab(LOO, "S29_leave_one_gene_out.csv")