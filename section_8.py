def repeated_cv(ids, targets, clinical_base=None, label="Molecular"):
    yy = Y.loc[ids]
    records = []

    if (
        not targets
        or len(ids) < 60
        or yy.E.sum() < 15
        or yy.E.value_counts().min() < N_FOLDS
    ):
        NOTICES.append(label + " CV skipped: insufficient data or no targets.")
        return pd.DataFrame()

    for seed in CV_SEEDS:
        splitter = StratifiedKFold(
            N_FOLDS, shuffle=True, random_state=seed
        )
        numerator = {feature: 0. for feature in targets}
        denominator = {feature: 0 for feature in targets}
        successes = {feature: 0 for feature in targets}
        baseline_numerator = 0.
        baseline_denominator = 0
        baseline_successes = 0

        for fold, (itr, ite) in enumerate(splitter.split(yy, yy.E), 1):
            train = yy.index[itr]
            test = yy.index[ite]
            held_patients = {pid(sample) for sample in test}
            ns = [
                sample for sample in N
                if pid(sample) not in held_patients
            ]

            try:
                ref = fit_ref(T[train], N[ns])
                trscore = make_scores(T[train], ref)
                tescore = make_scores(T[test], ref)
                raw_trscore = make_scores(T[train], ref, False)
            except Exception as error:
                for feature in targets:
                    CV_LOG.append({
                        "Label": label,
                        "Seed": seed,
                        "Fold": fold,
                        "Feature": feature,
                        "N_reference_adjacent": len(ns),
                        "N_train": len(train),
                        "N_test": len(test),
                        "Success": False,
                        "Error": str(error),
                        "Selection_policy": MODEL_SELECTION_POLICY,
                        "Interpretation": CV_INTERPRETATION,
                    })
                continue

            if clinical_base is not None:
                clinical_train = clinical_base.loc[train].copy()
                clinical_test = clinical_base.loc[test].copy()
                variable_cols = clinical_train.columns[
                    clinical_train.nunique() > 1
                ]
                clinical_train = clinical_train[variable_cols]
                clinical_test = clinical_test[variable_cols]

                for col in clinical_train:
                    if clinical_train[col].nunique() > 2:
                        mean = clinical_train[col].mean()
                        sd = clinical_train[col].std()
                        clinical_train[col] = (
                            clinical_train[col] - mean
                        ) / sd
                        clinical_test[col] = (
                            clinical_test[col] - mean
                        ) / sd

                if len(variable_cols):
                    baseline_data = yy.loc[train].join(clinical_train)
                    baseline_model, _ = fit_cox(
                        baseline_data,
                        f"{label}_{seed}_{fold}_baseline",
                    )

                    if baseline_model is not None:
                        baseline_risk = (
                            baseline_model.predict_log_partial_hazard(
                                clinical_test
                            ).to_numpy()
                        )
                        baseline_c = cindex(
                            yy.loc[test, "T"],
                            yy.loc[test, "E"],
                            baseline_risk,
                        )
                    else:
                        baseline_c = np.nan
                else:
                    baseline_c = cindex(
                        yy.loc[test, "T"],
                        yy.loc[test, "E"],
                        np.zeros(len(test)),
                    )

                weight = comparable(yy.loc[test, "T"], yy.loc[test, "E"])
                if np.isfinite(baseline_c) and weight:
                    baseline_numerator += weight * baseline_c
                    baseline_denominator += weight
                    baseline_successes += 1

            for feature in targets:
                tx = (
                    T.loc[feature, train]
                    if feature in genes else trscore[feature]
                )
                vx = (
                    T.loc[feature, test]
                    if feature in genes else tescore[feature]
                )
                mean, sd = tx.mean(), tx.std()

                audit = {
                    "Label": label,
                    "Seed": seed,
                    "Fold": fold,
                    "Feature": feature,
                    "N_reference_adjacent": len(ns),
                    "N_train": len(train),
                    "N_test": len(test),
                    "N_train_imputed": (
                        int(raw_trscore[feature].isna().sum())
                        if feature in [
                            "Coordination Deviation", "Coordination weighted"
                        ] else 0
                    ),
                    "Success": False,
                    "Selection_policy": MODEL_SELECTION_POLICY,
                    "Interpretation": CV_INTERPRETATION,
                }

                if (
                    not np.isfinite(sd)
                    or sd <= EPS
                    or not np.isfinite(tx.to_numpy()).all()
                    or not np.isfinite(vx.to_numpy()).all()
                ):
                    CV_LOG.append(audit)
                    continue

                train_data = yy.loc[train].copy()
                test_data = pd.DataFrame(index=test)

                if clinical_base is not None:
                    train_data = train_data.join(clinical_train)
                    test_data = test_data.join(clinical_test)

                train_data["X"] = (tx - mean) / sd
                test_data["X"] = (vx - mean) / sd

                model, penalty = fit_cox(
                    train_data, f"{label}_{seed}_{fold}_{feature}"
                )

                if model is None:
                    CV_LOG.append(audit)
                    continue

                risk = model.predict_log_partial_hazard(test_data).to_numpy()
                value = cindex(yy.loc[test, "T"], yy.loc[test, "E"], risk)
                weight = comparable(yy.loc[test, "T"], yy.loc[test, "E"])

                if np.isfinite(value) and weight:
                    numerator[feature] += weight * value
                    denominator[feature] += weight
                    successes[feature] += 1
                    audit.update({
                        "Success": True,
                        "Cindex": value,
                        "Penalty": penalty,
                    })

                CV_LOG.append(audit)

        baseline_value = (
            baseline_numerator / baseline_denominator
            if baseline_successes == N_FOLDS and baseline_denominator
            else np.nan
        )

        for feature in targets:
            value = (
                numerator[feature] / denominator[feature]
                if successes[feature] == N_FOLDS and denominator[feature]
                else np.nan
            )
            records.append({
                "Label": label,
                "Seed": seed,
                "Feature": feature,
                "Cindex": value,
                "Clinical_baseline_Cindex": baseline_value,
                "Difference_vs_baseline": value - baseline_value,
                "Complete_folds": successes[feature],
                "Selection_policy": MODEL_SELECTION_POLICY,
                "Interpretation": CV_INTERPRETATION,
            })

    return pd.DataFrame(records)


def cv_summary(df):
    if df.empty:
        return pd.DataFrame()

    rows = []

    for feature, subset in df.groupby("Feature"):
        complete = subset.Cindex.notna().sum() == len(CV_SEEDS)
        complete_baseline = (
            subset.Clinical_baseline_Cindex.notna().sum()
            == len(CV_SEEDS)
        )
        complete_delta = (
            subset.Difference_vs_baseline.notna().sum()
            == len(CV_SEEDS)
        )

        rows.append({
            "Feature": feature,
            "CV_mean": subset.Cindex.mean() if complete else np.nan,
            "CV_split_SD": (
                subset.Cindex.std(ddof=1) if complete else np.nan
            ),
            "Complete_repeats": int(subset.Cindex.notna().sum()),
            "Baseline_mean": (
                subset.Clinical_baseline_Cindex.mean()
                if complete_baseline else np.nan
            ),
            "Descriptive_delta_C": (
                subset.Difference_vs_baseline.mean()
                if complete_delta else np.nan
            ),
            "Selection_policy": MODEL_SELECTION_POLICY,
            "Interpretation": CV_INTERPRETATION,
        })

    return pd.DataFrame(rows)


# All outcome-independent score definitions, plus outcome-selected genes.
MOLECULAR_CV_FEATURES = list(dict.fromkeys(
    MAIN_SCORES + ALT_SCORES + TOP_GENES
))

CV = repeated_cv(
    tids,
    MOLECULAR_CV_FEATURES,
    label="All_scores_and_selected_genes",
)
tab(CV, "CV_molecular_repeats.csv")

CVSUM = cv_summary(CV)
tab(CVSUM, "S20_CV_molecular_summary.csv")

comparison = pd.concat([
    SCORE_OS.assign(Type="Fixed module/component definition"),
    GENE_OS.assign(Type="Retained individual gene"),
], ignore_index=True)

comparison["Selected_for_adjusted_Cox_and_CV"] = (
    comparison.Feature.isin(MODEL_FEATURES)
)
comparison["Followup_policy"] = MODEL_SELECTION_POLICY

if not CVSUM.empty:
    comparison = comparison.merge(
        CVSUM[
            ["Feature", "CV_mean", "CV_split_SD", "Complete_repeats",
             "Baseline_mean", "Descriptive_delta_C", "Interpretation"]
        ],
        on="Feature",
        how="left",
    )

tab(comparison, "Table03_OS_comparison.csv")

if (
    MODEL_FEATURES
    and len(BASE_CLIN) >= 60
    and Y.loc[BASE_CLIN.index, "E"].sum() >= 15
):
    CVC = repeated_cv(
        BASE_CLIN.index.tolist(),
        MODEL_FEATURES,
        BASE_CLIN,
        "Clinical_plus_selected_feature",
    )
    tab(CVC, "S18_clinical_plus_score_CV_repeats.csv")
    tab(cv_summary(CVC), "S18_clinical_plus_score_CV_summary.csv")