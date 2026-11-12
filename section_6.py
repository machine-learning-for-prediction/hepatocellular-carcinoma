def coordination(
    Z, R, zthr=COORD_Z, rthr=COORD_RHO, minpairs=COORD_PAIRS
):
    ii, jj = np.triu_indices(Z.shape[0], 1)
    rv = R[ii, jj]
    keep_pairs = np.isfinite(rv) & (np.abs(rv) >= rthr)
    ii, jj, rv = ii[keep_pairs], jj[keep_pairs], rv[keep_pairs]

    a, b = Z[ii], Z[jj]
    eligible = (np.abs(a) >= zthr) & (np.abs(b) >= zthr)
    signs = np.sign(a * b)
    eligible &= signs != 0
    violation = eligible & (signs != np.sign(rv)[:, None])
    count = eligible.sum(axis=0)
    good = count >= minpairs

    default = np.full(Z.shape[1], np.nan)
    weighted = default.copy()
    default[good] = violation[:, good].sum(axis=0) / count[good]

    denominator = (eligible * np.abs(rv)[:, None]).sum(axis=0)
    weighted_good = good & (denominator > EPS)
    weighted[weighted_good] = (
        violation * np.abs(rv)[:, None]
    ).sum(axis=0)[weighted_good] / denominator[weighted_good]

    continuous = (
        np.mean(
            -np.tanh(a) * np.tanh(b) * np.sign(rv)[:, None],
            axis=0,
        )
        if len(ii) else np.full(Z.shape[1], np.nan)
    )

    return default, weighted, continuous, count, len(ii)


def make_scores(xx, ref, impute=True):
    gg = ref["genes"]
    values = xx.loc[gg]
    zt = values.sub(ref["tm"], axis=0).div(ref["ts"], axis=0)
    zn = values.sub(ref["nm"], axis=0).div(ref["ns"], axis=0)
    out = pd.DataFrame(index=xx.columns)

    def module_mean(gl):
        available = [gene for gene in gl if gene in gg]
        return zt.loc[available].mean(axis=0) if available else np.nan

    out["Core Synthesis"] = module_mean(MODULES["Core synthesis"])
    out["Catabolism"] = module_mean(MODULES["Catabolism"])
    out["Transport-associated"] = module_mean(MODULES["Transport-associated"])
    out["Synthesis Support"] = module_mean(SUPPORT)
    out["Homeostatic Brake"] = module_mean(BRAKE)
    out["Balance Index"] = (
        out["Synthesis Support"] - out["Homeostatic Brake"]
    )
    out["Dysregulation Burden"] = zn.abs().mean(axis=0)
    out["Activation Burden"] = zn.clip(lower=0).mean(axis=0)
    out["Signed Burden"] = zn.mean(axis=0)

    default, weighted, continuous, count, _ = coordination(
        zn.to_numpy(), ref["R"]
    )
    out["Coordination Deviation"] = default
    out["Coordination weighted"] = weighted
    out["Coordination continuous"] = continuous
    out["Informative pairs"] = count

    if impute:
        for feature, value in ref["impute"].items():
            out[feature] = out[feature].fillna(value)

    return out


def fit_ref(tt, nn):
    eligible = [
        gene for gene in tt.index
        if gene in nn.index
        and tt.loc[gene].std() > EPS
        and nn.loc[gene].std() > EPS
        and (tt.loc[gene] > 0).mean() >= .10
    ]

    if len(eligible) < 6:
        raise RuntimeError("Too few score-reference genes.")

    a, b = tt.loc[eligible], nn.loc[eligible]
    ref = {
        "genes": eligible,
        "tm": a.mean(axis=1),
        "ts": a.std(axis=1),
        "nm": b.mean(axis=1),
        "ns": b.std(axis=1),
        "R": corr(b),
    }

    rawscore = make_scores(tt, ref, impute=False)
    ref["impute"] = {}

    for feature in ["Coordination Deviation", "Coordination weighted"]:
        values = rawscore[feature].dropna()
        if not len(values):
            raise RuntimeError("No evaluable " + feature)
        ref["impute"][feature] = float(values.median())

    return ref


REF = fit_ref(T, N)
RAW_SCORES = make_scores(T, REF, False)
SCORES = make_scores(T, REF)

tab(SCORES.reset_index(names="sample"), "patient_scores.csv")

imputation = pd.DataFrame([
    {
        "Feature": feature,
        "N_imputed": int(RAW_SCORES[feature].isna().sum()),
        "Fraction_imputed": float(RAW_SCORES[feature].isna().mean()),
        "Median": REF["impute"][feature],
    }
    for feature in ["Coordination Deviation", "Coordination weighted"]
])
tab(imputation, "S14_imputation.csv")

dependence = []
dependence_features = MAIN_SCORES + ["Informative pairs"]

for i, a in enumerate(dependence_features):
    for b in dependence_features[i + 1:]:
        rho, pvalue, count = sp(SCORES[a], SCORES[b])
        dependence.append({
            "Feature1": a,
            "Feature2": b,
            "rho": rho,
            "p": pvalue,
            "N": count,
        })

dep = pd.DataFrame(dependence)
dep["q"] = bh(dep.p)
tab(dep, "S09_score_dependence.csv")


def fit_cox(data, name, strata=None):
    for penalty in [0.] + RIDGES:
        try:
            with warnings.catch_warnings(record=True) as recorded:
                warnings.simplefilter("always")
                model = CoxPHFitter(penalizer=penalty)
                model.fit(data, "T", "E", strata=strata)

            convergence_warnings = [
                str(item.message) for item in recorded
                if issubclass(item.category, ConvergenceWarning)
            ]
            if convergence_warnings:
                raise RuntimeError("; ".join(convergence_warnings))
            if not np.isfinite(model.params_).all():
                raise RuntimeError("Nonfinite coefficient.")

            LOG.append({
                "Model": name,
                "Penalty": penalty,
                "Success": True,
                "Warnings": [str(item.message) for item in recorded],
            })
            return model, penalty

        except Exception as error:
            LOG.append({
                "Model": name,
                "Penalty": penalty,
                "Success": False,
                "Error": str(error),
            })

    return None, np.nan


def association(
    feature, values, base=None, name="Univariable", y=None
):
    outcomes = Y if y is None else y
    data = outcomes[["T", "E"]].copy()

    if base is not None:
        data = data.join(base)

    data["X"] = pd.Series(values).reindex(data.index)
    data = data.replace([np.inf, -np.inf], np.nan).dropna()

    selected_adjusted = (
        feature in MODEL_FEATURES
        and (
            name.startswith("Age_sex_stage")
            or name.startswith("Clinical_plus_")
        )
    )

    row = {
        "Feature": feature,
        "Model": name,
        "N": len(data),
        "Events": int(data.E.sum()),
        "HR": np.nan,
        "CI_low": np.nan,
        "CI_high": np.nan,
        "p": np.nan,
        "PH_X_p": np.nan,
        "C_apparent": np.nan,
        "Penalty": np.nan,
        "Feature_selection_policy": (
            MODEL_SELECTION_POLICY
            if selected_adjusted
            else "Not restricted by the adjusted-model shortlist"
        ),
        "Inference_scope": (
            "Exploratory same-cohort post-selection inference"
            if selected_adjusted
            else "See model-specific test family"
        ),
    }

    if (
        len(data) < 50
        or data.E.sum() < 10
        or not np.isfinite(data.X.std())
        or data.X.std() <= EPS
    ):
        row["Note"] = "Insufficient evaluable data"
        return row, None

    data["X"] = (data.X - data.X.mean()) / data.X.std()

    nuisance = [
        col for col in data
        if col not in ["T", "E", "X"]
    ]
    constant = [col for col in nuisance if data[col].nunique() <= 1]
    data = data.drop(columns=constant)

    for col in data:
        if col not in ["T", "E", "X"] and data[col].nunique() > 2:
            data[col] = (data[col] - data[col].mean()) / data[col].std()

    model, penalty = fit_cox(data, name + "_" + feature)

    if model is None:
        row["Note"] = "Cox fitting failed"
        return row, None

    result = model.summary.loc["X"]
    risk = model.predict_log_partial_hazard(data).to_numpy()

    row.update({
        "HR": result["exp(coef)"],
        "CI_low": result["exp(coef) lower 95%"],
        "CI_high": result["exp(coef) upper 95%"],
        "p": result["p"],
        "C_apparent": cindex(data["T"], data.E, risk),
        "Penalty": penalty,
        "Note": (
            "Unpenalised Wald inference"
            if penalty == 0 else "Ridge fallback; inference exploratory"
        ),
    })

    try:
        ph = proportional_hazard_test(
            model, data, time_transform="rank"
        )
        row["PH_X_p"] = float(ph.summary.loc["X", "p"])
        tab(
            ph.summary.reset_index(names="Term"),
            "PH_"
            + re.sub(r"[^A-Za-z0-9]+", "_", name + "_" + feature)
            + ".csv",
        )
    except Exception as error:
        row["PH_note"] = str(error)

    return row, (model, data, risk)


GENE_OS = pd.DataFrame([
    association(gene, T.loc[gene])[0] for gene in genes
])
GENE_OS["q_gene_family"] = bh(GENE_OS.p)
tab(GENE_OS, "S08_gene_OS_complete.csv")

SCORE_OS = pd.DataFrame([
    association(feature, SCORES[feature])[0]
    for feature in MAIN_SCORES
])
SCORE_OS["q_module_component_family"] = bh(SCORE_OS.p)
SCORE_OS["PH_q_within_family"] = bh(SCORE_OS.PH_X_p)
tab(SCORE_OS, "S15_module_component_OS_complete.csv")


def select_leading_features(results, n):
    eligible = results.loc[
        results["Penalty"].eq(0)
        & np.isfinite(results["p"])
        & np.isfinite(results["HR"])
    ].copy()

    qcol = (
        "q_gene_family"
        if "q_gene_family" in eligible.columns
        else "q_module_component_family"
    )
    eligible["selection_q"] = eligible[qcol]
    eligible = eligible.sort_values(
        ["selection_q", "p", "Feature"], kind="stable"
    ).head(n)
    eligible["Selection_rank"] = np.arange(1, len(eligible) + 1)
    return eligible


leading_scores = select_leading_features(SCORE_OS, TOP_N_SCORES)
leading_genes = select_leading_features(GENE_OS, TOP_N_GENES)

TOP_SCORES = leading_scores.Feature.tolist()
TOP_GENES = leading_genes.Feature.tolist()
MODEL_FEATURES = TOP_SCORES + TOP_GENES + ALT_SCORES

selected_features = pd.concat([
    leading_scores.assign(Selection_family="Main score"),
    leading_genes.assign(Selection_family="Gene"),
], ignore_index=True)
selected_features["Selection_policy"] = MODEL_SELECTION_POLICY
selected_features["Selection_requires_q_lt_05"] = False

tab(selected_features, "Table06_outcome_selected_model_features.csv")

tab(pd.DataFrame([
    {
        "Feature": feature,
        "Type": "Gene" if feature in genes else "Module/component",
        "Initial_univariable_OS_eligible": True,
        "Selected_for_adjusted_Cox_and_CV": feature in MODEL_FEATURES,
        "Selection_policy": MODEL_SELECTION_POLICY,
    }
    for feature in FOLLOWUP_FEATURES
]), "Table07_feature_followup_policy.csv")

print("\nSelected leading scores:", TOP_SCORES)
print("Selected leading genes:", TOP_GENES)

if len(TOP_SCORES) < TOP_N_SCORES or len(TOP_GENES) < TOP_N_GENES:
    NOTICES.append(
        "available eligible features retained."
    )

NOTICES.append(CV_INTERPRETATION)

lioness = []

for i in range(T.shape[1]):
    reduced = corr(np.delete(T.to_numpy(), i, axis=1))
    single = T.shape[1] * RA - (T.shape[1] - 1) * reduced
    lioness.append(np.mean(np.abs(single[IU] - RN[IU])))

LIONESS = pd.Series(lioness, index=tids)

alt_rows = [
    association(feature, SCORES[feature])[0]
    for feature in ["Coordination Deviation"] + ALT_SCORES
]
alt_rows.append(
    association("LIONESS-style descriptive", LIONESS)[0]
)

ALT_OS = pd.DataFrame(alt_rows)
ALT_OS["q_alternative_family"] = bh(ALT_OS.p)
tab(ALT_OS, "S11_alternative_coordination_complete.csv")