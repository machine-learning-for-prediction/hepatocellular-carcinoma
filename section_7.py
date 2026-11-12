def stage_number(value):
    text = str(value).upper().strip()
    text = re.sub(r"^STAGE\s*", "", text)
    if text in ["1", "2", "3", "4"]:
        return float(text)
    match = re.match(r"^(IV|III|II|I)(?:[ABC]|\s|$)", text)
    return {"I": 1., "II": 2., "III": 3., "IV": 4.}.get(
        match.group(1) if match else "", np.nan
    )


def grade_number(value):
    match = re.search(r"(?:G|GRADE\s*)([1-4])", str(value).upper())
    return float(match.group(1)) if match else np.nan


def cirrhosis_binary(value):
    if pd.isna(value):
        return np.nan
    text = str(value).upper().strip()
    if text in ["NAN", "NONE", ""]:
        return np.nan
    if text in ["6", "6.0"] or ("6" in text and "CIRRHOSIS" in text):
        return 1.0
    if text in ["0", "1", "2", "3", "4", "5"]:
        return 0.0
    if any(item in text for item in [
        "0 -", "1,2", "3,4", "5 -", "NO FIBROSIS",
        "PORTAL", "SEPTA", "NODULAR",
    ]):
        return 0.0
    return np.nan


def any_treatment(row):
    yes_values = {"YES", "Y", "1"}
    no_values = {"NO", "N", "0"}
    saw_yes = False
    saw_no = False

    for col in [
        "radiation_therapy",
        "additional_pharmaceutical_therapy",
        "post_op_ablation_embolization_tx",
        "postoperative_rx_tx",
        "history_of_neoadjuvant_treatment",
    ]:
        value = str(row[col]).upper().strip()
        saw_yes |= value in yes_values
        saw_no |= value in no_values

    if saw_yes:
        return 1.0
    if saw_no:
        return 0.0
    return np.nan


clinical_audit = []
C = pd.DataFrame(index=[pid(sample) for sample in tids])

try:
    cp = download(CLIN_URL, CACHE / "LIHC_clinicalMatrix")
    cc = pd.read_csv(cp, sep="\t", index_col=0)
    cc.index = cc.index.astype(str).str[:12]

    tab(
        pd.DataFrame({"Column": cc.columns}),
        "clinical_source_columns.csv",
    )

    candidates = {
        "age": ["age_at_initial_pathologic_diagnosis"],
        "sex": ["gender"],
        "stage": ["pathologic_stage"],
        "grade": ["neoplasm_histologic_grade"],
        "AFP": ["fetoprotein_outcome_value"],
        "cirrhosis": ["fibrosis_ishak_score"],
        "radiation_therapy": ["radiation_therapy"],
        "additional_pharmaceutical_therapy": [
            "additional_pharmaceutical_therapy"
        ],
        "post_op_ablation_embolization_tx": [
            "post_op_ablation_embolization_tx"
        ],
        "postoperative_rx_tx": ["postoperative_rx_tx"],
        "history_of_neoadjuvant_treatment": [
            "history_of_neoadjuvant_treatment"
        ],
    }

    mapped = pd.DataFrame(index=cc.index)

    for target, possible in candidates.items():
        source = next((col for col in possible if col in cc.columns), None)
        clinical_audit.append({"Target": target, "Source": source})
        mapped[target] = cc[source] if source else np.nan

    rows = []

    for patient, subset in mapped.groupby(level=0):
        result = {"patient": patient}
        for col in mapped.columns:
            values = subset[col].dropna().astype(str)
            values = values[~values.str.lower().isin([
                "nan", "na", "[not available]", "[not applicable]", ""
            ])].unique()
            result[col] = values[0] if len(values) == 1 else np.nan
        rows.append(result)

    C = pd.DataFrame(rows).set_index("patient").reindex(C.index)

except Exception as error:
    NOTICES.append("Clinical download/mapping unavailable: " + str(error))

for col in [
    "age", "sex", "stage", "grade", "AFP", "cirrhosis",
    "radiation_therapy", "additional_pharmaceutical_therapy",
    "post_op_ablation_embolization_tx", "postoperative_rx_tx",
    "history_of_neoadjuvant_treatment",
]:
    if col not in C.columns:
        C[col] = np.nan

C["age"] = pd.to_numeric(C["age"], errors="coerce")
C["sex_binary"] = C["sex"].astype(str).str.lower().map({
    "female": 0., "male": 1.
})
C["stage_numeric"] = C["stage"].map(stage_number)
C["grade_numeric"] = C["grade"].map(grade_number)
C["AFP"] = pd.to_numeric(C["AFP"], errors="coerce")
C["cirrhosis"] = C["cirrhosis"].map(cirrhosis_binary)
C["treatment"] = C.apply(any_treatment, axis=1)

tab(C.reset_index(names="patient"), "clinical_patient_table.csv")
tab(
    pd.DataFrame(clinical_audit, columns=["Target", "Source"]),
    "clinical_mapping_audit.csv",
)

tab(pd.DataFrame([
    {
        "Variable": col,
        "N_available": int(C[col].notna().sum()),
        "N_missing": int(C[col].isna().sum()),
    }
    for col in [
        "age", "sex_binary", "stage_numeric", "grade_numeric",
        "AFP", "cirrhosis", "treatment",
    ]
]), "S16_clinical_missingness.csv")

CS = C.copy()
CS.index = tids


def clinical_design(frame, model="base"):
    required = ["age", "sex_binary", "stage_numeric"]
    if model in ["extended", "full"]:
        required.append("grade_numeric")
    if model == "full":
        required += ["AFP", "cirrhosis", "treatment"]

    values = frame[required].copy()
    valid_cases = values.notna().all(axis=1)
    categorical = pd.Categorical(
        values.loc[valid_cases, "stage_numeric"],
        categories=[1., 2., 3., 4.],
    )
    dummies = pd.get_dummies(
        categorical, prefix="stage", drop_first=True, dtype=float
    )
    dummies.index = values.index[valid_cases]

    result = values.loc[valid_cases, ["age", "sex_binary"]].join(dummies)

    if model in ["extended", "full"]:
        result["grade"] = values.loc[valid_cases, "grade_numeric"]

    if model == "full":
        for col in ["AFP", "cirrhosis", "treatment"]:
            result[col] = values.loc[valid_cases, col]

    return result.loc[:, result.nunique() > 1]


BASE_CLIN = clinical_design(CS, "base")
EXT_CLIN = clinical_design(CS, "extended")
FULL_CLIN = clinical_design(CS, "full")

adjusted = []

for model_name, base in [
    ("Age_sex_stage", BASE_CLIN),
    ("Age_sex_stage_grade", EXT_CLIN),
    ("Age_sex_stage_grade_AFP_cirrhosis_treatment", FULL_CLIN),
]:
    for feature in MODEL_FEATURES:
        row, _ = association(
            feature, feature_values(feature), base, model_name
        )
        adjusted.append(row)

ADJ_OS = pd.DataFrame(adjusted)

if not ADJ_OS.empty:
    ADJ_OS["q_within_model_postselection_exploratory"] = (
        ADJ_OS.groupby("Model")["p"].transform(
            lambda values: bh(values.to_numpy())
        )
    )

tab(ADJ_OS, "S17_adjusted_OS.csv")
