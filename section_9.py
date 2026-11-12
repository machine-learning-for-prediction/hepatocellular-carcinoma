ALL_IDS = dids + nids
PROXIES = pd.DataFrame(index=ALL_IDS)
proxy_defs = []

for name, gl in PROXY_MARKERS.items():
    gg = [
        gene for gene in gl
        if gene in X.index
        and gene not in CURATED
        and np.isfinite(X.loc[gene, ALL_IDS]).all()
        and X.loc[gene, nids].std() > EPS
    ]

    proxy_defs.append({
        "Proxy": name,
        "Genes": ";".join(gg),
        "N_genes": len(gg),
        "Interpretation": "Expression-marker proxy, not cell fraction",
    })

    if len(gg) < 2:
        PROXIES[name] = np.nan
    else:
        values = X.loc[gg, ALL_IDS]
        standardized = values.sub(
            values[nids].mean(axis=1), axis=0
        ).div(values[nids].std(axis=1), axis=0)
        PROXIES[name] = standardized.mean(axis=0)

tab(pd.DataFrame(proxy_defs), "S19_proxy_definitions.csv")
tab(PROXIES.reset_index(names="sample"), "S19_proxy_values.csv")

PROXY_T = pd.DataFrame(index=tids)

for row in proxy_defs:
    gg = row["Genes"].split(";") if row["Genes"] else []
    if len(gg) < 2:
        PROXY_T[row["Proxy"]] = np.nan
        continue
    values = X.loc[gg, tids]
    PROXY_T[row["Proxy"]] = values.sub(
        X.loc[gg, nids].mean(axis=1), axis=0
    ).div(
        X.loc[gg, nids].std(axis=1), axis=0
    ).mean(axis=0)

if PU_ABSOLUTE_LOCAL.exists():
    pu_absolute_path = PU_ABSOLUTE_LOCAL
    INPUTS.append({
        "url": PU_ABSOLUTE_URL,
        "file": str(pu_absolute_path),
        "sha256": sha(pu_absolute_path),
        "retrieved_utc": (
            "User-supplied local copy; download date not established"
        ),
        "source": "GDC Pan-Cancer Atlas ABSOLUTE purity/ploidy file",
    })
else:
    pu_absolute_path = download(
        PU_ABSOLUTE_URL, CACHE / PU_ABSOLUTE_FILENAME
    )

pu_actual_md5 = pu_md5(pu_absolute_path)

if pu_actual_md5 != PU_ABSOLUTE_EXPECTED_MD5:
    raise RuntimeError(
        "The ABSOLUTE file does not match the supplied manifest checksum.\n"
        f"Expected MD5: {PU_ABSOLUTE_EXPECTED_MD5}\n"
        f"Observed MD5: {pu_actual_md5}\n"
        "Check the GDC source before proceeding."
    )

pu_absolute_raw = pd.read_csv(
    pu_absolute_path, sep="\t", low_memory=False
)

tab(
    pd.DataFrame({"Column": pu_absolute_raw.columns.astype(str)}),
    "S34_ABSOLUTE_source_columns.csv",
)

pu_column_normalized = {
    re.sub(r"[^a-z0-9]", "", str(col).lower()): col
    for col in pu_absolute_raw.columns
}
pu_purity_column = pu_column_normalized.get("purity")

if pu_purity_column is None:
    purity_candidates = [
        col for col in pu_absolute_raw.columns
        if re.sub(r"[^a-z0-9]", "", str(col).lower())
        in ["absolutepurity", "tumorpurity", "tumourpurity"]
    ]
    if len(purity_candidates) == 1:
        pu_purity_column = purity_candidates[0]

if pu_purity_column is None:
    raise RuntimeError(
        "Could not identify an unambiguous ABSOLUTE purity column. "
        "See S34_ABSOLUTE_source_columns.csv."
    )

barcode_candidates = []

for col in pu_absolute_raw.columns:
    values = pu_absolute_raw[col].astype(str).str.strip()
    count = int(values.str.match(
        r"^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}"
    ).sum())
    if count:
        barcode_candidates.append((count, str(col)))

if not barcode_candidates:
    raise RuntimeError("No TCGA sample-barcode column in ABSOLUTE file.")

sample_named = [
    col for col in pu_absolute_raw.columns
    if str(col).strip().lower() == "sample"
]
valid_named = (
    bool(sample_named)
    and pu_absolute_raw[sample_named[0]].astype(str).str.strip().str.match(
        r"^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}"
    ).any()
)

if valid_named:
    pu_barcode_column = sample_named[0]
else:
    barcode_candidates.sort(key=lambda item: (-item[0], item[1]))
    pu_barcode_column = barcode_candidates[0][1]

pu_source = pd.DataFrame({
    "Source_barcode": (
        pu_absolute_raw[pu_barcode_column]
        .astype(str).str.strip().str.upper()
    ),
    "ABSOLUTE_purity": pd.to_numeric(
        pu_absolute_raw[pu_purity_column], errors="coerce"
    ),
})

pu_source = pu_source.loc[
    pu_source.Source_barcode.str.match(
        r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}"
    )
].copy()

pu_source["sample15"] = pu_source.Source_barcode.str[:15]
pu_source["patient"] = pu_source.Source_barcode.str[:12]
pu_source["sample_type"] = pu_source.Source_barcode.str[13:15]
pu_source["Valid_fraction"] = (
    pu_source.ABSOLUTE_purity.notna()
    & np.isfinite(pu_source.ABSOLUTE_purity)
    & pu_source.ABSOLUTE_purity.between(0, 1)
)

pu_source_primary = pu_source.loc[
    pu_source.sample_type == "01"
].copy()

tab(
    pu_source_primary.loc[pu_source_primary.sample15.isin(tids)],
    "S34_ABSOLUTE_matching_source_records.csv",
)

purity_mapping_rows = []

for sample in tids:
    subset = pu_source_primary.loc[
        pu_source_primary.sample15 == sample
    ]
    values = subset.loc[
        subset.Valid_fraction, "ABSOLUTE_purity"
    ].to_numpy(float)
    purity_value = np.nan

    if not len(subset):
        status = "No exact selected sample15 match"
    elif not len(values):
        status = "Matched records; no valid purity fraction"
    elif np.max(values) - np.min(values) > 1e-6:
        status = "Conflicting source purities; left missing"
    else:
        purity_value = float(values[0])
        status = "Exact sample15 match"

    purity_mapping_rows.append({
        "sample": sample,
        "patient": pid(sample),
        "ABSOLUTE_purity": purity_value,
        "Mapping_status": status,
        "N_source_records": int(len(subset)),
        "N_valid_source_records": int(len(values)),
        "Source_barcodes": ";".join(subset.Source_barcode.tolist()),
        "Matching_policy": (
            "Exact selected sample15; no patient-level fallback"
        ),
    })

PU_ABSOLUTE_MAP = pd.DataFrame(purity_mapping_rows)
tab(PU_ABSOLUTE_MAP, "S34_ABSOLUTE_purity_mapping.csv")

purity = (
    PU_ABSOLUTE_MAP.set_index("sample")["ABSOLUTE_purity"]
    .reindex(tids).rename("Purity")
)

if int(purity.notna().sum()) < PU_MIN_MATCHED:
    raise RuntimeError(
        "Too few exact selected tumour samples have ABSOLUTE purity.\n"
        f"Matched: {int(purity.notna().sum())}; "
        f"required: {PU_MIN_MATCHED}.\n"
        "Inspect S34_ABSOLUTE_purity_mapping.csv."
    )

pu_disjoint_purity = purity.reindex(dids)

PU_STATUS.append({
    "Analysis": "ABSOLUTE purity mapping",
    "Status": "Completed",
    "N_tumours": len(tids),
    "N_matched_tumours": int(purity.notna().sum()),
    "N_missing_tumours": int(purity.isna().sum()),
    "N_disjoint_tumours": len(dids),
    "N_matched_disjoint_tumours": int(pu_disjoint_purity.notna().sum()),
    "N_high_purity_disjoint_ge_070": int(
        (pu_disjoint_purity >= .70).sum()
    ),
    "Note": (
        "Published DNA-based estimate; exact sample15 matching. "
        "DNA and RNA may derive from different aliquots. "
        "Missing estimates are not imputed."
    ),
})

js({
    "Source": "GDC Pan-Cancer Atlas ABSOLUTE purity/ploidy file",
    "Source_page": "https://gdc.cancer.gov/about-data/publications/pancanatlas",
    "Download_URL": PU_ABSOLUTE_URL,
    "GDC_file_UUID": PU_ABSOLUTE_UUID,
    "Published_filename": PU_ABSOLUTE_FILENAME,
    "Expected_MD5": PU_ABSOLUTE_EXPECTED_MD5,
    "Observed_MD5": pu_actual_md5,
    "Expected_size_bytes": PU_ABSOLUTE_EXPECTED_BYTES,
    "Observed_size_bytes": pu_absolute_path.stat().st_size,
    "Barcode_column": pu_barcode_column,
    "Purity_column": pu_purity_column,
    "Matching_policy": (
        "Exact selected TCGA sample15 barcode; primary tumours only; "
        "no patient-level fallback and no purity imputation"
    ),
    "N_available_tumours": int(purity.notna().sum()),
    "N_available_disjoint_tumours": int(pu_disjoint_purity.notna().sum()),
}, "ABSOLUTE_purity_provenance.json")

print(
    "\nABSOLUTE purity mapping complete:",
    int(purity.notna().sum()), "of", len(tids), "selected tumours;",
    int(pu_disjoint_purity.notna().sum()), "of", len(dids),
    "patient-disjoint tumours.",
)

if PU_ADD_ESTIMATE:
    try:
        if shutil.which("Rscript") is None:
            print("\nInstalling R for ESTIMATE...")
            subprocess.run(["apt-get", "update", "-qq"], check=True)
            subprocess.run([
                "apt-get", "install", "-y", "-qq",
                "r-base-core", "r-base-dev",
                "libcurl4-openssl-dev", "libssl-dev", "libxml2-dev",
            ], check=True)

        estimate_samples = tids + nids

        if len(set(estimate_samples)) != len(estimate_samples):
            raise RuntimeError("Duplicate ESTIMATE input sample IDs.")

        estimate_input = META / "ESTIMATE_full_expression_input.tsv"
        X.loc[:, estimate_samples].to_csv(
            estimate_input, sep="\t", index=True, index_label="Gene"
        )

        estimate_output = TABLES / "S35_ESTIMATE_scores.csv"
        estimate_gene_output = TABLES / "S35_ESTIMATE_retained_common_genes.csv"
        estimate_session = META / "ESTIMATE_R_sessionInfo.txt"
        estimate_r = META / "ESTIMATE_analysis.R"
        estimate_log = META / "ESTIMATE_R_execution.log"

        r_code = r"""
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("tidyestimate", quietly = TRUE)) {
    install.packages("tidyestimate", dependencies = NA, Ncpus = 2)
}
if (!requireNamespace("tidyestimate", quietly = TRUE)) {
    stop("tidyestimate installation failed.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
    stop("Expected input, score output, gene output and session output.")
}

x <- read.delim(
    args[1],
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

if (anyDuplicated(rownames(x))) {
    stop("Duplicate gene identifiers in ESTIMATE input.")
}
if (anyDuplicated(colnames(x))) {
    stop("Duplicate sample identifiers in ESTIMATE input.")
}
if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("Non-numeric sample column in ESTIMATE input.")
}

filtered <- tidyestimate::filter_common_genes(
    x,
    id = "hgnc_symbol",
    tidy = FALSE,
    tell_missing = TRUE,
    find_alias = FALSE
)

if (nrow(filtered) < 6000L) {
    stop("Fewer than 6000 common genes matched.")
}

expression_values <- as.matrix(filtered[, -1, drop = FALSE])
if (any(!is.finite(expression_values))) {
    stop("Non-finite retained ESTIMATE expression values.")
}

scores <- tidyestimate::estimate_score(
    filtered,
    is_affymetrix = FALSE
)

required_columns <- c("sample", "stromal", "immune", "estimate")
if (!all(required_columns %in% names(scores))) {
    stop(
        "Unexpected tidyestimate output columns: ",
        paste(names(scores), collapse = ", ")
    )
}

scores <- scores[, required_columns, drop = FALSE]
write.csv(scores, args[2], row.names = FALSE)

write.csv(
    data.frame(Gene = as.character(filtered[[1]])),
    args[3],
    row.names = FALSE
)

capture.output(
    {
        cat(
            "tidyestimate version:",
            as.character(utils::packageVersion("tidyestimate")),
            "\n"
        )
        cat("is_affymetrix: FALSE\n")
        cat("Gene aliases matched: FALSE\n")
        cat("Retained common genes:", nrow(filtered), "\n")
        print(sessionInfo())
    },
    file = args[4]
)

cat("\nESTIMATE completed.\n")
cat("Retained common genes:", nrow(filtered), "\n")
cat("Scored samples:", nrow(scores), "\n")
"""
        estimate_r.write_text(r_code, encoding="utf-8")

        with open(estimate_log, "w") as handle:
            result = subprocess.run([
                "Rscript",
                str(estimate_r),
                str(estimate_input),
                str(estimate_output),
                str(estimate_gene_output),
                str(estimate_session),
            ], stdout=handle, stderr=subprocess.STDOUT, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                "ESTIMATE R calculation failed.\n"
                + estimate_log.read_text(errors="replace")[-12000:]
            )

        estimate_df = pd.read_csv(
            estimate_output, dtype={"sample": str}
        )

        if estimate_df["sample"].duplicated().any():
            raise RuntimeError("Duplicate ESTIMATE output sample IDs.")

        estimate_df = estimate_df.set_index("sample").rename(columns={
            "stromal": "ESTIMATE Stromal",
            "immune": "ESTIMATE Immune",
            "estimate": "ESTIMATE Score",
        })

        estimate_names = [
            "ESTIMATE Stromal", "ESTIMATE Immune", "ESTIMATE Score"
        ]

        for name in estimate_names:
            estimate_df[name] = pd.to_numeric(
                estimate_df[name], errors="coerce"
            )

        PU_ESTIMATE = estimate_df.reindex(estimate_samples)[estimate_names]

        if not np.isfinite(PU_ESTIMATE.to_numpy()).all():
            raise RuntimeError("Missing or non-finite ESTIMATE scores.")

        annotated = PU_ESTIMATE.reset_index(names="sample")
        annotated["patient"] = annotated["sample"].map(pid)
        annotated["Tissue"] = np.where(
            annotated["sample"].isin(tids),
            "tumour", "adjacent non-tumour",
        )
        annotated["Interpretation"] = (
            "RNA-seq relative enrichment score; "
            "not a calibrated tumour-cell fraction"
        )

        tab(annotated, "S35_ESTIMATE_scores.csv")
        PU_ESTIMATE_AVAILABLE = True

        PU_STATUS.append({
            "Analysis": "ESTIMATE score calculation",
            "Status": "Completed",
            "N_tumours": len(tids),
            "N_adjacent": len(nids),
            "Note": (
                "tidyestimate; common-gene filtering; "
                "is_affymetrix=False; no purity-fraction conversion"
            ),
        })

        js({
            "Implementation": "R package tidyestimate",
            "Input": "Full available Xena expression for selected samples",
            "Common_gene_filtering": True,
            "Gene_identifier_type": "HGNC symbols",
            "Alias_matching": False,
            "is_affymetrix": False,
            "Purity_conversion": "Not performed",
            "Retained_common_genes_file": (
                "S35_ESTIMATE_retained_common_genes.csv"
            ),
            "R_session_file": str(estimate_session),
            "R_log_file": str(estimate_log),
            "Documentation": (
                "https://search.r-project.org/CRAN/refmans/"
                "tidyestimate/html/estimate_score.html"
            ),
        }, "ESTIMATE_provenance.json")

        print("ESTIMATE completed for", len(PU_ESTIMATE), "samples.")

    except Exception as error:
        PU_STATUS.append({
            "Analysis": "ESTIMATE score calculation",
            "Status": "Failed",
            "Reason": str(error),
        })
        tab(
            pd.DataFrame(PU_STATUS),
            "S39_purity_ESTIMATE_analysis_status.csv",
        )
        NOTICES.append("ESTIMATE addition failed: " + str(error))
        if PU_REQUIRE_ESTIMATE:
            raise
else:
    PU_STATUS.append({
        "Analysis": "ESTIMATE score calculation",
        "Status": "Disabled by user",
    })

tab(pd.DataFrame(PU_STATUS), "S39_purity_ESTIMATE_analysis_status.csv")

covs = PROXY_T.join(
    CS[["stage_numeric", "grade_numeric"]]
).join(purity)

cov_associations = []

for feature in FOLLOWUP_FEATURES:
    values = feature_values(feature)
    for covariate in covs:
        rho, pvalue, count = sp(values, covs[covariate])
        cov_associations.append({
            "Feature": feature,
            "Covariate": covariate,
            "rho": rho,
            "p": pvalue,
            "N": count,
        })

ca = pd.DataFrame(cov_associations)
ca["q_exploratory_family"] = bh(ca.p)
tab(ca, "S21_score_confounder_associations.csv")


def residual_corr(values, covariates):
    ranked = stats.rankdata(np.asarray(values, float), axis=1).T
    ranked_covariates = stats.rankdata(
        np.asarray(covariates, float), axis=0
    )
    design = np.column_stack([
        np.ones(len(ranked_covariates)), ranked_covariates
    ])
    residual = (
        ranked
        - design @ np.linalg.lstsq(design, ranked, rcond=None)[0]
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.corrcoef(residual.T)


adjust_names = [
    "Immune proxy", "Stromal proxy", "Hepatocyte differentiation proxy"
]

if PROXIES[adjust_names].notna().all().all():
    ct = PROXIES.loc[dids, adjust_names].to_numpy()
    cn = PROXIES.loc[nids, adjust_names].to_numpy()
    art = residual_corr(TD.to_numpy(), ct)
    arn = residual_corr(N.to_numpy(), cn)
    aad = art - arn

    rng = np.random.default_rng(46)
    aboot = []

    for _ in tqdm(
        range(N_ADJUSTED_BOOT),
        desc="Expression-proxy adjusted-correlation sensitivity",
    ):
        it = rng.integers(0, len(dids), len(dids))
        inn = rng.integers(0, len(nids), len(nids))
        aboot.append((
            residual_corr(TD.to_numpy()[:, it], ct[it])
            - residual_corr(N.to_numpy()[:, inn], cn[inn])
        )[IU])

    aboot = np.asarray(aboot)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        al, ah = np.nanpercentile(aboot, [2.5, 97.5], axis=0)

    out = []

    for a, b in ALL_EDGE_PAIRS:
        k = edge_key[frozenset((a, b))]
        i, j = genes.index(a), genes.index(b)
        out.append({
            "Gene1": a,
            "Gene2": b,
            "Unadjusted_delta": D[i, j],
            "Adjusted_partial_rank_delta": aad[i, j],
            "Adjusted_boot_CI_low": al[k],
            "Adjusted_boot_CI_high": ah[k],
            "Target_policy": "All retained pairs",
            "Note": (
                "Pointwise sensitivity CI; no BW test on adjusted residuals"
            ),
        })

    tab(pd.DataFrame(out), "S22_adjusted_edge_sensitivity.csv")

high = [
    sample for sample in dids
    if pd.notna(purity[sample]) and purity[sample] >= .70
]

if len(high) >= 30 and FOLLOWUP_EDGE_PAIRS:
    rh = corr(T[high])
    delta = rh - RN
    hl, hh, hv = edge_intervals(
        T[high].to_numpy(),
        N.to_numpy(),
        N_ADJUSTED_BOOT,
        47,
        pair_indices=selected_iu,
    )

    out = []

    for selected_k, (a, b) in enumerate(FOLLOWUP_EDGE_PAIRS):
        i, j = genes.index(a), genes.index(b)
        enough_valid = hv[selected_k] >= max(
            20, int(np.ceil(.90 * N_ADJUSTED_BOOT))
        )

        out.append({
            "Gene1": a,
            "Gene2": b,
            "N_high_purity": len(high),
            "Delta_high_purity": delta[i, j],
            "CI_low": hl[selected_k] if enough_valid else np.nan,
            "CI_high": hh[selected_k] if enough_valid else np.nan,
            "Bootstrap_valid": hv[selected_k],
            "Target_policy": EDGE_FOLLOWUP_POLICY,
            "Interpretation": (
                "Same-cohort post-selection sensitivity; pointwise CI"
            ),
        })

    tab(pd.DataFrame(out), "S23_high_purity_edges.csv")

proxy_adjust = []

for extra in [
    ["Hepatocyte differentiation proxy", "Proliferation proxy"],
    ["Immune proxy", "Stromal proxy"],
    ["Purity"],
]:
    design = BASE_CLIN.join(covs[extra]).replace(
        [np.inf, -np.inf], np.nan
    ).dropna()

    if len(design) < 60 or Y.loc[design.index, "E"].sum() < 15:
        continue

    for feature in MODEL_FEATURES:
        row, _ = association(
            feature,
            feature_values(feature),
            design,
            "Clinical_plus_" + "_".join(extra),
        )
        proxy_adjust.append(row)

if proxy_adjust:
    pa = pd.DataFrame(proxy_adjust)
    pa["q_sensitivity_family_postselection_exploratory"] = bh(pa.p)
    tab(pa, "S24_confounder_adjusted_OS_sensitivity.csv")