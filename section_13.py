external_status = []

if DOWNLOAD_EXTERNAL:
    try:
        geo = GEOparse.get_GEO(
            geo=EXTERNAL_GSE,
            destdir=str(CACHE),
            annotate_gpl=True,
        )

        metadata_rows = []

        for gsm, obj in geo.gsms.items():
            row = {"sample": gsm}
            for key, values in obj.metadata.items():
                row[key] = " | ".join(map(str, values))
            metadata_rows.append(row)

        tab(
            pd.DataFrame(metadata_rows),
            "external_GEO_metadata_to_review.csv",
        )

        matrices, annotation_audit = [], []

        for gplid, gpl in geo.gpls.items():
            samples = [
                sample for sample, obj in geo.gsms.items()
                if gplid in obj.metadata.get("platform_id", [])
            ]
            if not samples:
                continue

            annotation = gpl.table.copy()
            symbol_col = next((
                col for col in [
                    "Gene Symbol", "GENE_SYMBOL", "Symbol",
                    "SYMBOL", "gene_symbol",
                ]
                if col in annotation.columns
            ), None)
            id_col = next((
                col for col in ["ID", "ID_REF"]
                if col in annotation.columns
            ), None)

            annotation_audit.append({
                "Platform": gplid,
                "Symbol_column": symbol_col,
                "ID_column": id_col,
                "Columns": ";".join(map(str, annotation.columns)),
            })

            if symbol_col is None or id_col is None:
                continue

            symbols = (
                annotation[symbol_col].fillna("").astype(str).str.strip()
            )
            good = (
                symbols.str.match(r"^[A-Za-z0-9_-]+$")
                & ~symbols.isin(["", "NA", "nan", "---"])
            )
            mapping = pd.Series(
                symbols[good].values,
                index=annotation.loc[good, id_col].astype(str),
            )
            mapping = mapping.loc[~mapping.index.duplicated()]

            expression_columns = {}

            for sample in samples:
                sample_table = geo.gsms[sample].table
                if not {"ID_REF", "VALUE"}.issubset(sample_table.columns):
                    continue

                values = pd.Series(
                    pd.to_numeric(
                        sample_table.VALUE, errors="coerce"
                    ).values,
                    index=sample_table.ID_REF.astype(str),
                )
                expression_columns[sample] = values.loc[
                    ~values.index.duplicated()
                ]

            if not expression_columns:
                continue

            expression = pd.DataFrame(expression_columns)
            expression = expression.loc[
                expression.index.intersection(mapping.index)
            ]
            expression["gene"] = mapping.reindex(expression.index)
            expression = expression.groupby("gene").median()
            matrices.append((gplid, expression))

        tab(
            pd.DataFrame(annotation_audit),
            "external_annotation_audit.csv",
        )

        if not matrices:
            raise RuntimeError("No usable automatically mapped GEO platform.")

        platform_id, EX = max(
            matrices, key=lambda item: item[1].shape[1]
        )

        js({
            "GSE": EXTERNAL_GSE,
            "Platform": platform_id,
            "Expression_quantiles": np.nanpercentile(
                EX.to_numpy(), [0, 1, 50, 99, 100]
            ),
            "Applied_log2": EXTERNAL_LOG2_TRANSFORM,
            "Note": "Inspect GEO processing metadata before transformation.",
        }, "external_expression_scale.json")

        if EXTERNAL_LOG2_TRANSFORM:
            if np.nanmin(EX.to_numpy()) < 0:
                raise RuntimeError("Cannot log2-transform negative intensities.")
            EX = np.log2(EX + 1)

        tab(pd.DataFrame({
            "sample": EX.columns,
            "patient": "",
            "tissue": "",
            "OS_time": "",
            "OS_event": "",
        }), "external_metadata_TEMPLATE.csv")

        if not EXTERNAL_METADATA_VERIFIED or not EXT_META_FILE.exists():
            external_status.append({
                "Status": "Downloaded, inference not performed",
                "Reason": "Verified sample/patient/tissue metadata required",
            })
        else:
            em = pd.read_csv(
                EXT_META_FILE, dtype={"sample": str, "patient": str}
            )
            required = {
                "sample", "patient", "tissue", "OS_time", "OS_event"
            }
            if not required.issubset(em.columns):
                raise RuntimeError("External metadata missing required columns.")
            if em["sample"].duplicated().any():
                raise RuntimeError("Duplicate external samples.")

            unmatched_external_samples = em.loc[
                ~em["sample"].isin(EX.columns), "sample"
            ].tolist()
            em = em.loc[em["sample"].isin(EX.columns)].copy()

            if (
                em.patient.isna().any()
                or (em.patient.str.strip() == "").any()
            ):
                raise RuntimeError("Verified external patient IDs required.")
            if not em.tissue.isin(["tumour", "adjacent"]).all():
                raise RuntimeError("External tissue must be tumour/adjacent.")
            if em.duplicated(["patient", "tissue"]).any():
                raise RuntimeError(
                    "Multiple external samples per patient/tissue."
                )

            em = em.set_index("sample")
            et = em.index[em.tissue == "tumour"].tolist()
            en = em.index[em.tissue == "adjacent"].tolist()
            adjacent_patients = set(em.loc[en, "patient"])
            edt = [
                sample for sample in et
                if em.loc[sample, "patient"] not in adjacent_patients
            ]

            eg = [
                gene for gene in genes
                if gene in EX.index
                and np.isfinite(EX.loc[gene, et + en]).all()
                and EX.loc[gene, en].std() > EPS
            ]

            tab(pd.DataFrame({
                "Gene": genes,
                "Available_external": [gene in eg for gene in genes],
            }), "S30_external_gene_coverage.csv")

            external_status.append({
                "Status": "Verified metadata loaded",
                "Platform": platform_id,
                "N_tumour": len(et),
                "N_adjacent": len(en),
                "N_disjoint_tumour": len(edt),
                "N_available_genes": len(eg),
                "N_metadata_samples_not_on_selected_platform": len(
                    unmatched_external_samples
                ),
                "N_discovery_rewired_target_pairs": len(FOLLOWUP_EDGE_PAIRS),
            })

            if len(edt) >= 20 and len(en) >= 10 and FOLLOWUP_EDGE_PAIRS:
                external_pairs = [
                    (a, b) for a, b in FOLLOWUP_EDGE_PAIRS
                    if a in eg and b in eg
                ]
                eb_rows = []

                if external_pairs:
                    external_genes = [
                        gene for gene in genes
                        if any(gene in pair for pair in external_pairs)
                    ]
                    external_index = {
                        gene: i for i, gene in enumerate(external_genes)
                    }
                    external_pair_indices = (
                        np.array([
                            external_index[a] for a, b in external_pairs
                        ], dtype=int),
                        np.array([
                            external_index[b] for a, b in external_pairs
                        ], dtype=int),
                    )

                    external_tvalues = EX.loc[
                        external_genes, edt
                    ].to_numpy()
                    external_nvalues = EX.loc[
                        external_genes, en
                    ].to_numpy()
                    external_rt = corr(external_tvalues)
                    external_rn = corr(external_nvalues)

                    elow, ehigh, evalid = edge_intervals(
                        external_tvalues,
                        external_nvalues,
                        N_EDGE_BOOT,
                        51,
                        pair_indices=external_pair_indices,
                    )

                    for k, (a, b) in enumerate(external_pairs):
                        ei, ej = external_index[a], external_index[b]
                        i, j = genes.index(a), genes.index(b)
                        rt, rn = external_rt[ei, ej], external_rn[ei, ej]
                        enough_valid = evalid[k] >= max(
                            20, int(np.ceil(.90 * N_EDGE_BOOT))
                        )

                        eb_rows.append({
                            "Gene1": a,
                            "Gene2": b,
                            "External_rho_adjacent": rn,
                            "External_rho_tumour": rt,
                            "External_delta": rt - rn,
                            "External_CI_low": (
                                elow[k] if enough_valid else np.nan
                            ),
                            "External_CI_high": (
                                ehigh[k] if enough_valid else np.nan
                            ),
                            "External_bootstrap_valid": evalid[k],
                            "External_p": float(
                                bw(rt, rn, len(edt), len(en))
                            ),
                            "TCGA_delta": D[i, j],
                            "Same_direction": (
                                np.sign(rt - rn) == np.sign(D[i, j])
                            ),
                            "N_tumour_disjoint": len(edt),
                            "N_adjacent": len(en),
                            "Target_policy": (
                                "Primary TCGA rewired pairs with "
                                "available external expression"
                            ),
                        })

                eb = pd.DataFrame(eb_rows)

                if not eb.empty:
                    eb["q_among_available_discovery_rewired_pairs"] = bh(
                        eb.External_p
                    )
                    eb["External_BW_and_CI_supported"] = (
                        eb.q_among_available_discovery_rewired_pairs.lt(.05)
                        & (
                            eb.External_CI_low.gt(0)
                            | eb.External_CI_high.lt(0)
                        )
                    )
                    eb["Directionally_supported_replication"] = (
                        eb.External_BW_and_CI_supported & eb.Same_direction
                    )

                tab(eb, "S31_external_edge_replication.csv")
                external_status.append({
                    "Status": "Discovery-rewired external edge testing completed",
                    "N_available_target_pairs": len(external_pairs),
                })
            else:
                external_status.append({
                    "Status": "External edge inference not performed",
                    "Reason": (
                        "Insufficient patient-disjoint samples "
                        "or no primary TCGA rewired pairs"
                    ),
                })

            if len(et) >= 30 and len(en) >= 10 and len(eg) >= 6:
                er = fit_ref(EX.loc[eg, et], EX.loc[eg, en])
                es = make_scores(EX.loc[eg, et], er)

                ey = pd.DataFrame(index=et)
                ey["T"] = pd.to_numeric(
                    em.loc[et, "OS_time"], errors="coerce"
                )
                ey["E"] = pd.to_numeric(
                    em.loc[et, "OS_event"], errors="coerce"
                )
                ey = ey.loc[
                    ey["T"].notna()
                    & np.isfinite(ey["T"])
                    & ey["T"].gt(0)
                    & ey["E"].isin([0, 1])
                ]

                erows, feature_coverage_rows = [], []
                reference_genes = set(er["genes"])

                for feature in FOLLOWUP_FEATURES:
                    if feature in genes:
                        needed = [feature]
                        eligible = feature in eg
                        reason = (
                            "Available"
                            if eligible else "Gene unavailable externally"
                        )
                    else:
                        needed = SCORE_GENE_REQUIREMENTS[feature]
                        eligible = (
                            len(needed) > 0
                            and set(needed).issubset(reference_genes)
                        )
                        reason = (
                            "Complete retained constituent coverage"
                            if eligible
                            else "Incomplete retained constituent/reference coverage"
                        )

                    feature_coverage_rows.append({
                        "Feature": feature,
                        "Eligible": eligible,
                        "Required_discovery_genes": ";".join(needed),
                        "Missing_reference_genes": ";".join(
                            gene for gene in needed
                            if gene not in reference_genes
                        ),
                        "Selected_for_TCGA_adjusted_Cox_CV": (
                            feature in MODEL_FEATURES
                        ),
                        "Reason": reason,
                    })

                    if not eligible:
                        continue

                    values = (
                        EX.loc[feature, ey.index]
                        if feature in genes else es.loc[ey.index, feature]
                    )
                    row, _ = association(
                        feature, values, name="External association", y=ey
                    )
                    row["Reference_policy"] = (
                        "External outcome-independent reference"
                    )
                    row["Target_policy"] = (
                        "All eligible main scores and retained genes; "
                        "external univariable OS not restricted by "
                        "the TCGA adjusted-Cox/CV shortlist"
                    )
                    erows.append(row)

                tab(
                    pd.DataFrame(feature_coverage_rows),
                    "S32_external_feature_eligibility.csv",
                )

                external_os = pd.DataFrame(erows)
                if not external_os.empty:
                    external_os["q_external_target_family"] = bh(
                        external_os.p
                    )
                tab(external_os, "S32_external_survival_associations.csv")

            INPUTS.append({
                "file": str(EXT_META_FILE),
                "sha256": sha(EXT_META_FILE),
            })

        for cached_geo_file in CACHE.glob(EXTERNAL_GSE + "*"):
            if cached_geo_file.is_file():
                INPUTS.append({
                    "file": str(cached_geo_file),
                    "sha256": sha(cached_geo_file),
                    "source": "GEOparse cached GEO file",
                })

    except Exception as error:
        external_status.append({
            "Status": "Unavailable/error",
            "Reason": str(error),
        })
        NOTICES.append("External cohort: " + str(error))
else:
    external_status.append({"Status": "Disabled by user"})

tab(
    pd.DataFrame(external_status),
    "external_validation_status.csv",
)

js({
    "ABSOLUTE_available_tumours": int(purity.notna().sum()),
    "ABSOLUTE_missing_tumours": int(purity.isna().sum()),
    "ABSOLUTE_available_disjoint_tumours": int(
        purity.reindex(dids).notna().sum()
    ),
    "High_purity_disjoint_tumours_ge_070": pu_high_count,
    "High_purity_file_created": (
        TABLES / "S23_high_purity_edges.csv"
    ).exists(),
    "ESTIMATE_available": PU_ESTIMATE_AVAILABLE,
    "ESTIMATE_adjusted_edge_file_created": (
        TABLES / "S37_ESTIMATE_adjusted_edge_sensitivity.csv"
    ).exists(),
    "ESTIMATE_OS_file_created": (
        TABLES / "S38_ESTIMATE_adjusted_OS_sensitivity.csv"
    ).exists(),
    "Additional_edge_bootstrap_seed": 260101,
    "Additional_edge_bootstrap_resamples": PU_N_ESTIMATE_BOOT,
    "Adjusted_Cox_CV_targets": MODEL_FEATURES,
    "Selected_edge_followup_targets": FOLLOWUP_EDGE_PAIRS,
    "Feature_selection_policy": MODEL_SELECTION_POLICY,
    "Edge_followup_policy": EDGE_FOLLOWUP_POLICY,
    "Outcome_based_feature_selection": True,
    "Association_based_edge_followup_selection": True,
    "CV_interpretation": CV_INTERPRETATION,
    "Status_records": PU_STATUS,
}, "purity_ESTIMATE_results_summary.json")

versions = {"Python": platform.python_version()}

for package in PACKAGES.values():
    try:
        versions[package] = version(package)
    except Exception:
        versions[package] = "unknown"

summary = {
    "N_tumour": len(tids),
    "N_adjacent": len(nids),
    "N_disjoint_tumour": len(dids),
    "N_pairs": len(paired),
    "N_events": int(Y.E.sum()),
    "N_genes": len(genes),
    "Unavailable_genes": [gene for gene in CURATED if gene not in genes],
    "Formal_disjoint_DE_up": int((
        DE.q_disjoint.lt(.05)
        & DE.Mean_difference_disjoint.gt(0)
    ).sum()),
    "Formal_disjoint_DE_down": int((
        DE.q_disjoint.lt(.05)
        & DE.Mean_difference_disjoint.lt(0)
    ).sum()),
    "N_rewired": int(EDGES.Rewired.sum()),
    "Global_mean_abs_delta": obs_global,
    "All_initial_features": FOLLOWUP_FEATURES,
    "Top_scores": TOP_SCORES,
    "Top_genes": TOP_GENES,
    "Adjusted_Cox_CV_targets": MODEL_FEATURES,
    "N_all_discovery_pairs": len(ALL_EDGE_PAIRS),
    "N_selected_edge_followup_pairs": len(FOLLOWUP_EDGE_PAIRS),
    "Selected_edge_followup_targets": FOLLOWUP_EDGE_PAIRS,
    "Outcome_based_feature_selection": True,
    "Association_based_edge_followup_selection": True,
    "Feature_selection_policy": MODEL_SELECTION_POLICY,
    "Edge_followup_policy": EDGE_FOLLOWUP_POLICY,
    "CV_interpretation": CV_INTERPRETATION,
    "Coordination_imputed_N": int(
        RAW_SCORES["Coordination Deviation"].isna().sum()
    ),
    "Burden_identity_max_error": float(np.max(np.abs(
        SCORES["Activation Burden"]
        - .5 * (
            SCORES["Signed Burden"] + SCORES["Dysregulation Burden"]
        )
    ))),
    "Successful_penalized_fits": sum(
        record.get("Success", False)
        and record.get("Penalty", 0) > 0
        for record in LOG
    ),
}

js(summary, "results_summary.json")

js({
    "UTC": datetime.now(timezone.utc).isoformat(),
    "Dataset_identifiers": {
        "Expression": EXPR_ID,
        "Survival": SURV_ID,
        "Clinical": CLIN_ID,
    },
    "Inputs": INPUTS,
    "Versions": versions,
    "Model_log": LOG,
    "CV_audit": CV_LOG,
    "Duplicate_gene_rows_removed": duplicate_gene_rows,
    "Selected_tumours": tids,
    "Selected_adjacent": nids,
    "All_initial_features": FOLLOWUP_FEATURES,
    "Adjusted_Cox_CV_targets": MODEL_FEATURES,
    "All_discovery_edge_pairs": ALL_EDGE_PAIRS,
    "Selected_edge_followup_targets": FOLLOWUP_EDGE_PAIRS,
    "Selection_policy": {
        "Outcome_based_feature_selection": True,
        "Association_based_edge_followup_selection": True,
        "Feature_policy": MODEL_SELECTION_POLICY,
        "Edge_policy": EDGE_FOLLOWUP_POLICY,
        "Restricted_edge_analyses": [
            "Balanced-size sensitivity",
            "High-purity edge sensitivity",
            "ESTIMATE-adjusted edge bootstrap",
            "External edge testing",
        ],
        "Expression_proxy_adjusted_edge_policy": "All retained pairs",
        "Balanced_BH_family": "All original retained gene pairs",
        "External_edge_BH_family": (
            "Available discovery-rewired pairs"
        ),
        "External_univariable_OS_policy": (
            "All eligible main scores and retained genes"
        ),
        "CV_note": CV_INTERPRETATION,
        "Adjusted_Cox_inference": (
            "Exploratory same-cohort post-selection inference; "
            "shortlist BH adjustment does not correct selection bias"
        ),
    },
    "Parameters": {
        "Seed": SEED,
        "Edge_boot": N_EDGE_BOOT,
        "Network_permutations": N_NETWORK_PERM,
        "Graph_permutations": N_GRAPH_PERM,
        "Reference_boot": N_REFERENCE_BOOT,
        "Balanced_repeats": N_BALANCED_REPEATS,
        "Balanced_boot_per_repeat": N_BALANCED_BOOT,
        "Adjusted_boot": N_ADJUSTED_BOOT,
        "Simulation_repeats": N_SIM_REPEATS,
        "CV_seeds": CV_SEEDS,
        "CV_folds": N_FOLDS,
        "Ridge_fallbacks": RIDGES,
        "Coordination_z": COORD_Z,
        "Coordination_rho": COORD_RHO,
        "Coordination_minpairs": COORD_PAIRS,
        "Top_n_scores": TOP_N_SCORES,
        "Top_n_genes": TOP_N_GENES,
        "Feature_ranking": (
            "Original-family univariable OS BH q, then p, then name"
        ),
        "Feature_selection_requires_q_lt_05": False,
        "Feature_selection_requires_unpenalised_fit": True,
        "ESTIMATE_enabled": PU_ADD_ESTIMATE,
        "ESTIMATE_required": PU_REQUIRE_ESTIMATE,
        "ESTIMATE_edge_boot": PU_N_ESTIMATE_BOOT,
        "ESTIMATE_edge_seed": 260101,
        "External_GSE": EXTERNAL_GSE,
        "External_metadata_verified_flag": EXTERNAL_METADATA_VERIFIED,
        "External_log2_transform": EXTERNAL_LOG2_TRANSFORM,
    },
}, "reproducibility.json")

try:
    try:
        ipython = get_ipython()
    except NameError:
        ipython = None

    source = None
    source_origin = None

    if ipython is not None:
        history = ipython.history_manager.input_hist_raw
        count = ipython.execution_count
        candidate = history[count] if count < len(history) else history[-1]
        if (
            "S37_ESTIMATE_adjusted_edge_sensitivity.csv" in candidate
            and "MODEL_SELECTION_POLICY" in candidate
        ):
            source = candidate
            source_origin = "IPython executed cell"

    if source is None and "__file__" in globals():
        source_path = Path(__file__)
        if source_path.is_file():
            source = source_path.read_text(encoding="utf-8")
            source_origin = str(source_path)

    if source is not None:
        (OUT / "analysis_executed_with_selection.py").write_text(
            source, encoding="utf-8"
        )
        js({
            "Source_origin": source_origin,
            "Source_SHA256": hashlib.sha256(
                source.encode("utf-8")
            ).hexdigest(),
            "ESTIMATE_enabled": PU_ADD_ESTIMATE,
            "ESTIMATE_failure_stops_execution": PU_REQUIRE_ESTIMATE,
            "ABSOLUTE_minimum_matched_tumours": PU_MIN_MATCHED,
            "Outcome_based_feature_selection": True,
            "Association_based_edge_followup_selection": True,
            "Adjusted_Cox_CV_targets": MODEL_FEATURES,
            "Selected_edge_followup_targets": FOLLOWUP_EDGE_PAIRS,
            "CV_interpretation": CV_INTERPRETATION,
        }, "source_execution_provenance.json")
    else:
        NOTICES.append(
            "Automatic source export unavailable."
        )

except Exception as error:
    NOTICES.append("Source export unavailable: " + str(error))

js({"Notices": NOTICES}, "analysis_notices.json")

inventory = pd.DataFrame([
    {
        "File": str(path.relative_to(OUT)),
        "Bytes": path.stat().st_size,
    }
    for path in OUT.rglob("*") if path.is_file()
])

inventory.to_csv(OUT / "output_inventory.csv", index=False)
archive = shutil.make_archive(str(OUT), "zip", root_dir=OUT)


print("COMPLETE")
print("Initial main-score and gene features:", len(FOLLOWUP_FEATURES))
print("Adjusted Cox / CV scores:", TOP_SCORES)
print("Adjusted Cox / CV genes:", TOP_GENES)
print("All discovery gene pairs:", len(ALL_EDGE_PAIRS))
print("Rewired edge follow-up pairs:", len(FOLLOWUP_EDGE_PAIRS))
print("Outcome-informed feature selection was used.")
print("Association-informed edge follow-up selection was used.")
print(CV_INTERPRETATION)
print("Output directory:", OUT)
print("Output archive:", archive)
