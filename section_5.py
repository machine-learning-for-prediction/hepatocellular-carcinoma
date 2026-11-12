broad = {frozenset((a, b)) for a, b, _ in RELATIONS}
narrow = {frozenset(pair) for pair in NARROW_PAIRS}

relationship_rows = []

for a, b, kind in RELATIONS:
    is_narrow = frozenset((a, b)) in narrow
    links = (
        ";".join(
            "https://www.kegg.jp/entry/ec:" + EC[gene]
            for gene in [a, b] if gene in EC
        )
        if is_narrow else ""
    )

    relationship_rows.append({
        "Gene1": a,
        "Gene2": b,
        "Relationship_type": kind,
        "Available": a in genes and b in genes,
        "Narrow_reaction_linked": is_narrow,
        "Source_pointer": links,
        "Verification_status": (
            "Reaction-level EC pointers"
            if is_narrow
            else "refer to S00 for sources"
        ),
        "Interpretation": (
            "Shared metabolite/reaction connection; not physical interaction"
            if is_narrow else "Curated broader association"
        ),
    })

tab(pd.DataFrame(relationship_rows), "S01_relationship_provenance.csv")

rng = np.random.default_rng(45)
graph_perms = np.array([
    rng.permutation(G) for _ in range(N_GRAPH_PERM)
])

graph_results = []
absD = np.abs(D)
allvals = absD[IU]
total = allvals.sum()

for name, relation_set in [
    ("Broad curated", broad),
    ("Narrow reaction-linked", narrow),
]:
    mask = np.array([
        frozenset((genes[i], genes[j])) in relation_set
        for i, j in zip(*IU)
    ])

    if not mask.any() or mask.all():
        continue

    ii, jj = IU[0][mask], IU[1][mask]
    nr = int(mask.sum())
    no = len(mask) - nr
    observed = allvals[mask].mean() - allvals[~mask].mean()

    null = []
    for permutation in graph_perms:
        value = absD[permutation[ii], permutation[jj]].sum()
        null.append(value / nr - (total - value) / no)

    graph_results.append({
        "Set": name,
        "N_related": nr,
        "N_other": no,
        "Mean_related": allvals[mask].mean(),
        "Mean_other": allvals[~mask].mean(),
        "Difference": observed,
        "Graph_label_p": upperp(observed, null),
    })

gr = pd.DataFrame(graph_results)
if len(gr):
    gr["q_across_two_sets"] = bh(gr.Graph_label_p)
tab(gr, "S12_relationship_sensitivity.csv")

module_cons = []

for cohort, matrix in [("Adjacent", N), ("Disjoint tumour", TD)]:
    for module, gl in all_modules.items():
        gg = [gene for gene in gl if gene in matrix.index]
        if not gg:
            continue

        values = matrix.loc[gg].to_numpy()
        count = len(gg)

        if count >= 2:
            correlations = corr(values)
            upper = correlations[np.triu_indices(count, 1)]
            standardized = (
                values - values.mean(axis=1, keepdims=True)
            ) / values.std(axis=1, ddof=1, keepdims=True)

            pc = (
                PCA().fit(standardized.T).explained_variance_ratio_[0]
                if np.isfinite(standardized).all() else np.nan
            )
            mean_r = float(upper.mean())
            mean_abs = float(np.abs(upper).mean())
        else:
            pc = mean_r = mean_abs = np.nan

        module_cons.append({
            "Cohort": cohort,
            "Module": module,
            "N_genes": count,
            "Mean_signed_Spearman": mean_r,
            "Mean_absolute_Spearman": mean_abs,
            "PC1_variance_fraction": pc,
            "Note": (
                "Single-gene module: consistency not estimable"
                if count == 1
                else "Functional annotation; consistency is descriptive"
            ),
        })

tab(pd.DataFrame(module_cons), "S13_module_consistency.csv")