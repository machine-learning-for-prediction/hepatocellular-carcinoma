xp = download(EXPR_URL, CACHE / "HiSeqV2.gz")
yp = download(SURV_URL, CACHE / "LIHC_survival.txt")

raw = pd.read_csv(xp, sep="\t", index_col=0)
raw.index = raw.index.astype(str)
duplicate_gene_rows = int(raw.index.duplicated().sum())
raw = raw.loc[~raw.index.duplicated(keep="first")].apply(
    pd.to_numeric, errors="coerce"
)

bm = pd.DataFrame([
    {
        "position": i,
        "original": str(sample),
        "sample": str(sample)[:15],
        "patient": pid(sample),
        "type": stype(sample),
        "vial": str(sample)[15:16] or "~",
    }
    for i, sample in enumerate(raw.columns)
])

selected = bm.sort_values(
    ["sample", "vial", "original", "position"], kind="stable"
).drop_duplicates("sample")

X = raw.iloc[:, selected.position.to_numpy()].copy()
X.columns = selected["sample"].tolist()
tab(bm, "barcode_mapping.csv")

sv = pd.read_csv(yp, sep="\t")
if not {"sample", "OS.time", "OS"}.issubset(sv.columns):
    raise RuntimeError("Unexpected survival columns.")

sv["sample15"] = sv["sample"].astype(str).str[:15]
sv["OS.time"] = pd.to_numeric(sv["OS.time"], errors="coerce")
sv["OS"] = pd.to_numeric(sv["OS"], errors="coerce")

for sample, subset in sv.groupby("sample15"):
    if len(subset[["OS.time", "OS"]].drop_duplicates()) > 1:
        raise RuntimeError("Conflicting survival records: " + sample)

sv = sv.drop_duplicates("sample15").set_index("sample15")

tumour_all = [sample for sample in X if stype(sample) == "01"]
normal_all = [sample for sample in X if stype(sample) == "11"]
common = sorted(set(tumour_all) & set(sv.index))

valid = [
    sample for sample in common
    if pd.notna(sv.loc[sample, "OS.time"])
    and sv.loc[sample, "OS.time"] > 0
    and sv.loc[sample, "OS"] in [0, 1]
]


def dedup(ids):
    frame = selected.set_index("sample").loc[ids].reset_index()
    return frame.sort_values(
        ["patient", "vial", "original"], kind="stable"
    ).drop_duplicates("patient")["sample"].tolist()


tids = dedup(valid)
nids = dedup(normal_all)

if len(tids) < 10 or len(nids) < 5:
    raise RuntimeError("Insufficient selected tumour or adjacent samples.")

Xt = X[tids]
Xn = X[nids]

keep = (Xt.var(axis=1) > 0) & ((Xt > 0).mean(axis=1) >= .10)
filtered = X.loc[keep]
genes = [gene for gene in CURATED if gene in filtered.index]

if len(genes) < 6:
    raise RuntimeError("Fewer than six pathway genes remain.")

T = filtered.loc[genes, tids]
N = filtered.loc[genes, nids]

if (
    not np.isfinite(T.to_numpy()).all()
    or not np.isfinite(N.to_numpy()).all()
):
    raise RuntimeError("Nonfinite retained pathway expression.")

if (N.std(axis=1) <= EPS).any():
    raise RuntimeError("Zero adjacent-liver SD in retained pathway genes.")

tb = {pid(sample): sample for sample in tids}
nb = {pid(sample): sample for sample in nids}
paired = sorted(set(tb) & set(nb))
pt = [tb[patient] for patient in paired]
pn = [nb[patient] for patient in paired]
dids = [sample for sample in tids if pid(sample) not in nb]
TD = T[dids]

if len(dids) < 5:
    raise RuntimeError("Insufficient patient-disjoint tumour samples.")

Y = pd.DataFrame({
    "T": sv.loc[tids, "OS.time"].to_numpy(float),
    "E": sv.loc[tids, "OS"].to_numpy(int),
}, index=tids)

G = len(genes)
IU = np.triu_indices(G, 1)
FOLLOWUP_FEATURES = MAIN_SCORES + genes
ALL_EDGE_PAIRS = [(genes[i], genes[j]) for i, j in zip(*IU)]

SCORE_GENE_REQUIREMENTS = {
    "Core Synthesis": [
        gene for gene in MODULES["Core synthesis"] if gene in genes
    ],
    "Catabolism": [
        gene for gene in MODULES["Catabolism"] if gene in genes
    ],
    "Transport-associated": [
        gene for gene in MODULES["Transport-associated"] if gene in genes
    ],
    "Synthesis Support": [gene for gene in SUPPORT if gene in genes],
    "Homeostatic Brake": [gene for gene in BRAKE if gene in genes],
    "Balance Index": sorted((set(SUPPORT) | set(BRAKE)) & set(genes)),
    "Dysregulation Burden": list(genes),
    "Activation Burden": list(genes),
    "Signed Burden": list(genes),
    "Coordination Deviation": list(genes),
}


def feature_values(feature):
    return T.loc[feature] if feature in genes else SCORES[feature]


tab(
    pd.DataFrame(ALL_EDGE_PAIRS, columns=["Gene1", "Gene2"]),
    "Table04_all_discovery_edge_pairs.csv",
)

flow = pd.DataFrame([
    ("Unique sample15 expression columns", X.shape[1]),
    ("Primary tumours", len(tumour_all)),
    ("Adjacent non-tumour samples", len(normal_all)),
    ("Primary tumours with OS record", len(common)),
    ("After valid OS filtering", len(valid)),
    ("Tumours after patient deduplication", len(tids)),
    ("Adjacent samples after patient deduplication", len(nids)),
    ("Matched patient pairs", len(paired)),
    ("Patient-disjoint tumours", len(dids)),
    ("Genes after duplicate-row handling", len(raw)),
    ("Genes after tumour variance/detection filter", len(filtered)),
    ("Curated genes", len(CURATED)),
    ("Available curated genes", len(genes)),
    ("OS events", int(Y.E.sum())),
], columns=["Step", "Count"])

tab(flow, "Table01_cohort_flow.csv")
print(flow.to_string(index=False))


all_modules = dict(MODULES)
all_modules["Synthesis Support"] = SUPPORT
all_modules["Homeostatic Brake"] = BRAKE

tab(pd.DataFrame([
    {
        "Module": module,
        "Curated_genes": ";".join(gl),
        "Available_genes": ";".join(gene for gene in gl if gene in genes),
        "N_available": sum(gene in genes for gene in gl),
    }
    for module, gl in all_modules.items()
]), "S02_modules.csv")

js({
    "quantiles": np.nanpercentile(
        X.to_numpy(), [0, 1, 25, 50, 75, 99, 100]
    ),
    "note": "Distribution compatibility check only.",
}, "expression_scale.json")

de = []

for gene in genes:
    a = TD.loc[gene].to_numpy()
    b = N.loc[gene].to_numpy()

    u, p = stats.mannwhitneyu(
        a, b, alternative="two-sided", method="asymptotic"
    )
    _, pa = stats.mannwhitneyu(
        T.loc[gene], b, alternative="two-sided", method="asymptotic"
    )

    difference = (
        T.loc[gene, pt].to_numpy()
        - N.loc[gene, pn].to_numpy()
    )
    pp = (
        1.0 if len(difference) and np.all(difference == 0)
        else float(stats.wilcoxon(difference).pvalue)
        if len(difference) else np.nan
    )

    de.append({
        "Gene": gene,
        "Mean_difference_disjoint": a.mean() - b.mean(),
        "Cliffs_delta_disjoint": 2 * u / (len(a) * len(b)) - 1,
        "p_disjoint": p,
        "Mean_difference_all_descriptive": T.loc[gene].mean() - b.mean(),
        "p_all_exploratory": pa,
        "N_pairs": len(difference),
        "Mean_paired_difference": (
            difference.mean() if len(difference) else np.nan
        ),
        "p_paired": pp,
    })

DE = pd.DataFrame(de)

for family in ["disjoint", "all_exploratory", "paired"]:
    DE["q_" + family] = bh(DE["p_" + family])

tab(DE, "S03_S04_differential_expression.csv")

RN = corr(N)
RT = corr(TD)
RA = corr(T)
D = RT - RN

if not all(np.isfinite(matrix).all() for matrix in [RN, RT, RA]):
    raise RuntimeError("Undefined primary correlation matrix.")


def edge_intervals(a, b, nboot, seed, pair_indices=None):
    pair_indices = IU if pair_indices is None else pair_indices
    number_pairs = len(pair_indices[0])

    if number_pairs == 0:
        return np.array([]), np.array([]), np.array([], dtype=int)

    rng = np.random.default_rng(seed)
    out = np.full((nboot, number_pairs), np.nan)

    for repeat in range(nboot):
        aa = a[:, rng.integers(0, a.shape[1], a.shape[1])]
        bb = b[:, rng.integers(0, b.shape[1], b.shape[1])]
        with np.errstate(invalid="ignore", divide="ignore"):
            out[repeat] = (corr(aa) - corr(bb))[pair_indices]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        low, high = np.nanpercentile(out, [2.5, 97.5], axis=0)

    return low, high, np.isfinite(out).sum(axis=0)


lo, hi, nvalid = edge_intervals(
    TD.to_numpy(), N.to_numpy(), N_EDGE_BOOT, SEED
)

ep = bw(RT[IU], RN[IU], TD.shape[1], N.shape[1])
eq = bh(ep)


def pattern(rn, rt):
    if abs(rn) >= .30 and abs(rt) < .15:
        return "loss of coupling"
    if abs(rn) < .15 and abs(rt) >= .30:
        return "gain of coupling"
    if np.sign(rn) != np.sign(rt):
        return "sign reversal"
    return "attenuation" if abs(rt) < abs(rn) else "strengthening"


EDGES = pd.DataFrame([
    {
        "Gene1": genes[i],
        "Gene2": genes[j],
        "rho_adjacent": RN[i, j],
        "rho_tumour_disjoint": RT[i, j],
        "delta_rho": D[i, j],
        "abs_delta_rho": abs(D[i, j]),
        "BW_p": ep[k],
        "BW_q": eq[k],
        "CI_low": lo[k],
        "CI_high": hi[k],
        "Bootstrap_valid": nvalid[k],
        "Rewired": bool(
            eq[k] < .05
            and nvalid[k] >= .9 * N_EDGE_BOOT
            and (lo[k] > 0 or hi[k] < 0)
        ),
        "Pattern": pattern(RN[i, j], RT[i, j]),
    }
    for k, (i, j) in enumerate(zip(*IU))
]).sort_values("abs_delta_rho", ascending=False)

tab(EDGES, "S05_all_edges.csv")
tab(EDGES.loc[EDGES.Rewired], "Table02_rewired_edges.csv")

rewired_keys = {
    frozenset((row.Gene1, row.Gene2))
    for row in EDGES.loc[EDGES.Rewired].itertuples()
}

FOLLOWUP_EDGE_PAIRS = [
    pair for pair in ALL_EDGE_PAIRS
    if frozenset(pair) in rewired_keys
]

selected_iu = (
    np.array(
        [genes.index(a) for a, b in FOLLOWUP_EDGE_PAIRS], dtype=int
    ),
    np.array(
        [genes.index(b) for a, b in FOLLOWUP_EDGE_PAIRS], dtype=int
    ),
)

selected_edge_table = pd.DataFrame(
    FOLLOWUP_EDGE_PAIRS, columns=["Gene1", "Gene2"]
)
selected_edge_table["Target_policy"] = EDGE_FOLLOWUP_POLICY
tab(selected_edge_table, "Table05_selected_edge_followup_pairs.csv")

print("\nRewired pairs selected for follow-up:", len(FOLLOWUP_EDGE_PAIRS))

if not FOLLOWUP_EDGE_PAIRS:
    NOTICES.append(
        "No primary rewired pairs; selected edge follow-up analyses skipped."
    )

edge_key = {
    frozenset((genes[i], genes[j])): k
    for k, (i, j) in enumerate(zip(*IU))
}

mods = np.array([GENE_MODULE[gene] for gene in genes])
labels = np.array([
    " | ".join(sorted([mods[i], mods[j]]))
    for i, j in zip(*IU)
])
block_names = sorted(set(labels))
bix = [np.flatnonzero(labels == block) for block in block_names]