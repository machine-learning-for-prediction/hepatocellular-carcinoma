order = DE.sort_values("Mean_difference_disjoint").Gene.tolist()
dd = DE.set_index("Gene").loc[order]
fig, axes = plt.subplots(1, 2, figsize=(14, 12), sharey=True)
positions = np.arange(len(order))

for ax, column, qcolumn, title in [
    (
        axes[0], "Mean_difference_disjoint", "q_disjoint",
        f"A: disjoint tumours n={len(dids)}; adjacent n={len(nids)}",
    ),
    (
        axes[1], "Mean_paired_difference", "q_paired",
        f"B: matched patient pairs n={len(paired)}",
    ),
]:
    ax.hlines(positions, 0, dd[column], color="lightgray")
    ax.scatter(
        dd[column], positions,
        c=np.where(dd[qcolumn] < .05, "#20639B", "#BBBBBB"),
        s=45,
    )
    ax.axvline(0, color="black", lw=.8)
    ax.set_title(title, loc="left", fontsize=11)
    ax.set_xlabel("Mean log-scale expression difference")

axes[0].set_yticks(positions)
axes[0].set_yticklabels(order, fontsize=11)
fig.suptitle(
    "Coloured points: BH q < 0.05 within the corresponding test family"
)
fig.tight_layout()
figsave(fig, "Fig01_DE")

fig, axes = plt.subplots(1, 3, figsize=(27, 9))

for ax, matrix, title, limit in [
    (axes[0], RN, f"Adjacent non-tumour n={len(nids)}", 1),
    (axes[1], RT, f"Patient-disjoint tumour n={len(dids)}", 1),
    (axes[2], D, "Delta rho: tumour minus adjacent", 2),
]:
    sns.heatmap(
        pd.DataFrame(matrix, index=genes, columns=genes),
        ax=ax, cmap="RdBu_r", center=0,
        vmin=-limit, vmax=limit, square=True,
        cbar_kws={"label": "rho" if limit == 1 else "delta rho"},
    )
    ax.set_title(title, loc="left")
    ax.tick_params(labelsize=10)

fig.tight_layout()
figsave(fig, "Fig02_correlations")

fig, axes = plt.subplots(1, 2, figsize=(17, 10))
ee = EDGES.head(20).iloc[::-1].reset_index(drop=True)

for i, row in ee.iterrows():
    color = "#C43C39" if row.Rewired else "#999999"
    axes[0].plot([row.CI_low, row.CI_high], [i, i], color=color)
    axes[0].scatter(row.delta_rho, i, color=color)

axes[0].set_yticks(range(len(ee)))
axes[0].set_yticklabels(ee.Gene1 + "–" + ee.Gene2)
axes[0].axvline(0, color="black", lw=.8)
axes[0].set_xlabel("Delta rho; pointwise bootstrap 95% CI")
axes[0].set_title("A: largest correlation changes", loc="left")

if not bs.empty:
    bb = bs.sort_values("Direction_frequency")
    pair_labels = bb.Gene1 + "–" + bb.Gene2
    axes[1].barh(
        pair_labels, bb.Direction_frequency,
        label="Same direction", alpha=.75,
    )
    axes[1].scatter(
        bb.BW_q_frequency, pair_labels,
        color="black", label="BH significance frequency",
    )
    axes[1].set_xlim(0, 1)
    axes[1].legend()
    axes[1].tick_params(
        axis="y",
        labelsize=max(4, min(10, 250 / len(bb))),
    )
else:
    axes[1].axis("off")
    axes[1].text(
        .5, .5, "No primary rewired pairs",
        ha="center", va="center",
    )

axes[1].set_title(
    "B: balanced-size stability",
    loc="left",
)
fig.tight_layout()
figsave(fig, "Fig03_edge_robustness")

forest = GENE_OS.dropna(subset=["HR"]).copy()
fig, ax = plt.subplots(figsize=(9, 11))

forest = forest.sort_values("HR", ascending=True)
for i, row in enumerate(forest.itertuples()):
    color = "#C43C39" if row.q_gene_family < .05 else "#999999"
    ax.plot([row.CI_low, row.CI_high], [i, i], color=color)
    ax.scatter(row.HR, i, color=color)

ax.set_yticks(range(len(forest)))
ax.set_yticklabels(forest.Feature)
ax.set_xscale("log")
ax.axvline(1, color="black", ls="--")
ax.set_xlabel("Univariable OS HR per SD; 95% CI")
ax.set_title("Colour: individual-gene-family BH q < 0.05", loc="left")
fig.tight_layout()
figsave(fig, "Fig04_gene_OS")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
sns.boxplot(
    data=SIM,
    x="Mean_shift",
    y="Median_default_evaluable",
    hue="Covariance",
    ax=axes[0],
)
axes[0].tick_params(axis="x", rotation=20)
axes[0].set_title("A: score responds to means and covariance")
axes[1].scatter(nsd, np.abs(ZN).mean(axis=1))

for i, gene in enumerate(genes):
    axes[1].annotate(
        gene, (nsd.iloc[i], np.abs(ZN).mean(axis=1)[i]), fontsize=7
    )

axes[1].set_xlabel("Adjacent-reference SD")
axes[1].set_ylabel("Mean absolute standardized contribution")
axes[1].set_title("B: reference-variance influence")
fig.tight_layout()
figsave(fig, "FigS_construct_and_reference")