def network_stats(delta):
    absolute = np.abs(delta).copy()
    np.fill_diagonal(absolute, 0)
    values = absolute[IU]
    return (
        values.mean(),
        absolute.sum(axis=1),
        np.array([values[index].mean() for index in bix]),
    )


obs_global, obs_node, obs_block = network_stats(D)


def permutation_network(a, b, uniform, nperm, seed):
    if uniform:
        a = (stats.rankdata(a, axis=1) - .5) / a.shape[1]
        b = (stats.rankdata(b, axis=1) - .5) / b.shape[1]

    pool = np.hstack([b, a])
    nn = b.shape[1]
    rng = np.random.default_rng(seed)
    global_stats, node_stats, block_stats = [], [], []
    invalid = 0
    attempts = 0

    while len(global_stats) < nperm:
        attempts += 1
        if attempts > 20 * nperm:
            raise RuntimeError("Excess undefined permutation matrices.")

        index = rng.permutation(pool.shape[1])
        with np.errstate(invalid="ignore", divide="ignore"):
            delta = (
                corr(pool[:, index[nn:]])
                - corr(pool[:, index[:nn]])
            )

        if not np.isfinite(delta).all():
            invalid += 1
            continue

        global_value, node_value, block_value = network_stats(delta)
        global_stats.append(global_value)
        node_stats.append(node_value)
        block_stats.append(block_value)

    return (
        np.array(global_stats),
        np.array(node_stats),
        np.array(block_stats),
        invalid,
    )


perm_rows, local_rows = [], []

for uniform, label, seed in [
    (True, "Rank-uniformised approximate", 43),
    (False, "Conventional raw-expression sensitivity", 44),
]:
    gs, ns, block_null, invalid = permutation_network(
        TD.to_numpy(), N.to_numpy(), uniform, N_NETWORK_PERM, seed
    )

    perm_rows.append({
        "Method": label,
        "Observed_mean_abs_delta": obs_global,
        "Null_mean": gs.mean(),
        "p": upperp(obs_global, gs),
        "N_valid": len(gs),
        "Invalid_attempts": invalid,
    })

    npv = np.array([
        upperp(obs_node[i], ns[:, i]) for i in range(G)
    ])
    bpv = np.array([
        upperp(obs_block[i], block_null[:, i])
        for i in range(len(block_names))
    ])

    for names, values, pvalues, level in [
        (genes, obs_node, npv, "Node"),
        (block_names, obs_block, bpv, "Module block"),
    ]:
        for name, value, pvalue, qvalue in zip(
            names, values, pvalues, bh(pvalues)
        ):
            local_rows.append({
                "Method": label,
                "Level": level,
                "Name": name,
                "Statistic": value,
                "p": pvalue,
                "q_within_level": qvalue,
            })

tab(pd.DataFrame(perm_rows), "S06_global_permutation_methods.csv")
tab(pd.DataFrame(local_rows), "S07_node_module_tests.csv")