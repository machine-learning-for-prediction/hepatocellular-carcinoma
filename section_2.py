def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def js(value, name):
    (META / name).write_text(
        json.dumps(clean(value), indent=2, allow_nan=False),
        encoding="utf-8",
    )


def tab(df, name):
    df.to_csv(TABLES / name, index=False)


def figsave(fig, name):
    fig.savefig(FIGURES / (name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / (name + ".pdf"), bbox_inches="tight")
    if SHOW_FIGURES:
        plt.show()
    plt.close(fig)


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    values = p[ok]
    order = np.argsort(values)
    adjusted = np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    out[ok] = restored
    return out


def pid(value):
    return str(value)[:12]


def stype(value):
    return str(value)[13:15]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pu_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


session = requests.Session()
session.mount(
    "https://",
    HTTPAdapter(max_retries=Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
    )),
)


def download(url, path):
    path = Path(path)
    if not path.exists():
        with session.get(url, stream=True, timeout=(30, 300)) as response:
            response.raise_for_status()
            temporary = Path(str(path) + ".part")
            with open(temporary, "wb") as handle:
                for block in response.iter_content(1024 * 1024):
                    if block:
                        handle.write(block)
            temporary.replace(path)
        Path(str(path) + ".retrieved.txt").write_text(
            datetime.now(timezone.utc).isoformat()
        )

    retrieved = Path(str(path) + ".retrieved.txt")
    INPUTS.append({
        "url": url,
        "file": str(path),
        "sha256": sha(path),
        "retrieved_utc": (
            retrieved.read_text()
            if retrieved.exists()
            else "Unknown: cached before retrieval logging"
        ),
    })
    return path


def corr(values):
    return np.corrcoef(stats.rankdata(np.asarray(values, float), axis=1))


def sp(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if (
        ok.sum() < 5
        or np.std(a[ok]) < EPS
        or np.std(b[ok]) < EPS
    ):
        return np.nan, np.nan, int(ok.sum())
    rho, p = stats.spearmanr(a[ok], b[ok])
    return float(rho), float(p), int(ok.sum())


def ci(values, minimum_valid=20):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) < minimum_valid:
        return [np.nan, np.nan]
    return np.percentile(values, [2.5, 97.5]).tolist()


def bw(rt, rn, nt, nn):
    rt = np.clip(rt, -.999999, .999999)
    rn = np.clip(rn, -.999999, .999999)
    variance = (
        (1 + rt**2 / 2) / (nt - 3)
        + (1 + rn**2 / 2) / (nn - 3)
    )
    return 2 * stats.norm.sf(
        np.abs((np.arctanh(rt) - np.arctanh(rn)) / np.sqrt(variance))
    )


def upperp(observed, null):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or not len(null):
        return np.nan
    return (1 + np.sum(null >= observed)) / (len(null) + 1)


def cindex(t, e, risk):
    try:
        return float(concordance_index(t, -np.asarray(risk), e))
    except Exception:
        return np.nan


def comparable(t, e):
    t = np.asarray(t)
    e = np.asarray(e)
    i, j = np.triu_indices(len(t), 1)
    return int((
        ((t[i] < t[j]) & (e[i] == 1))
        | ((t[j] < t[i]) & (e[j] == 1))
        | ((t[i] == t[j]) & (e[i] != e[j]))
    ).sum())

if not EXT_META_FILE.exists():
    print(
        "\nWARNING: external_metadata_reviewed.csv is not present.\n"
        "External inference will not be performed without it.\n"
        "The main TCGA, ABSOLUTE and ESTIMATE analyses do not require it.\n"
    )