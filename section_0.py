import sys
import subprocess
import importlib.util

PACKAGES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "sklearn": "scikit-learn",
    "lifelines": "lifelines",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "requests": "requests",
    "tqdm": "tqdm",
    "GEOparse": "GEOparse",
}

for module, package in PACKAGES.items():
    if importlib.util.find_spec(module) is None:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", package]
        )

import os
import re
import json
import hashlib
import warnings
import shutil
import platform
from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index, ConvergenceWarning
import matplotlib.pyplot as plt
import seaborn as sns
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm.auto import tqdm
import GEOparse

if Path("/content").exists():
    os.chdir("/content")

SEED = 42
CV_SEEDS = [42, 43, 44, 45, 46]
N_FOLDS = 5

N_EDGE_BOOT = 5000
N_NETWORK_PERM = 5000
N_GRAPH_PERM = 10000
N_REFERENCE_BOOT = 200
N_BALANCED_REPEATS = 200
N_BALANCED_BOOT = 200
N_ADJUSTED_BOOT = 500
N_SIM_REPEATS = 200

COORD_Z = 1.0
COORD_RHO = 0.20
COORD_PAIRS = 5
EPS = 1e-12
RIDGES = [0.01, 0.05, 0.10]
SHOW_FIGURES = True

TOP_N_SCORES = 6
TOP_N_GENES = 6

DOWNLOAD_EXTERNAL = True
EXTERNAL_GSE = "GSE76427"
EXTERNAL_METADATA_VERIFIED = True
EXTERNAL_LOG2_TRANSFORM = False

PU_ADD_ESTIMATE = True
PU_REQUIRE_ESTIMATE = True
PU_MIN_MATCHED = 30
PU_N_ESTIMATE_BOOT = 500

BASE = Path("/content") if Path("/content").exists() else Path.cwd()
OUT = BASE / (
    "polyamine_revision_"
    + datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
)
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
META = OUT / "metadata"
CACHE = BASE / "polyamine_revision_cache"

for path in [OUT, TABLES, FIGURES, META, CACHE]:
    path.mkdir(parents=True, exist_ok=True)

EXT_META_FILE = BASE / "external_metadata_reviewed.csv"

PU_ABSOLUTE_LOCAL = BASE / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
PU_ABSOLUTE_UUID = "4f277128-f793-4354-a13d-30cc7fe9f6b5"
PU_ABSOLUTE_URL = "https://api.gdc.cancer.gov/data/" + PU_ABSOLUTE_UUID
PU_ABSOLUTE_FILENAME = "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
PU_ABSOLUTE_EXPECTED_MD5 = "8ea2ca92c8ae58350538999dfa1174da"
PU_ABSOLUTE_EXPECTED_BYTES = 901812

HUB = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
EXPR_ID = "TCGA.LIHC.sampleMap/HiSeqV2"
SURV_ID = "survival/LIHC_survival.txt"
CLIN_ID = "TCGA.LIHC.sampleMap/LIHC_clinicalMatrix"

EXPR_URL = HUB + "TCGA.LIHC.sampleMap%2FHiSeqV2.gz"
SURV_URL = HUB + "survival%2FLIHC_survival.txt"
CLIN_URL = HUB + "TCGA.LIHC.sampleMap%2FLIHC_clinicalMatrix"

LOG, NOTICES, INPUTS, CV_LOG = [], [], [], []
PU_STATUS = []
PU_ESTIMATE = None
PU_ESTIMATE_AVAILABLE = False

MODEL_FEATURES = []
TOP_SCORES = []
TOP_GENES = []
FOLLOWUP_EDGE_PAIRS = []

MODEL_SELECTION_POLICY = (
    "Outcome-informed discovery shortlist: up to 6 main scores and 6 genes ranked by original-family univariable OS BH q, then p-value, then feature name; successful unpenalised fits only"
)
EDGE_FOLLOWUP_POLICY = (
    "Primary TCGA rewired pairs only: BW BH q < 0.05, at least "
    "90% valid bootstraps, and pointwise bootstrap CI excluding zero"
)
CV_INTERPRETATION = (
    "Post-selection descriptive CV: feature identities were selected "
    "using full-cohort OS outcomes before CV; selection leakage may cause optimistic performance estimates"
)

sns.set_theme(style="whitegrid")
plt.rcParams["pdf.fonttype"] = 42