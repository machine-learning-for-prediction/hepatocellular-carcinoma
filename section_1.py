MODULES = {
    "Core synthesis": ["ODC1", "AMD1", "SRM", "SMS"],
    "Antizyme inhibitors": ["AZIN1", "AZIN2"],
    "Antizymes": ["OAZ1", "OAZ2", "OAZ3"],
    "Catabolism": ["SAT1", "PAOX", "SMOX"],
    "Arginine/ornithine supply": ["ARG1", "ARG2", "ASS1", "ASL"],
    "NOS / arginine competition": ["NOS1", "NOS2", "NOS3"],
    "Methionine/SAM/MTA support": [
        "MAT1A", "MAT2A", "MTAP", "APIP", "ENOPH1"
    ],
    "Transport-associated": [
        "SLC3A2", "SLC7A1", "ATP13A2", "ATP13A3", "ATP13A4", "ATP13A5"
    ],
}

GENE_MODULE = {gene: module for module, gl in MODULES.items() for gene in gl}
CURATED = sorted(GENE_MODULE)

SUPPORT = sorted(set(
    MODULES["Core synthesis"]
    + MODULES["Antizyme inhibitors"]
    + MODULES["Arginine/ornithine supply"]
    + MODULES["Methionine/SAM/MTA support"]
))

BRAKE = MODULES["Antizymes"] + MODULES["Catabolism"]

MAIN_SCORES = [
    "Core Synthesis",
    "Catabolism",
    "Transport-associated",
    "Synthesis Support",
    "Homeostatic Brake",
    "Balance Index",
    "Dysregulation Burden",
    "Activation Burden",
    "Signed Burden",
    "Coordination Deviation",
]

ALT_SCORES = ["Coordination weighted", "Coordination continuous"]

RELATIONS = [
    ("ARG1", "ODC1", "substrate supply"),
    ("ARG2", "ODC1", "substrate supply"),
    ("ASS1", "ASL", "consecutive reactions"),
    ("ASL", "ARG1", "substrate supply"),
    ("ASL", "ARG2", "substrate supply"),
    ("ODC1", "SRM", "substrate supply"),
    ("SRM", "SMS", "substrate supply"),
    ("AMD1", "SRM", "donor supply"),
    ("AMD1", "SMS", "donor supply"),
    ("MAT1A", "AMD1", "substrate supply"),
    ("MAT2A", "AMD1", "substrate supply"),
    ("SRM", "MTAP", "by-product salvage"),
    ("SMS", "MTAP", "by-product salvage"),
    ("MTAP", "APIP", "indirect pathway connection"),
    ("APIP", "ENOPH1", "consecutive reactions"),
    ("OAZ1", "ODC1", "regulatory annotation"),
    ("OAZ2", "ODC1", "regulatory annotation"),
    ("OAZ3", "ODC1", "regulatory annotation"),
    ("AZIN1", "OAZ1", "regulatory annotation"),
    ("AZIN1", "OAZ2", "regulatory annotation"),
    ("AZIN1", "OAZ3", "regulatory annotation"),
    ("AZIN1", "ODC1", "indirect regulation"),
    ("SAT1", "PAOX", "substrate supply"),
    ("SAT1", "SRM", "product handling"),
    ("SAT1", "SMS", "product handling"),
    ("SMOX", "SMS", "product handling"),
    ("SAT1", "SLC3A2", "transport-associated annotation"),
    ("NOS1", "ARG1", "substrate competition"),
    ("NOS2", "ARG1", "substrate competition"),
    ("NOS3", "ARG1", "substrate competition"),
    ("NOS1", "ARG2", "substrate competition"),
    ("NOS2", "ARG2", "substrate competition"),
    ("NOS3", "ARG2", "substrate competition"),
    ("SLC7A1", "ARG1", "substrate uptake"),
    ("SLC7A1", "ARG2", "substrate uptake"),
]

EC = {
    "ARG1": "3.5.3.1",
    "ARG2": "3.5.3.1",
    "ODC1": "4.1.1.17",
    "ASS1": "6.3.4.5",
    "ASL": "4.3.2.1",
    "SRM": "2.5.1.16",
    "SMS": "2.5.1.22",
    "AMD1": "4.1.1.50",
    "MAT1A": "2.5.1.6",
    "MAT2A": "2.5.1.6",
    "SAT1": "2.3.1.57",
    "PAOX": "1.5.3.13",
    "SMOX": "1.5.3.16",
}

NARROW_PAIRS = [
    ("ARG1", "ODC1"),
    ("ARG2", "ODC1"),
    ("ASS1", "ASL"),
    ("ASL", "ARG1"),
    ("ASL", "ARG2"),
    ("ODC1", "SRM"),
    ("SRM", "SMS"),
    ("AMD1", "SRM"),
    ("AMD1", "SMS"),
    ("MAT1A", "AMD1"),
    ("MAT2A", "AMD1"),
    ("SAT1", "PAOX"),
]

PROXY_MARKERS = {
    "Immune proxy": ["PTPRC", "CD3D", "CD3E", "CD247", "LST1", "TYROBP"],
    "Stromal proxy": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "Hepatocyte differentiation proxy": [
        "ALB", "TTR", "APOA1", "APOA2", "CYP3A4", "CYP2E1", "HNF4A"
    ],
    "Proliferation proxy": ["MKI67", "PCNA", "MCM2", "MCM4", "TOP2A"],
}
