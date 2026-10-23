from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
OUT = PROCESSED / "analysis"

OUT.mkdir(parents=True, exist_ok=True)

MISSING_SENTINEL = 5.397605e-79


def read_xpt(path):
    df = pd.read_sas(path, format="xport")

    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(
                lambda x: x.decode("utf-8").strip()
                if isinstance(x, bytes)
                else x
            )

    for col in df.select_dtypes(include=[np.number]).columns:
        df.loc[
            np.isclose(
                df[col],
                MISSING_SENTINEL,
                rtol=0,
                atol=1e-100,
            ),
            col,
        ] = np.nan

    return df


def keep_existing(df, columns):
    return df[
        [col for col in columns if col in df.columns]
    ].copy()


def load_cycle(cycle, suffix):

    base = RAW / cycle

    # --------------------------------------------------------
    # Core dataset
    # --------------------------------------------------------

    core = pd.read_csv(
        PROCESSED / "core_sleep_glycemic.csv"
    )

    core = core.loc[
        core["cycle"] == cycle
    ].copy()

    # --------------------------------------------------------
    # Demographics
    # --------------------------------------------------------

    demo = read_xpt(
        base / f"DEMO_{suffix}.XPT"
    )

    demo = keep_existing(
        demo,
        [
            "SEQN",
            "RIDAGEYR",
            "RIAGENDR",
            "RIDRETH3",
            "DMDEDUC2",
            "DMDMARTL",
            "INDFMPIR",
        ],
    )

    # --------------------------------------------------------
    # Body measurements
    # --------------------------------------------------------

    bmx = read_xpt(
        base / f"BMX_{suffix}.XPT"
    )

    bmx = keep_existing(
        bmx,
        [
            "SEQN",
            "BMXBMI",
            "BMXWAIST",
        ],
    )

    # --------------------------------------------------------
    # Blood pressure
    # --------------------------------------------------------

    bpx = read_xpt(
        base / f"BPX_{suffix}.XPT"
    )

    bpx = keep_existing(
        bpx,
        [
            "SEQN",
            "BPXSY1",
            "BPXDI1",
            "BPXSY2",
            "BPXDI2",
            "BPXSY3",
            "BPXDI3",
            "BPXSY4",
            "BPXDI4",
        ],
    )

    # --------------------------------------------------------
    # Physical activity
    # --------------------------------------------------------

    paq = read_xpt(
        base / f"PAQ_{suffix}.XPT"
    )

    paq = keep_existing(
        paq,
        [
            "SEQN",
            "PAQ605",
            "PAQ620",
            "PAQ635",
            "PAQ650",
            "PAQ665",
            "PAD680",
        ],
    )

    # --------------------------------------------------------
    # Alcohol
    # --------------------------------------------------------

    alq = read_xpt(
        base / f"ALQ_{suffix}.XPT"
    )

    alq = keep_existing(
        alq,
        [
            "SEQN",
            "ALQ101",
            "ALQ110",
            "ALQ120Q",
            "ALQ120U",
            "ALQ130",
        ],
    )

    # --------------------------------------------------------
    # Smoking
    # --------------------------------------------------------

    smq = read_xpt(
        base / f"SMQ_{suffix}.XPT"
    )

    smq = keep_existing(
        smq,
        [
            "SEQN",
            "SMQ020",
            "SMQ040",
            "SMQ050Q",
            "SMQ050U",
        ],
    )

    # --------------------------------------------------------
    # Depression
    # --------------------------------------------------------

    dpq = read_xpt(
        base / f"DPQ_{suffix}.XPT"
    )

    dpq_cols = [
        "SEQN",
        "DPQ010",
        "DPQ020",
        "DPQ030",
        "DPQ040",
        "DPQ050",
        "DPQ060",
        "DPQ070",
        "DPQ080",
        "DPQ090",
        "DPQ100",
    ]

    dpq = keep_existing(
        dpq,
        dpq_cols,
    )

    # --------------------------------------------------------
    # Diabetes
    # --------------------------------------------------------

    diq = read_xpt(
        base / f"DIQ_{suffix}.XPT"
    )

    diq = keep_existing(
        diq,
        [
            "SEQN",
            "DIQ010",
            "DID040",
            "DIQ160",
            "DIQ170",
            "DIQ180",
        ],
    )

    # --------------------------------------------------------
    # Merge
    # --------------------------------------------------------

    datasets = [
        demo,
        bmx,
        bpx,
        paq,
        alq,
        smq,
        dpq,
        diq,
    ]

    df = core.copy()

    for extra in datasets:
        if "SEQN" not in extra.columns:
            continue

        df = df.merge(
            extra,
            on="SEQN",
            how="left",
            validate="one_to_one",
        )

    # --------------------------------------------------------
    # Recompute circular sleep timing
    # --------------------------------------------------------

    df["midpoint_sin"] = np.sin(
        2 * np.pi *
        df["sleep_midpoint_min"] /
        1440
    )

    df["midpoint_cos"] = np.cos(
        2 * np.pi *
        df["sleep_midpoint_min"] /
        1440
    )

    # --------------------------------------------------------
    # Convert depression score
    # --------------------------------------------------------

    dpq_cols = [
        col for col in df.columns
        if col.startswith("DPQ")
        and col != "DPQ_TOTAL"
    ]

    if dpq_cols:
        df["DPQ_TOTAL"] = (
            df[dpq_cols]
            .replace(
                {
                    7: np.nan,
                    9: np.nan,
                }
            )
            .sum(
                axis=1,
                min_count=len(dpq_cols),
            )
        )

    # --------------------------------------------------------
    # Blood pressure summary
    # --------------------------------------------------------

    systolic_cols = [
        col for col in [
            "BPXSY1",
            "BPXSY2",
            "BPXSY3",
            "BPXSY4",
        ]
        if col in df.columns
    ]

    diastolic_cols = [
        col for col in [
            "BPXDI1",
            "BPXDI2",
            "BPXDI3",
            "BPXDI4",
        ]
        if col in df.columns
    ]

    if systolic_cols:
        df["SBP_mean"] = df[
            systolic_cols
        ].mean(axis=1)

    if diastolic_cols:
        df["DBP_mean"] = df[
            diastolic_cols
        ].mean(axis=1)

    # --------------------------------------------------------
    # Primary analysis flags
    # --------------------------------------------------------

    df["sleep_exposure_available"] = (
        df["valid_sleep_timing"] &
        df["valid_sleep_duration"]
    )

    df["hba1c_available"] = (
        df["LBXGH"].notna()
    )

    df["primary_candidate"] = (
        df["sleep_exposure_available"] &
        df["hba1c_available"]
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    output = (
        OUT /
        f"analysis_{cycle.replace('-', '_')}.csv"
    )

    df.to_csv(
        output,
        index=False,
    )

    return df, output


# ============================================================
# Build both cycles
# ============================================================

all_data = []

for cycle, suffix in [
    ("2015-2016", "I"),
    ("2017-2018", "J"),
]:

    print("\n" + "=" * 80)
    print(cycle)
    print("=" * 80)

    df, output = load_cycle(
        cycle,
        suffix,
    )

    print(f"Rows: {len(df):,}")
    print(f"Columns: {len(df.columns):,}")
    print(f"Saved: {output}")

    print("\nPrimary candidate:")
    print(
        df["primary_candidate"]
        .value_counts(dropna=False)
        .to_string()
    )

    all_data.append(df)


# ============================================================
# Combined analysis dataset
# ============================================================

combined = pd.concat(
    all_data,
    ignore_index=True,
)

assert combined["SEQN"].is_unique

combined_output = (
    OUT /
    "analysis_combined.csv"
)

combined.to_csv(
    combined_output,
    index=False,
)

print("\n" + "=" * 80)
print("COMBINED ANALYSIS DATASET")
print("=" * 80)

print(
    f"Rows: {len(combined):,}"
)

print(
    f"Columns: {len(combined.columns):,}"
)

print(
    f"Saved: {combined_output}"
)

# ------------------------------------------------------------
# Candidate-variable missingness
# ------------------------------------------------------------

candidate_columns = [
    "RIDAGEYR",
    "RIAGENDR",
    "RIDRETH3",
    "DMDEDUC2",
    "DMDMARTL",
    "INDFMPIR",
    "BMXBMI",
    "BMXWAIST",
    "SBP_mean",
    "DBP_mean",
    "PAD680",
    "ALQ130",
    "SMQ020",
    "DPQ_TOTAL",
    "DIQ010",
    "sleep_midpoint_min",
    "sleep_duration_derived",
    "SLD012",
    "LBXGH",
]

candidate_columns = [
    col for col in candidate_columns
    if col in combined.columns
]

missingness = (
    combined[candidate_columns]
    .isna()
    .mean()
    .mul(100)
    .sort_values(
        ascending=False
    )
    .round(2)
)

print("\nCandidate variable missingness (%):")
print(
    missingness.to_string()
)

missingness.to_csv(
    OUT / "candidate_missingness.csv"
)

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
