from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"

OUT.mkdir(parents=True, exist_ok=True)


MISSING_SENTINEL = 5.397605e-79


def read_xpt(path):
    df = pd.read_sas(path, format="xport")

    # Decode byte-string columns
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(
                lambda x: x.decode("utf-8").strip()
                if isinstance(x, bytes)
                else x
            )

    # Convert SAS missing sentinel to NaN
    for col in df.select_dtypes(include=[np.number]).columns:
        df.loc[
            np.isclose(
                df[col],
                MISSING_SENTINEL,
                rtol=0,
                atol=1e-100
            ),
            col
        ] = np.nan

    return df


def parse_time(value):
    if pd.isna(value):
        return np.nan

    value = str(value).strip()

    if not value or value.lower() == "nan":
        return np.nan

    try:
        hour, minute = value.split(":")
        hour = int(hour)
        minute = int(minute)

        if hour == 24:
            hour = 0

        if not (0 <= hour <= 23):
            return np.nan

        if not (0 <= minute <= 59):
            return np.nan

        return hour * 60 + minute

    except Exception:
        return np.nan


def circular_midpoint(sleep_min, wake_min):
    """
    Calculate midpoint of a sleep interval crossing midnight.

    Example:
        22:00 -> 06:00
        midpoint = 02:00
    """

    if pd.isna(sleep_min) or pd.isna(wake_min):
        return np.nan

    wake_adjusted = wake_min

    if wake_adjusted <= sleep_min:
        wake_adjusted += 1440

    midpoint = (sleep_min + wake_adjusted) / 2

    return midpoint % 1440


def sleep_interval_hours(sleep_min, wake_min):

    if pd.isna(sleep_min) or pd.isna(wake_min):
        return np.nan

    wake_adjusted = wake_min

    if wake_adjusted <= sleep_min:
        wake_adjusted += 1440

    duration = (wake_adjusted - sleep_min) / 60

    return duration


def load_cycle(cycle, suffix):

    base = RAW / cycle

    demo = read_xpt(base / f"DEMO_{suffix}.XPT")
    slq = read_xpt(base / f"SLQ_{suffix}.XPT")
    ghb = read_xpt(base / f"GHB_{suffix}.XPT")

    # Primary glucose outcome
    glu = read_xpt(base / f"GLU_{suffix}.XPT")

    # Fasting examination information
    fast = read_xpt(base / f"FASTQX_{suffix}.XPT")

    # Insulin is secondary for now
    ins = read_xpt(base / f"INS_{suffix}.XPT")

    # Keep only variables needed at this stage

    demo = demo[
        [
            "SEQN",
            "RIDAGEYR",
            "RIDSTATR",
            "RIAGENDR",
            "RIDRETH3",
            "DMDEDUC2",
            "DMDMARTL",
            "INDFMPIR",
            "WTINT2YR",
            "WTMEC2YR",
            "SDMVPSU",
            "SDMVSTRA",
        ]
    ]

    slq_cols = [
        "SEQN",
        "SLQ300",
        "SLQ310",
        "SLD012",
    ]

    # Weekend variables exist only in 2017-2018
    for col in ["SLQ320", "SLQ330", "SLD013"]:
        if col in slq.columns:
            slq_cols.append(col)

    slq = slq[slq_cols]

    ghb = ghb[
        [
            "SEQN",
            "LBXGH",
        ]
    ]

    glu = glu[
        [
            "SEQN",
            "WTSAF2YR",
            "LBXGLU",
        ]
    ]

    fast = fast[
        [
            "SEQN",
            "PHAFSTHR",
            "PHAFSTMN",
            "PHDSESN",
        ]
    ]

    ins_cols = [
        "SEQN",
        "WTSAF2YR",
        "LBXIN",
    ]

    ins = ins[ins_cols]

    # Merge one-to-one on SEQN
    df = demo.merge(
        slq,
        on="SEQN",
        how="left",
        validate="one_to_one",
    )

    df = df.merge(
        ghb,
        on="SEQN",
        how="left",
        validate="one_to_one",
    )

    df = df.merge(
        glu,
        on="SEQN",
        how="left",
        validate="one_to_one",
        suffixes=("", "_GLU"),
    )

    df = df.merge(
        fast,
        on="SEQN",
        how="left",
        validate="one_to_one",
    )

    df = df.merge(
        ins,
        on="SEQN",
        how="left",
        validate="one_to_one",
        suffixes=("", "_INS"),
    )

    # ------------------------------------------------------------
    # Adult MEC cohort
    # ------------------------------------------------------------

    df = df.loc[
        (df["RIDAGEYR"] >= 18) &
        (df["RIDSTATR"] == 2)
    ].copy()

    # ------------------------------------------------------------
    # Parse weekday sleep timing
    # ------------------------------------------------------------

    df["sleep_onset_min"] = df["SLQ300"].apply(parse_time)
    df["wake_time_min"] = df["SLQ310"].apply(parse_time)

    # Sleep duration derived directly from timing
    df["sleep_duration_derived"] = [
        sleep_interval_hours(s, w)
        for s, w in zip(
            df["sleep_onset_min"],
            df["wake_time_min"]
        )
    ]

    # Sleep midpoint
    df["sleep_midpoint_min"] = [
        circular_midpoint(s, w)
        for s, w in zip(
            df["sleep_onset_min"],
            df["wake_time_min"]
        )
    ]

    df["sleep_midpoint_hour"] = (
        df["sleep_midpoint_min"] / 60
    )

    # ------------------------------------------------------------
    # Compare reported vs derived duration
    # ------------------------------------------------------------

    df["sleep_duration_difference"] = (
        df["SLD012"] -
        df["sleep_duration_derived"]
    )

    # ------------------------------------------------------------
    # Plausibility flags
    # ------------------------------------------------------------

    df["valid_sleep_timing"] = (
        df["sleep_onset_min"].notna() &
        df["wake_time_min"].notna()
    )

    df["valid_sleep_duration"] = (
        df["SLD012"].between(3, 14)
    )

    df["valid_derived_duration"] = (
        df["sleep_duration_derived"].between(3, 14)
    )

    # Timing-derived duration should approximately agree
    df["duration_consistent"] = (
        df["sleep_duration_difference"].abs() <= 1
    )

    # HbA1c available
    df["hba1c_available"] = df["LBXGH"].notna()

    # Fasting glucose available
    df["fasting_glucose_available"] = (
        df["LBXGLU"].notna()
    )

    # Insulin available
    df["insulin_available"] = (
        df["LBXIN"].notna()
    )

    df["cycle"] = cycle

    return df


datasets = []

for cycle, suffix in [
    ("2015-2016", "I"),
    ("2017-2018", "J"),
]:

    print("\n" + "=" * 80)
    print(cycle)
    print("=" * 80)

    df = load_cycle(cycle, suffix)

    print(f"Adult MEC cohort: {len(df):,}")
    print(
        "Valid sleep timing:",
        f"{df['valid_sleep_timing'].sum():,}"
    )
    print(
        "Valid sleep duration:",
        f"{df['valid_sleep_duration'].sum():,}"
    )
    print(
        "HbA1c available:",
        f"{df['hba1c_available'].sum():,}"
    )
    print(
        "Fasting glucose available:",
        f"{df['fasting_glucose_available'].sum():,}"
    )

    print("\nSleep duration comparison:")
    print(
        df[
            [
                "SLD012",
                "sleep_duration_derived",
                "sleep_duration_difference",
            ]
        ].describe().to_string()
    )

    print("\nSleep midpoint:")
    print(
        df["sleep_midpoint_hour"]
        .describe()
        .to_string()
    )

    datasets.append(df)


final = pd.concat(
    datasets,
    ignore_index=True
)

# Global uniqueness check
assert final["SEQN"].is_unique is True

output = OUT / "core_sleep_glycemic.csv"

final.to_csv(
    output,
    index=False
)

print("\n" + "=" * 80)
print("CORE DATASET CREATED")
print("=" * 80)

print(f"Rows: {len(final):,}")
print(f"Columns: {len(final.columns):,}")
print(f"Saved: {output}")

print("\nRows by cycle:")
print(
    final["cycle"]
    .value_counts()
    .sort_index()
)

print("\nMissingness:")
print(
    final[
        [
            "sleep_onset_min",
            "wake_time_min",
            "sleep_midpoint_min",
            "SLD012",
            "LBXGH",
            "LBXGLU",
            "LBXIN",
        ]
    ]
    .isna()
    .mean()
    .mul(100)
    .round(2)
)
