from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "processed" / "core_sleep_glycemic.csv"
OUT = ROOT / "data" / "processed" / "sleep_qc"

OUT.mkdir(parents=True, exist_ok=True)


df = pd.read_csv(INPUT)


# ------------------------------------------------------------
# 1. Basic sleep-duration QC
# ------------------------------------------------------------

df["derived_duration_valid"] = (
    df["sleep_duration_derived"].between(3, 14)
)

df["reported_duration_valid"] = (
    df["SLD012"].between(3, 14)
)

df["duration_consistent_30m"] = (
    df["sleep_duration_difference"].abs() <= 0.5
)

df["duration_consistent_1h"] = (
    df["sleep_duration_difference"].abs() <= 1.0
)


# ------------------------------------------------------------
# 2. Identify extreme / suspicious sleep intervals
# ------------------------------------------------------------

extreme_duration = df.loc[
    (df["sleep_duration_derived"] < 3) |
    (df["sleep_duration_derived"] > 14),
    [
        "SEQN",
        "cycle",
        "SLQ300",
        "SLQ310",
        "SLD012",
        "sleep_onset_min",
        "wake_time_min",
        "sleep_duration_derived",
        "sleep_duration_difference",
    ]
].sort_values(
    "sleep_duration_derived"
)

extreme_duration.to_csv(
    OUT / "extreme_sleep_duration.csv",
    index=False
)


# ------------------------------------------------------------
# 3. Largest discrepancies between reported and derived duration
# ------------------------------------------------------------

largest_discrepancies = (
    df.loc[
        df["SLD012"].notna() &
        df["sleep_duration_derived"].notna(),
        [
            "SEQN",
            "cycle",
            "SLQ300",
            "SLQ310",
            "SLD012",
            "sleep_duration_derived",
            "sleep_duration_difference",
        ]
    ]
    .assign(
        abs_difference=lambda x:
            x["sleep_duration_difference"].abs()
    )
    .sort_values(
        "abs_difference",
        ascending=False
    )
)

largest_discrepancies.head(100).to_csv(
    OUT / "largest_duration_discrepancies.csv",
    index=False
)


# ------------------------------------------------------------
# 4. Duration agreement by cycle
# ------------------------------------------------------------

agreement = (
    df.groupby("cycle")
    .agg(
        n=("SEQN", "size"),
        reported_duration_available=(
            "SLD012",
            lambda x: x.notna().sum()
        ),
        derived_duration_available=(
            "sleep_duration_derived",
            lambda x: x.notna().sum()
        ),
        consistent_30m=(
            "duration_consistent_30m",
            "sum"
        ),
        consistent_1h=(
            "duration_consistent_1h",
            "sum"
        ),
        median_abs_difference=(
            "sleep_duration_difference",
            lambda x: x.abs().median()
        ),
        mean_abs_difference=(
            "sleep_duration_difference",
            lambda x: x.abs().mean()
        ),
    )
)

agreement["pct_consistent_30m"] = (
    agreement["consistent_30m"] /
    agreement["reported_duration_available"]
    * 100
)

agreement["pct_consistent_1h"] = (
    agreement["consistent_1h"] /
    agreement["reported_duration_available"]
    * 100
)

agreement.to_csv(
    OUT / "duration_agreement_by_cycle.csv"
)


# ------------------------------------------------------------
# 5. Sleep midpoint distribution
# ------------------------------------------------------------

midpoint = df.loc[
    df["sleep_midpoint_min"].notna(),
    [
        "SEQN",
        "cycle",
        "sleep_midpoint_min",
        "sleep_midpoint_hour",
    ]
].copy()


# Circular representation
midpoint["midpoint_sin"] = np.sin(
    2 * np.pi *
    midpoint["sleep_midpoint_min"] /
    1440
)

midpoint["midpoint_cos"] = np.cos(
    2 * np.pi *
    midpoint["sleep_midpoint_min"] /
    1440
)


midpoint.to_csv(
    OUT / "sleep_midpoint_circular.csv",
    index=False
)


# ------------------------------------------------------------
# 6. Primary exposure eligibility
# ------------------------------------------------------------

df["primary_sleep_exposure_available"] = (
    df["valid_sleep_timing"] &
    df["valid_sleep_duration"] &
    df["duration_consistent_1h"]
)

df["primary_outcome_available"] = (
    df["LBXGH"].notna()
)

df["primary_analysis_candidate"] = (
    df["primary_sleep_exposure_available"] &
    df["primary_outcome_available"]
)


# ------------------------------------------------------------
# 7. Eligibility by cycle
# ------------------------------------------------------------

eligibility = (
    df.groupby("cycle")
    .agg(
        adults=("SEQN", "size"),
        valid_sleep_timing=(
            "valid_sleep_timing",
            "sum"
        ),
        valid_sleep_duration=(
            "valid_sleep_duration",
            "sum"
        ),
        duration_consistent=(
            "duration_consistent_1h",
            "sum"
        ),
        hba1c_available=(
            "hba1c_available",
            "sum"
        ),
        primary_exposure_available=(
            "primary_sleep_exposure_available",
            "sum"
        ),
        primary_analysis_candidate=(
            "primary_analysis_candidate",
            "sum"
        ),
    )
)

for col in [
    "valid_sleep_timing",
    "valid_sleep_duration",
    "duration_consistent",
    "hba1c_available",
    "primary_exposure_available",
    "primary_analysis_candidate",
]:
    eligibility[f"{col}_pct"] = (
        eligibility[col] /
        eligibility["adults"] *
        100
    )

eligibility.to_csv(
    OUT / "eligibility_by_cycle.csv"
)


# ------------------------------------------------------------
# 8. Sleep midpoint descriptive statistics
# ------------------------------------------------------------

midpoint_summary = (
    df.groupby("cycle")["sleep_midpoint_hour"]
    .describe()
)

midpoint_summary.to_csv(
    OUT / "sleep_midpoint_summary.csv"
)


# ------------------------------------------------------------
# 9. Print results
# ------------------------------------------------------------

print("\n" + "=" * 80)
print("SLEEP QC")
print("=" * 80)

print("\nExtreme derived durations:")
print(f"Count: {len(extreme_duration):,}")

if len(extreme_duration):
    print(
        extreme_duration[
            [
                "SEQN",
                "cycle",
                "SLQ300",
                "SLQ310",
                "SLD012",
                "sleep_duration_derived",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


print("\n" + "-" * 80)
print("DURATION AGREEMENT")
print("-" * 80)

print(
    agreement[
        [
            "n",
            "reported_duration_available",
            "consistent_30m",
            "pct_consistent_30m",
            "consistent_1h",
            "pct_consistent_1h",
            "median_abs_difference",
            "mean_abs_difference",
        ]
    ].round(3).to_string()
)


print("\n" + "-" * 80)
print("PRIMARY ANALYSIS ELIGIBILITY")
print("-" * 80)

print(
    eligibility.round(2).to_string()
)


print("\n" + "-" * 80)
print("SLEEP MIDPOINT")
print("-" * 80)

print(
    midpoint_summary.round(3).to_string()
)


print("\n" + "-" * 80)
print("CIRCULAR MIDPOINT RANGE")
print("-" * 80)

print(
    midpoint[
        [
            "midpoint_sin",
            "midpoint_cos",
        ]
    ].describe().round(4).to_string()
)


print("\n" + "=" * 80)
print("QC FILES CREATED")
print("=" * 80)

for path in sorted(OUT.glob("*")):
    print(path)
