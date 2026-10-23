from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed" / "audit"

OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {
    "2015-2016": RAW / "2015-2016",
    "2017-2018": RAW / "2017-2018",
}

CANDIDATE_PREFIXES = [
    "SLQ",
    "DEMO",
    "GHB",
    "GLU",
    "FASTQX",
    "INS",
    "BMX",
    "BPX",
    "BPQ",
    "PAQ",
    "ALQ",
    "SMQ",
    "DIQ",
    "MCQ",
    "DPQ",
    "TCHOL",
    "TRIGLY",
    "HDL",
]

audit = []
schemas = {}

for cycle, directory in CYCLES.items():

    print("\n" + "=" * 80)
    print(cycle)
    print("=" * 80)

    if not directory.exists():
        print(f"WARNING: directory missing: {directory}")
        continue

    files = sorted(directory.glob("*.XPT"))

    if not files:
        print("WARNING: no XPT files found")
        continue

    for path in files:

        print(f"\n--- {path.name} ---")

        try:
            df = pd.read_sas(path, format="xport")
        except Exception as e:
            print(f"ERROR: {e}")
            continue

        cols = list(df.columns)

        seqn_exists = "SEQN" in df.columns
        duplicate_seqn = None

        if seqn_exists:
            duplicate_seqn = int(df["SEQN"].duplicated().sum())

        missing = (
            df.isna()
            .mean()
            .mul(100)
            .sort_values(ascending=False)
        )

        print(f"Rows: {len(df):,}")
        print(f"Columns: {len(df):,}")
        print(f"SEQN present: {seqn_exists}")

        if seqn_exists:
            print(f"Duplicate SEQN: {duplicate_seqn}")

        print("\nVariables:")
        print(" ".join(cols))

        print("\nTop missingness:")
        print(missing.head(15).round(2).to_string())

        prefix = path.stem.split("_")[0]

        candidate_cols = [
            c for c in cols
            if any(c.startswith(p) for p in CANDIDATE_PREFIXES)
        ]

        schemas.setdefault(cycle, {})[path.name] = {
            "columns": cols,
            "shape": [int(df.shape[0]), int(df.shape[1])],
            "duplicate_seqn": duplicate_seqn,
            "candidate_columns": candidate_cols,
        }

        audit.append({
            "cycle": cycle,
            "file": path.name,
            "rows": int(df.shape[0]),
            "columns": int(df.shape[1]),
            "seqn_present": seqn_exists,
            "duplicate_seqn": duplicate_seqn,
        })


audit_df = pd.DataFrame(audit)

audit_df.to_csv(
    OUT / "file_audit.csv",
    index=False
)

with open(OUT / "schemas.json", "w") as f:
    json.dump(schemas, f, indent=2, default=str)


# ------------------------------------------------------------------
# Cross-cycle comparison
# ------------------------------------------------------------------

def file_prefix(filename):
    return filename.rsplit(".", 1)[0].rsplit("_", 1)[0]


files_by_cycle = {}

for cycle in schemas:
    files_by_cycle[cycle] = {
        file_prefix(filename): filename
        for filename in schemas[cycle]
    }


common_prefixes = sorted(
    set(files_by_cycle.get("2015-2016", {}))
    &
    set(files_by_cycle.get("2017-2018", {}))
)

comparison = []

for prefix in common_prefixes:

    f_i = files_by_cycle["2015-2016"][prefix]
    f_j = files_by_cycle["2017-2018"][prefix]

    cols_i = set(schemas["2015-2016"][f_i]["columns"])
    cols_j = set(schemas["2017-2018"][f_j]["columns"])

    common = sorted(cols_i & cols_j)
    only_i = sorted(cols_i - cols_j)
    only_j = sorted(cols_j - cols_i)

    comparison.append({
        "prefix": prefix,
        "2015_file": f_i,
        "2017_file": f_j,
        "common_variables": len(common),
        "only_2015_2016": len(only_i),
        "only_2017_2018": len(only_j),
        "common": "|".join(common),
        "only_2015_2016": "|".join(only_i),
        "only_2017_2018": "|".join(only_j),
    })


comparison_df = pd.DataFrame(comparison)

comparison_df.to_csv(
    OUT / "cross_cycle_comparison.csv",
    index=False
)


# ------------------------------------------------------------------
# Explicit sleep-variable audit
# ------------------------------------------------------------------

sleep_variables = [
    "SEQN",
    "SLQ300",
    "SLQ310",
    "SLQ320",
    "SLQ330",
    "SLD012",
    "SLD013",
]

sleep_audit = []

for cycle in CYCLES:

    slq_files = [
        f for f in schemas.get(cycle, {})
        if f.startswith("SLQ_")
    ]

    if not slq_files:
        continue

    filename = slq_files[0]
    columns = set(schemas[cycle][filename]["columns"])

    for variable in sleep_variables:
        sleep_audit.append({
            "cycle": cycle,
            "variable": variable,
            "present": variable in columns,
        })

sleep_df = pd.DataFrame(sleep_audit)

sleep_df.to_csv(
    OUT / "sleep_variable_audit.csv",
    index=False
)


print("\n" + "=" * 80)
print("AUDIT COMPLETE")
print("=" * 80)

print(f"\nSaved:")
print(f"  {OUT / 'file_audit.csv'}")
print(f"  {OUT / 'schemas.json'}")
print(f"  {OUT / 'cross_cycle_comparison.csv'}")
print(f"  {OUT / 'sleep_variable_audit.csv'}")

print("\nSleep-variable comparison:")
print(
    sleep_df.to_string(index=False)
)
