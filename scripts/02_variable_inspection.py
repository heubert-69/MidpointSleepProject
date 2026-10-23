from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

FILES = {
    "2015-2016": {
        "DEMO": RAW / "2015-2016" / "DEMO_I.XPT",
        "SLQ": RAW / "2015-2016" / "SLQ_I.XPT",
        "GHB": RAW / "2015-2016" / "GHB_I.XPT",
        "GLU": RAW / "2015-2016" / "GLU_I.XPT",
        "FASTQX": RAW / "2015-2016" / "FASTQX_I.XPT",
        "INS": RAW / "2015-2016" / "INS_I.XPT",
    },
    "2017-2018": {
        "DEMO": RAW / "2017-2018" / "DEMO_J.XPT",
        "SLQ": RAW / "2017-2018" / "SLQ_J.XPT",
        "GHB": RAW / "2017-2018" / "GHB_J.XPT",
        "GLU": RAW / "2017-2018" / "GLU_J.XPT",
        "FASTQX": RAW / "2017-2018" / "FASTQX_J.XPT",
        "INS": RAW / "2017-2018" / "INS_J.XPT",
    },
}


VARIABLES = {
    "DEMO": [
        "SEQN",
        "RIDAGEYR",
        "RIDSTATR",
        "RIAGENDR",
        "RIDRETH1",
        "RIDRETH3",
        "DMDEDUC2",
        "DMDMARTL",
        "INDFMPIR",
        "WTINT2YR",
        "WTMEC2YR",
        "SDMVPSU",
        "SDMVSTRA",
    ],

    "SLQ": [
        "SEQN",
        "SLQ300",
        "SLQ310",
        "SLQ320",
        "SLQ330",
        "SLD012",
        "SLD013",
        "SLQ030",
        "SLQ040",
        "SLQ050",
        "SLQ120",
    ],

    "GHB": [
        "SEQN",
        "LBXGH",
    ],

    "GLU": [
        "SEQN",
        "WTSAF2YR",
        "LBXGLU",
        "LBDGLUSI",
    ],

    "FASTQX": [
        "SEQN",
        "PHAFSTHR",
        "PHAFSTMN",
        "PHDSESN",
    ],

    "INS": [
        "SEQN",
        "WTSAF2YR",
        "LBXIN",
        "LBDINSI",
        "LBDINLC",
        "PHAFSTHR",
        "PHAFSTMN",
    ],
}


def inspect_file(path, variables):

    print("\n" + "=" * 80)
    print(path.name)
    print("=" * 80)

    df = pd.read_sas(path, format="xport")

    available = [
        v for v in variables
        if v in df.columns
    ]

    missing_variables = [
        v for v in variables
        if v not in df.columns
    ]

    print(f"Shape: {df.shape}")

    if missing_variables:
        print("\nMissing requested variables:")
        print(missing_variables)

    for variable in available:

        s = df[variable]

        print("\n" + "-" * 60)
        print(variable)
        print("-" * 60)

        print("dtype:", s.dtype)
        print("missing:", f"{s.isna().mean() * 100:.2f}%")
        print("unique non-null:", s.nunique(dropna=True))

        if s.nunique(dropna=True) <= 30:

            print("\nValue counts:")
            print(
                s.value_counts(
                    dropna=False
                ).sort_index().to_string()
            )

        elif variable in [
            "SLQ300",
            "SLQ310",
            "SLQ320",
            "SLQ330",
            "SLD012",
            "SLD013",
            "RIDAGEYR",
            "LBXGH",
            "LBXGLU",
            "LBXIN",
            "PHAFSTHR",
            "PHAFSTMN",
        ]:

            print("\nDescriptive statistics:")
            print(
                s.describe(
                    percentiles=[
                        0.01,
                        0.05,
                        0.25,
                        0.50,
                        0.75,
                        0.95,
                        0.99,
                    ]
                ).to_string()
            )


for cycle, files in FILES.items():

    print("\n\n")
    print("#" * 100)
    print(f"# {cycle}")
    print("#" * 100)

    for name, path in files.items():

        if not path.exists():
            print(f"\nMISSING FILE: {path}")
            continue

        inspect_file(
            path,
            VARIABLES[name]
        )
