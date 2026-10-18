import csv
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrix
from statsmodels.genmod.cov_struct import Exchangeable

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT_ALL = BASE / "adjusted_motif_gee_all.tsv"
OUT_CLEAN = BASE / "adjusted_motif_gee_nonoverlap.tsv"
OUT_DIAG = BASE / "adjusted_motif_gee_diagnostics.tsv"

df = pd.read_csv(INPUT, sep="\t")

required = [
    "orthogroup", "species", "GC_percent",
    "start_codon_class", "overlap_flag",
    "has_candidate_spacer_4_12",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns: {missing}")

df["positive"] = (
    df["has_candidate_spacer_4_12"].astype(str).str.lower().eq("true")
).astype(int)

df["overlap"] = (
    df["overlap_flag"].astype(str).str.lower().eq("true")
).astype(int)

df["gc"] = pd.to_numeric(df["GC_percent"], errors="coerce")

# Only 10 records have "Other_or_unexpected" starts.
# Collapse all non-ATG starts to avoid a very sparse category.
df["start_group"] = np.where(
    df["start_codon_class"].eq("ATG"), "ATG", "Non_ATG"
)

df["species"] = df["species"].astype(str)
df["orthogroup"] = df["orthogroup"].astype(str)

if df[required].isna().any().any() or df["gc"].isna().any():
    raise SystemExit("Missing values detected; inspect before modeling.")

if df["orthogroup"].nunique() != 348:
    raise SystemExit("Expected 348 orthogroup clusters.")

# Center and scale GC percentage for numerical stability.
gc_mean = df["gc"].mean()
gc_sd = df["gc"].std(ddof=0)
if gc_sd == 0:
    raise SystemExit("GC has zero variance.")
df["gc_z"] = (df["gc"] - gc_mean) / gc_sd

# Set ZM4 as the reference species.
species_levels = sorted(df["species"].unique())
if "ZM4" not in species_levels:
    raise SystemExit("ZM4 reference species not found.")

diagnostics = []

def fit_gee(data, include_overlap, outfile, model_name):
    data = data.copy()
    formula = (
        "positive ~ C(species, Treatment(reference='ZM4'))"
        " + gc_z + C(start_group)"
    )
    if include_overlap:
        formula += " + overlap"

    # Check the fixed-effect design matrix for rank deficiency.
    X = dmatrix(formula.split("~", 1)[1], data, return_type="dataframe")
    rank = np.linalg.matrix_rank(X.to_numpy())
    condition_number = np.linalg.cond(X.to_numpy())

    diagnostics.append({
        "model": model_name,
        "n": len(data),
        "clusters": data["orthogroup"].nunique(),
        "design_columns": X.shape[1],
        "design_rank": rank,
        "design_full_rank": rank == X.shape[1],
        "design_condition_number": condition_number,
        "gc_mean_percent_original": gc_mean,
        "gc_sd_percent_original": gc_sd,
        "start_group_counts": str(data["start_group"].value_counts().to_dict()),
        "species_counts": str(data["species"].value_counts().to_dict()),
    })

    if rank < X.shape[1]:
        raise SystemExit(
            f"{model_name}: design matrix is rank deficient; "
            "do not interpret coefficients."
        )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = smf.gee(
            formula=formula,
            groups="orthogroup",
            data=data,
            family=sm.families.Binomial(),
            cov_struct=Exchangeable(),
        )
        result = model.fit(maxiter=200, cov_type="robust")

    print(f"\n{model_name}")
    print("Formula:", formula)
    print("Converged:", getattr(result, "converged", "not reported"))
    print("Observations:", int(result.nobs))
    print("Orthogroup clusters:", data["orthogroup"].nunique())
    print("Design condition number:", f"{condition_number:.3f}")
    print(result.summary())

    if caught:
        print("Warnings:")
        for w in caught:
            print(" -", str(w.message))

    rows = []
    ci = result.conf_int()
    for term in result.params.index:
        beta = float(result.params[term])
        se = float(result.bse[term])
        p = float(result.pvalues[term])
        lo, hi = map(float, ci.loc[term])
        rows.append({
            "model": model_name,
            "term": term,
            "coefficient_log_odds": beta,
            "robust_SE": se,
            "odds_ratio": np.exp(beta),
            "OR_CI95_lower": np.exp(lo),
            "OR_CI95_upper": np.exp(hi),
            "Wald_p_value": p,
            "converged": getattr(result, "converged", None),
            "n": int(result.nobs),
            "clusters": data["orthogroup"].nunique(),
            "design_condition_number": condition_number,
        })

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    print("Saved coefficient table:", outfile)

fit_gee(
    df,
    include_overlap=True,
    outfile=OUT_ALL,
    model_name="all_windows_adjusted",
)

clean = df.loc[df["overlap"] == 0].copy()
if clean["orthogroup"].nunique() < 100:
    raise SystemExit("Too few orthogroups remain in non-overlap analysis.")

fit_gee(
    clean,
    include_overlap=False,
    outfile=OUT_CLEAN,
    model_name="nonoverlap_adjusted",
)

with open(OUT_DIAG, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(diagnostics[0]), delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(diagnostics)

print("\nDiagnostics saved:", OUT_DIAG)
print("GC centering/scaling: mean =", gc_mean, "SD =", gc_sd)
