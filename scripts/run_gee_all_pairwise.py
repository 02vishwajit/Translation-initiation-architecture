import csv
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm
from patsy import dmatrix

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT = BASE / "adjusted_gee_all_30_pairwise.tsv"

df = pd.read_csv(INPUT, sep="\t")
df["positive"] = (
    df["has_candidate_spacer_4_12"].astype(str).str.lower().eq("true")
).astype(int)
df["overlap"] = (
    df["overlap_flag"].astype(str).str.lower().eq("true")
).astype(int)
df["gc"] = pd.to_numeric(df["GC_percent"], errors="raise")
df["gc_z"] = (df["gc"] - df["gc"].mean()) / df["gc"].std(ddof=0)
df["start_group"] = np.where(
    df["start_codon_class"].eq("ATG"), "ATG", "Non_ATG"
)
df["species"] = df["species"].astype(str)
df["orthogroup"] = df["orthogroup"].astype(str)

species = sorted(df["species"].unique())
if len(species) != 6:
    raise SystemExit("Expected exactly six species.")

models = [
    (
        "all_windows_adjusted",
        df.copy(),
        "positive ~ C(species, Treatment(reference='ZM4')) + gc_z + C(start_group) + overlap",
    ),
    (
        "nonoverlap_adjusted",
        df.loc[df["overlap"] == 0].copy(),
        "positive ~ C(species, Treatment(reference='ZM4')) + gc_z + C(start_group)",
    ),
]

results = []

for model_name, data, formula in models:
    fit = smf.gee(
        formula=formula,
        groups="orthogroup",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit(maxiter=200, cov_type="robust")

    if not fit.converged:
        raise SystemExit(f"{model_name} did not converge.")

    # Construct design information independently of Statsmodels' PandasData.
    design_formula = formula.split("~", 1)[1].strip()
    design = dmatrix(
        design_formula,
        data,
        return_type="dataframe"
    )
    design_info = design.design_info

    reference = data.iloc[[0]].copy()
    reference["species"] = "ZM4"
    reference["gc_z"] = 0.0
    reference["start_group"] = "ATG"
    reference["overlap"] = 0

    from patsy import build_design_matrices
    covariance = np.asarray(fit.cov_params())
    beta = np.asarray(fit.params)

    # Align Patsy columns to the fitted coefficient order.
    if list(design.columns) != list(fit.params.index):
        raise SystemExit(
            f"Design/coefficient mismatch in {model_name}: "
            f"{list(design.columns)} vs {list(fit.params.index)}"
        )

    for a, b in combinations(species, 2):
        row_a = reference.copy()
        row_b = reference.copy()
        row_a["species"] = a
        row_b["species"] = b

        xa = np.asarray(build_design_matrices([design_info], row_a)[0])[0]
        xb = np.asarray(build_design_matrices([design_info], row_b)[0])[0]
        contrast = xa - xb

        log_or = float(contrast @ beta)
        variance = float(contrast @ covariance @ contrast)
        if variance < -1e-10:
            raise SystemExit(f"Negative contrast variance for {a} vs {b}")
        se = np.sqrt(max(variance, 0.0))

        if se == 0:
            z = 0.0 if log_or == 0 else np.inf
            p = 1.0 if log_or == 0 else 0.0
        else:
            z = log_or / se
            p = float(2 * norm.sf(abs(z)))

        results.append({
            "model": model_name,
            "species_A": a,
            "species_B": b,
            "contrast_definition": "A relative to B; adjusted at mean GC, ATG start; overlap=0 for reference covariates",
            "log_odds_ratio_A_vs_B": log_or,
            "robust_SE": se,
            "odds_ratio_A_vs_B": np.exp(log_or),
            "OR_CI95_lower": np.exp(log_or - 1.96 * se),
            "OR_CI95_upper": np.exp(log_or + 1.96 * se),
            "Wald_p_value": p,
            "n_observations": int(fit.nobs),
            "n_orthogroups": int(data["orthogroup"].nunique()),
        })

# Benjamini-Hochberg adjustment across all 30 comparisons.
order = sorted(range(len(results)), key=lambda i: results[i]["Wald_p_value"])
m = len(results)
adjusted = [1.0] * m
running = 1.0

for j in range(m - 1, -1, -1):
    i = order[j]
    running = min(running, results[i]["Wald_p_value"] * m / (j + 1))
    adjusted[i] = min(1.0, running)

for i, row in enumerate(results):
    row["BH_FDR_across_30_tests"] = adjusted[i]

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(results)

print("ALL PAIRWISE ADJUSTED GEE CONTRASTS")
print("FDR correction: 30 tests across both models")
for r in sorted(results, key=lambda x: x["BH_FDR_across_30_tests"]):
    print(
        f"{r['model']}: {r['species_A']} vs {r['species_B']} | "
        f"OR={r['odds_ratio_A_vs_B']:.3f} "
        f"(95% CI {r['OR_CI95_lower']:.3f}-{r['OR_CI95_upper']:.3f}) | "
        f"p={r['Wald_p_value']:.4g} | "
        f"FDR={r['BH_FDR_across_30_tests']:.4g}"
    )

print("\nSaved:", OUT)
