import csv
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT_GLOBAL = BASE / "adjusted_gee_joint_species_tests.tsv"
OUT_PAIRS = BASE / "adjusted_gee_species_vs_ZM4_FDR.tsv"

df = pd.read_csv(INPUT, sep="\t")
df["positive"] = (
    df["has_candidate_spacer_4_12"].astype(str).str.lower().eq("true")
).astype(int)
df["overlap"] = (
    df["overlap_flag"].astype(str).str.lower().eq("true")
).astype(int)
df["gc"] = pd.to_numeric(df["GC_percent"], errors="raise")
df["gc_z"] = (
    df["gc"] - df["gc"].mean()
) / df["gc"].std(ddof=0)
df["start_group"] = np.where(
    df["start_codon_class"].eq("ATG"), "ATG", "Non_ATG"
)
df["species"] = df["species"].astype(str)
df["orthogroup"] = df["orthogroup"].astype(str)

formula_all = (
    "positive ~ C(species, Treatment(reference='ZM4'))"
    " + gc_z + C(start_group) + overlap"
)
formula_clean = (
    "positive ~ C(species, Treatment(reference='ZM4'))"
    " + gc_z + C(start_group)"
)

models = [
    ("all_windows_adjusted", df, formula_all),
    ("nonoverlap_adjusted", df.loc[df["overlap"] == 0].copy(), formula_clean),
]

global_rows = []
pair_rows = []

for name, data, formula in models:
    result = smf.gee(
        formula=formula,
        groups="orthogroup",
        data=data,
        family=sm.families.Binomial(),
        cov_struct=sm.cov_struct.Exchangeable(),
    ).fit(maxiter=200, cov_type="robust")

    if not getattr(result, "converged", False):
        raise SystemExit(f"{name}: model did not converge.")

    terms = [
        term for term in result.params.index
        if term.startswith("C(species, Treatment(reference='ZM4'))[T.")
    ]
    if len(terms) != 5:
        raise SystemExit(f"{name}: expected 5 species coefficients, found {len(terms)}")

    beta = result.params.loc[terms].to_numpy(dtype=float)
    covariance = result.cov_params().loc[terms, terms].to_numpy(dtype=float)

    if np.linalg.matrix_rank(covariance) != len(terms):
        raise SystemExit(f"{name}: species covariance matrix is rank deficient.")

    statistic = float(beta @ np.linalg.solve(covariance, beta))
    df_test = len(terms)
    p_global = float(chi2.sf(statistic, df_test))

    global_rows.append({
        "model": name,
        "n_observations": int(result.nobs),
        "n_orthogroup_clusters": int(data["orthogroup"].nunique()),
        "wald_chi_square": statistic,
        "df": df_test,
        "global_species_p_value": p_global,
        "converged": result.converged,
    })

    print(f"\n{name}")
    print(f"Joint species Wald chi-square({df_test}) = {statistic:.6f}")
    print(f"Global p-value = {p_global:.8g}")

    for term in terms:
        beta_i = float(result.params[term])
        lo, hi = map(float, result.conf_int().loc[term])
        pair_rows.append({
            "model": name,
            "term": term,
            "odds_ratio": float(np.exp(beta_i)),
            "OR_CI95_lower": float(np.exp(lo)),
            "OR_CI95_upper": float(np.exp(hi)),
            "Wald_p_value": float(result.pvalues[term]),
        })

# BH correction across all ten species-versus-ZM4 tests (both models).
order = sorted(
    range(len(pair_rows)),
    key=lambda i: pair_rows[i]["Wald_p_value"]
)
m = len(pair_rows)
adjusted = [1.0] * m
running = 1.0

for j in range(m - 1, -1, -1):
    i = order[j]
    running = min(
        running,
        pair_rows[i]["Wald_p_value"] * m / (j + 1)
    )
    adjusted[i] = min(1.0, running)

for i, row in enumerate(pair_rows):
    row["BH_FDR_across_10_tests"] = adjusted[i]

for path, rows in [
    (OUT_GLOBAL, global_rows),
    (OUT_PAIRS, pair_rows),
]:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

print("\nADJUSTED SPECIES-VERSUS-ZM4 COMPARISONS")
print("Model\tTerm\tOR\t95% CI\tP\tBH-FDR")
for r in pair_rows:
    print(
        f"{r['model']}\t{r['term']}\t{r['odds_ratio']:.3f}\t"
        f"{r['OR_CI95_lower']:.3f}-{r['OR_CI95_upper']:.3f}\t"
        f"{r['Wald_p_value']:.6g}\t"
        f"{r['BH_FDR_across_10_tests']:.6g}"
    )

print("\nSaved:")
print(OUT_GLOBAL)
print(OUT_PAIRS)
