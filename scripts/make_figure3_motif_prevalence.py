from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT = Path("14_expanded_study/results/gee_publication")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")
df["positive"] = (
    df["has_candidate_spacer_4_12"]
    .astype(str).str.strip().str.lower().eq("true")
).astype(int)

labels = {
    "CT_ATCC27405": "A. thermocellus",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "ZM4": "Z. mobilis ZM4",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
}
order = [
    "ZM4", "ZM_NCIMB11163", "NOVO", "RSP",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

rows = []
for species in order:
    sub = df[df["species"] == species]
    n = len(sub)
    k = int(sub["positive"].sum())
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    rows.append({
        "species": species,
        "label": labels[species],
        "positive_n": k,
        "total_n": n,
        "prevalence_percent": 100 * k / n,
        "Wilson_95CI_lower_percent": 100 * lo,
        "Wilson_95CI_upper_percent": 100 * hi,
    })

summary = pd.DataFrame(rows)
if not (summary["total_n"] == 348).all():
    raise SystemExit(
        "Expected 348 orthogroups per species. Check dataset and filtering."
    )

summary.to_csv(
    OUT / "Figure3_species_motif_prevalence.tsv",
    sep="\t", index=False
)

fig, ax = plt.subplots(figsize=(9, 5.8))
y = np.arange(len(summary))
x = summary["prevalence_percent"].to_numpy()
lo = summary["Wilson_95CI_lower_percent"].to_numpy()
hi = summary["Wilson_95CI_upper_percent"].to_numpy()

ax.errorbar(
    x, y,
    xerr=np.vstack([x - lo, hi - x]),
    fmt="o", capsize=4, linewidth=1.5, markersize=6
)

ax.set_yticks(y)
ax.set_yticklabels(summary["label"])
ax.set_xlabel("Candidate motif-positive upstream windows (%)")
ax.set_title("Candidate upstream motif prevalence across six bacterial genomes")
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.grid(axis="x", linestyle=":", alpha=0.45)
ax.invert_yaxis()

for i, row in summary.iterrows():
    ax.annotate(
        f'{int(row["positive_n"])}/{int(row["total_n"])} '
        f'({row["prevalence_percent"]:.1f}%)',
        (row["Wilson_95CI_upper_percent"], i),
        xytext=(5, 0), textcoords="offset points",
        va="center", fontsize=8
    )

fig.tight_layout()
fig.savefig(OUT / "Figure3_species_motif_prevalence.pdf", bbox_inches="tight")
fig.savefig(
    OUT / "Figure3_species_motif_prevalence.png",
    dpi=600, bbox_inches="tight"
)
plt.close(fig)

print("\nSpecies-wide motif prevalence")
print(summary.to_string(index=False))
print("\nSaved figure and table to:", OUT)
