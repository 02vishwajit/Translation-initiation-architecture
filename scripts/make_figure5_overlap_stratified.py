from pathlib import Path
import pandas as pd
import numpy as np
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
df["overlap_group"] = (
    df["overlap_flag"].astype(str).str.strip().str.lower()
    .map({"true": "Overlapping", "false": "Non-overlapping"})
)

if df["overlap_group"].isna().any():
    raise SystemExit("Unexpected overlap_flag values; inspect source data.")

labels = {
    "ZM4": "Z. mobilis ZM4",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "CT_ATCC27405": "A. thermocellus",
}
order = list(labels)

rows = []
for species in order:
    for group in ["Non-overlapping", "Overlapping"]:
        sub = df[(df["species"] == species) & (df["overlap_group"] == group)]
        n = len(sub)
        k = int(sub["positive"].sum())
        lo, hi = proportion_confint(k, n, method="wilson")
        rows.append({
            "species": species,
            "label": labels[species],
            "overlap_group": group,
            "positive_n": k,
            "total_n": n,
            "prevalence_percent": 100 * k / n,
            "Wilson_95CI_lower_percent": 100 * lo,
            "Wilson_95CI_upper_percent": 100 * hi,
        })

summary = pd.DataFrame(rows)
summary.to_csv(
    OUT / "Figure5_overlap_stratified_prevalence.tsv",
    sep="\t", index=False
)

fig, ax = plt.subplots(figsize=(10, 6.5))
ybase = np.arange(len(order))
offset = 0.13

for group, dy, marker in [
    ("Non-overlapping", -offset, "o"),
    ("Overlapping", offset, "s"),
]:
    sub = summary[summary["overlap_group"] == group].set_index("species").loc[order]
    x = sub["prevalence_percent"].to_numpy()
    lo = sub["Wilson_95CI_lower_percent"].to_numpy()
    hi = sub["Wilson_95CI_upper_percent"].to_numpy()

    ax.errorbar(
        x, ybase + dy,
        xerr=np.vstack([x - lo, hi - x]),
        fmt=marker, capsize=3, markersize=5,
        label=group
    )

    for i, (_, row) in enumerate(sub.iterrows()):
        ax.annotate(
            f'{row["positive_n"]}/{row["total_n"]}',
            (row["prevalence_percent"], ybase[i] + dy),
            xytext=(5, 0), textcoords="offset points",
            va="center", fontsize=7
        )

ax.set_yticks(ybase)
ax.set_yticklabels([labels[s] for s in order])
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.set_xlabel("Candidate motif-positive upstream windows (%)")
ax.set_title("Candidate motif prevalence by CDS-overlap status")
ax.grid(axis="x", linestyle=":", alpha=0.4)
ax.legend(frameon=False, loc="lower right")
ax.invert_yaxis()
fig.tight_layout()

fig.savefig(OUT / "Figure5_overlap_stratified_prevalence.pdf", bbox_inches="tight")
fig.savefig(
    OUT / "Figure5_overlap_stratified_prevalence.png",
    dpi=600, bbox_inches="tight"
)
plt.close(fig)

print(summary.to_string(index=False))
print("\nSaved Figure 5 PDF, PNG and TSV to:", OUT)
