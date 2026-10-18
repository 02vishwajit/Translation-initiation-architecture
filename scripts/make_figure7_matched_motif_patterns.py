from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT = Path("14_expanded_study/results/gee_publication")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")
df["positive"] = (
    df["has_candidate_spacer_4_12"]
    .astype(str).str.strip().str.lower().eq("true")
).astype(int)

species_order = [
    "ZM4", "ZM_NCIMB11163", "NOVO",
    "RSP", "TSA_JW_SL_YS485", "CT_ATCC27405"
]

wide = df.pivot(index="orthogroup", columns="species", values="positive")

if set(wide.columns) != set(species_order):
    raise SystemExit("Species columns do not match expected six-species panel.")
if len(wide) != 348 or wide.isna().any().any():
    raise SystemExit("Expected 348 complete orthogroups; check input data.")

wide = wide[species_order]
wide["motif_positive_species_count"] = wide.sum(axis=1).astype(int)

counts = (
    wide["motif_positive_species_count"]
    .value_counts()
    .reindex(range(7), fill_value=0)
    .rename_axis("positive_species_count")
    .reset_index(name="orthogroup_count")
)
counts["percentage"] = 100 * counts["orthogroup_count"] / len(wide)

counts.to_csv(
    OUT / "Figure7_shared_motif_count_distribution.tsv",
    sep="\t", index=False
)

# Heatmap: orthogroups grouped by how many species are positive.
wide = wide.sort_values("motif_positive_species_count", ascending=False)
matrix = wide[species_order].to_numpy()

fig, axes = plt.subplots(
    1, 2, figsize=(11, 8),
    gridspec_kw={"width_ratios": [1.1, 3.4]}
)

x = counts["positive_species_count"].to_numpy()
y = counts["orthogroup_count"].to_numpy()
axes[0].barh(x, y)
axes[0].set_yticks(range(7))
axes[0].set_ylabel("Species with candidate motif (out of 6)")
axes[0].set_xlabel("Orthogroup count")
axes[0].set_title("Distribution")

axes[1].imshow(
    matrix, aspect="auto", interpolation="nearest",
    cmap="Greys", vmin=0, vmax=1
)
axes[1].set_xticks(range(6))
axes[1].set_xticklabels(species_order, rotation=45, ha="right")
axes[1].set_ylabel("Orthogroups (sorted by motif-positive count)")
axes[1].set_title("Matched motif status")
axes[1].set_xlabel("Species")

fig.suptitle(
    "Candidate upstream motif status across matched single-copy orthogroups",
    y=0.99
)
fig.tight_layout()
fig.savefig(OUT / "Figure7_matched_motif_patterns.pdf", bbox_inches="tight")
fig.savefig(
    OUT / "Figure7_matched_motif_patterns.png",
    dpi=600, bbox_inches="tight"
)
plt.close(fig)

print(counts.to_string(index=False))
print("\nSaved Figure 7 PDF, PNG and summary TSV to:", OUT)
