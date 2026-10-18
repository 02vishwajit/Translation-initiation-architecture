from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_ushuffle_species_summary.tsv"
OUT = Path("14_expanded_study/results/gee_publication")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")

labels = {
    "ZM4": "Z. mobilis ZM4",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "CT_ATCC27405": "A. thermocellus",
}
order = list(labels)
df = df.set_index("species").loc[order].reset_index()
df["label"] = df["species"].map(labels)

n = df["n_sequences"].to_numpy()
observed = df["observed_prevalence_percent"].to_numpy()
null_mean = df["null_mean_prevalence_percent"].to_numpy()
null_lo = 100 * df["null_2_5_percentile_positive"].to_numpy() / n
null_hi = 100 * df["null_97_5_percentile_positive"].to_numpy() / n
y = np.arange(len(df))

fig, ax = plt.subplots(figsize=(10, 6.6))
offset = 0.13

# Null distribution: mean and central 95% range.
ax.errorbar(
    null_mean, y + offset,
    xerr=np.vstack([null_mean - null_lo, null_hi - null_mean]),
    fmt="o", capsize=3, markersize=5,
    label="Shuffled null mean (central 95% range)"
)

# Observed prevalence and explicit count/percentage labels.
ax.scatter(
    observed, y - offset,
    marker="D", s=48,
    label="Observed prevalence"
)

for i, row in df.iterrows():
    ax.annotate(
        f'{row["observed_prevalence_percent"]:.1f}% '
        f'({int(row["observed_positive"])}/{int(row["n_sequences"])})',
        (row["observed_prevalence_percent"], i - offset),
        xytext=(7, 0), textcoords="offset points",
        va="center", ha="left", fontsize=8.5
    )

ax.set_yticks(y)
ax.set_yticklabels(df["label"])
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.set_xlabel("Candidate motif-positive upstream windows (%)")
ax.set_title("Observed motif prevalence versus dinucleotide-preserving null")
ax.grid(axis="x", linestyle=":", alpha=0.4)
ax.legend(
    frameon=False, loc="upper center",
    bbox_to_anchor=(0.5, -0.16), ncol=1
)
fig.subplots_adjust(left=0.27, right=0.97, top=0.90, bottom=0.25)

fig.savefig(OUT / "Figure4_dinucleotide_null.pdf", bbox_inches="tight")
fig.savefig(
    OUT / "Figure4_dinucleotide_null.png",
    dpi=600, bbox_inches="tight"
)
plt.close(fig)

print("Updated Figure 4; source summary and statistics unchanged.")
print("Saved:")
print(OUT / "Figure4_dinucleotide_null.pdf")
print(OUT / "Figure4_dinucleotide_null.png")
