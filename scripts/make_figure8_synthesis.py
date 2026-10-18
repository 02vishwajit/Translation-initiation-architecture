from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path("14_expanded_study")
FEATURES = ROOT / "features/single_copy_upstream"
OUT = ROOT / "results/gee_publication"
OUT.mkdir(parents=True, exist_ok=True)

species_order = [
    "ZM4", "ZM_NCIMB11163", "NOVO", "RSP",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

labels = {
    "ZM4": "Z. mobilis ZM4",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "CT_ATCC27405": "A. thermocellus",
}

# -------------------------------
# Panel A: observed motif prevalence
# -------------------------------
feature_file = FEATURES / "single_copy_sequence_features_corrected.tsv"
df = pd.read_csv(feature_file, sep="\t")

required = {"species", "has_candidate_spacer_4_12"}
if not required.issubset(df.columns):
    raise ValueError(f"Missing feature columns: {required - set(df.columns)}")

df["motif_positive"] = (
    df["has_candidate_spacer_4_12"]
    .astype(str).str.strip().str.lower()
    .map({"true": 1, "false": 0})
)

if df["motif_positive"].isna().any():
    raise ValueError("Unexpected values in has_candidate_spacer_4_12.")

summary = (
    df.groupby("species")["motif_positive"]
    .agg(positive="sum", n="count", prevalence="mean")
    .reindex(species_order)
)

if summary.isna().any().any():
    raise ValueError("Species missing from feature summary.")

summary["prevalence_percent"] = summary["prevalence"] * 100
summary["label"] = [labels[s] for s in species_order]
summary.to_csv(OUT / "Figure8_motif_prevalence.tsv", sep="\t", index_label="species")

# -------------------------------
# Panel B: adjusted GEE vs ZM4
# -------------------------------
gee_file = FEATURES / "adjusted_gee_all_30_pairwise.tsv"
gee = pd.read_csv(gee_file, sep="\t")

required_gee = {
    "model", "species_A", "species_B", "odds_ratio_A_vs_B",
    "OR_CI95_lower", "OR_CI95_upper", "Wald_p_value",
    "BH_FDR_across_30_tests"
}
if not required_gee.issubset(gee.columns):
    raise ValueError(f"Missing GEE columns: {required_gee - set(gee.columns)}")

gee = gee[gee["model"] == "all_windows_adjusted"].copy()

# Keep rows that compare each target species directly against ZM4.
targets = [
    "ZM_NCIMB11163", "NOVO", "RSP",
    "CT_ATCC27405", "TSA_JW_SL_YS485"
]
rows = []

for sp in targets:
    direct = gee[(gee["species_A"] == sp) & (gee["species_B"] == "ZM4")]
    reverse = gee[(gee["species_A"] == "ZM4") & (gee["species_B"] == sp)]

    if len(direct) == 1:
        r = direct.iloc[0].copy()
        # Direct row already represents target species / ZM4.
    elif len(reverse) == 1:
        r = reverse.iloc[0].copy()
        # Reverse the log-OR and invert the CI to get target species / ZM4.
        old_or = float(r["odds_ratio_A_vs_B"])
        old_low = float(r["OR_CI95_lower"])
        old_high = float(r["OR_CI95_upper"])
        r["odds_ratio_A_vs_B"] = 1.0 / old_or
        r["OR_CI95_lower"] = 1.0 / old_high
        r["OR_CI95_upper"] = 1.0 / old_low
        r["species_A"], r["species_B"] = sp, "ZM4"
    else:
        raise ValueError(f"Could not uniquely find comparison for {sp} vs ZM4.")

    rows.append(r)

panel = pd.DataFrame(rows)
panel["species_label"] = panel["species_A"].map(labels)
panel.to_csv(OUT / "Figure8_adjusted_GEE.tsv", sep="\t", index=False)

# -------------------------------
# Draw two-panel figure
# -------------------------------
fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(14, 7),
    gridspec_kw={"width_ratios": [1.05, 1.2]}
)

# Panel A
y1 = np.arange(len(summary))
ax1.barh(y1, summary["prevalence_percent"])
ax1.set_yticks(y1)
ax1.set_yticklabels(summary["label"])
ax1.invert_yaxis()
ax1.set_xlim(0, 100)
ax1.set_xlabel("Candidate motif-positive upstream windows (%)")
ax1.set_title("A. Observed motif prevalence")

for i, (_, row) in enumerate(summary.iterrows()):
    ax1.text(
        row["prevalence_percent"] + 1.5, i,
        f'{int(row["positive"])}/{int(row["n"])}',
        va="center", fontsize=8
    )

# Panel B
y2 = np.arange(len(panel))
ors = panel["odds_ratio_A_vs_B"].astype(float).to_numpy()
low = panel["OR_CI95_lower"].astype(float).to_numpy()
high = panel["OR_CI95_upper"].astype(float).to_numpy()
fdr = panel["BH_FDR_across_30_tests"].astype(float).to_numpy()

ax2.errorbar(
    ors, y2,
    xerr=[ors - low, high - ors],
    fmt="o", capsize=3
)
ax2.axvline(1, color="grey", linestyle="--", linewidth=1)
ax2.set_xscale("log")
ax2.set_yticks(y2)
ax2.set_yticklabels(panel["species_label"])
ax2.invert_yaxis()
ax2.set_xlabel("Adjusted odds ratio vs ZM4 (log scale)")
ax2.set_title("B. Adjusted GEE comparisons")

for i, (_, row) in enumerate(panel.iterrows()):
    ax2.text(
        row["OR_CI95_upper"] * 1.08, i,
        f'OR {row["odds_ratio_A_vs_B"]:.2f}; FDR {row["BH_FDR_across_30_tests"]:.2g}',
        va="center", fontsize=7
    )

fig.suptitle(
    "Candidate translation-initiation-associated sequence features",
    fontsize=13
)
fig.tight_layout(rect=[0, 0, 1, 0.95])

fig.savefig(OUT / "Figure8_quantitative_synthesis.pdf", bbox_inches="tight")
fig.savefig(OUT / "Figure8_quantitative_synthesis.png", dpi=600, bbox_inches="tight")
plt.close(fig)

print("\nObserved motif prevalence:")
print(summary[["positive", "n", "prevalence_percent"]].round(2).to_string())

print("\nAdjusted comparisons against ZM4:")
print(panel[[
    "species_A", "odds_ratio_A_vs_B", "OR_CI95_lower",
    "OR_CI95_upper", "Wald_p_value", "BH_FDR_across_30_tests"
]].to_string(index=False))

print(f"\nSaved figure and tables in: {OUT}")
