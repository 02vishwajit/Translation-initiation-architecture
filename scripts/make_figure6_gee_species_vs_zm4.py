from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "adjusted_gee_all_30_pairwise.tsv"
OUT = Path("14_expanded_study/results/gee_publication")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")
df = df[df["model"] == "all_windows_adjusted"].copy()

# Select all five species-versus-ZM4 comparisons, reversing
# odds ratios and confidence limits where ZM4 is species A.
rows = []
comparators = [
    "CT_ATCC27405",
    "NOVO",
    "RSP",
    "TSA_JW_SL_YS485",
    "ZM_NCIMB11163",
]

for comparator in comparators:
    direct = df[
        (df["species_A"] == comparator) &
        (df["species_B"] == "ZM4")
    ]
    reverse = df[
        (df["species_A"] == "ZM4") &
        (df["species_B"] == comparator)
    ]

    if len(direct) == 1:
        row = direct.iloc[0].copy()
    elif len(reverse) == 1:
        row = reverse.iloc[0].copy()
        old_or = row["odds_ratio_A_vs_B"]
        old_lo = row["OR_CI95_lower"]
        old_hi = row["OR_CI95_upper"]

        row["species_A"] = comparator
        row["species_B"] = "ZM4"
        row["log_odds_ratio_A_vs_B"] = -row["log_odds_ratio_A_vs_B"]
        row["odds_ratio_A_vs_B"] = 1.0 / old_or
        row["OR_CI95_lower"] = 1.0 / old_hi
        row["OR_CI95_upper"] = 1.0 / old_lo
    else:
        raise SystemExit(
            f"Could not uniquely find comparison for {comparator} vs ZM4."
        )

    rows.append(row)

df = pd.DataFrame(rows)

labels = {
    "CT_ATCC27405": "A. thermocellus",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
}
df["label"] = df["species_A"].map(labels)

df["significant"] = df["BH_FDR_across_30_tests"] < 0.05
df = df.sort_values("odds_ratio_A_vs_B").reset_index(drop=True)

df.to_csv(
    OUT / "Figure6_species_vs_ZM4.tsv",
    sep="\t", index=False
)

fig, ax = plt.subplots(figsize=(8.5, 5.2))
y = np.arange(len(df))

for i, row in df.iterrows():
    color = "tab:red" if row["significant"] else "0.45"
    ax.errorbar(
        row["odds_ratio_A_vs_B"], i,
        xerr=[[
            row["odds_ratio_A_vs_B"] - row["OR_CI95_lower"]
        ], [
            row["OR_CI95_upper"] - row["odds_ratio_A_vs_B"]
        ]],
        fmt="o", color=color, ecolor=color,
        capsize=4, markersize=6, linewidth=1.6
    )

    ax.annotate(
        f'OR {row["odds_ratio_A_vs_B"]:.2f} '
        f'({row["OR_CI95_lower"]:.2f}–{row["OR_CI95_upper"]:.2f})',
        (row["OR_CI95_upper"], i),
        xytext=(6, 0), textcoords="offset points",
        va="center", fontsize=8
    )

ax.axvline(1, color="black", linestyle="--", linewidth=1)
ax.set_xscale("log")
ax.set_yticks(y)
ax.set_yticklabels(df["label"])
ax.set_xlabel("Adjusted odds ratio relative to Z. mobilis ZM4 (log scale)")
ax.set_title("Adjusted species effects on candidate motif prevalence")
ax.grid(axis="x", linestyle=":", alpha=0.4)

from matplotlib.lines import Line2D
legend_items = [
    Line2D([0], [0], marker="o", color="tab:red", linestyle="None",
           label="FDR < 0.05"),
    Line2D([0], [0], marker="o", color="0.45", linestyle="None",
           label="FDR ≥ 0.05"),
]
ax.legend(handles=legend_items, frameon=False, loc="lower right")

fig.tight_layout()
fig.savefig(OUT / "Figure6_species_vs_ZM4.pdf", bbox_inches="tight")
fig.savefig(
    OUT / "Figure6_species_vs_ZM4.png",
    dpi=600, bbox_inches="tight"
)
plt.close(fig)

print(df[[
    "label", "odds_ratio_A_vs_B", "OR_CI95_lower",
    "OR_CI95_upper", "Wald_p_value",
    "BH_FDR_across_30_tests", "significant"
]].to_string(index=False))

print("\nSaved Figure 6 PDF, PNG and TSV to:", OUT)
