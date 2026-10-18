from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "adjusted_gee_all_30_pairwise.tsv"
OUT = Path("14_expanded_study/results/gee_publication")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT, sep="\t")

required = [
    "model", "species_A", "species_B", "odds_ratio_A_vs_B",
    "OR_CI95_lower", "OR_CI95_upper", "Wald_p_value",
    "BH_FDR_across_30_tests"
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise SystemExit(f"Missing columns in input file: {missing}")

labels = {
    "CT_ATCC27405": "A. thermocellus (CT)",
    "NOVO": "N. aromaticivorans",
    "RSP": "C. sphaeroides",
    "TSA_JW_SL_YS485": "T. saccharolyticum",
    "ZM4": "Z. mobilis ZM4",
    "ZM_NCIMB11163": "Z. mobilis NCIMB 11163",
}

def prepare(frame):
    frame = frame.copy()
    frame["comparison"] = (
        frame["species_A"].map(labels).fillna(frame["species_A"])
        + " vs "
        + frame["species_B"].map(labels).fillna(frame["species_B"])
    )
    frame["significant_FDR_005"] = (
        frame["BH_FDR_across_30_tests"] < 0.05
    )
    return frame.sort_values(
        ["BH_FDR_across_30_tests", "Wald_p_value"]
    )

primary = prepare(df[df["model"] == "all_windows_adjusted"])
sensitivity = prepare(df[df["model"] == "nonoverlap_adjusted"])

if len(primary) != 15 or len(sensitivity) != 15:
    raise SystemExit(
        f"Expected 15 comparisons per model; got "
        f"{len(primary)} and {len(sensitivity)}."
    )

# Export complete model-specific tables.
primary.to_csv(OUT / "Table_primary_all_windows.tsv", sep="\t", index=False)
sensitivity.to_csv(
    OUT / "Table_sensitivity_nonoverlap.tsv", sep="\t", index=False
)

significant = pd.concat([primary, sensitivity], ignore_index=True)
significant = significant[significant["significant_FDR_005"]]
significant.to_csv(
    OUT / "Table_FDR_significant_comparisons.tsv",
    sep="\t", index=False
)

def forest_plot(data, title, filename):
    d = data.sort_values("BH_FDR_across_30_tests").copy()
    d = d.iloc[::-1].reset_index(drop=True)

    y = np.arange(len(d))
    OR = d["odds_ratio_A_vs_B"].to_numpy(float)
    lo = d["OR_CI95_lower"].to_numpy(float)
    hi = d["OR_CI95_upper"].to_numpy(float)
    fdr = d["BH_FDR_across_30_tests"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    for i in range(len(d)):
        ax.plot([lo[i], hi[i]], [y[i], y[i]], linewidth=1.5)
        ax.plot(OR[i], y[i], marker="o", markersize=5)

    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(d["comparison"], fontsize=8.5)
    ax.set_xlabel("Adjusted odds ratio (95% confidence interval)")
    ax.set_title(title)
    ax.grid(axis="x", linestyle=":", alpha=0.45)
    ax.tick_params(axis="y", length=0)

    # Mark FDR-significant comparisons in the right margin.
    xmax = max(hi.max(), 1.1)
    xmin = max(min(lo.min(), 0.5), 0.01)
    ax.set_xlim(xmin * 0.8, xmax * 2.5)

    for yi, q in zip(y, fdr):
        marker = "FDR < 0.05" if q < 0.05 else ""
        if marker:
            ax.text(
                xmax * 1.05, yi, marker,
                va="center", fontsize=7
            )

    fig.tight_layout()
    fig.savefig(OUT / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{filename}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

forest_plot(
    primary,
    "Adjusted pairwise motif-prevalence comparisons: all upstream windows",
    "Figure_GEE_primary_all_windows"
)
forest_plot(
    sensitivity,
    "Sensitivity analysis: non-overlapping upstream windows",
    "Figure_GEE_sensitivity_nonoverlap"
)

print("\nPublication outputs saved in:", OUT)
print("\nPrimary model: FDR-significant comparisons")
print(primary.loc[primary["significant_FDR_005"], [
    "comparison", "odds_ratio_A_vs_B", "OR_CI95_lower",
    "OR_CI95_upper", "Wald_p_value", "BH_FDR_across_30_tests"
]].to_string(index=False))

print("\nSensitivity model: FDR-significant comparisons")
print(sensitivity.loc[sensitivity["significant_FDR_005"], [
    "comparison", "odds_ratio_A_vs_B", "OR_CI95_lower",
    "OR_CI95_upper", "Wald_p_value", "BH_FDR_across_30_tests"
]].to_string(index=False))

print("\nFiles created:")
for p in sorted(OUT.iterdir()):
    print(p)
