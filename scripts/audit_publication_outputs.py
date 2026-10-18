from pathlib import Path
import pandas as pd
import numpy as np

base = Path("14_expanded_study")
features = base / "features/single_copy_upstream"
out = base / "results/gee_publication"

def read_tsv(path):
    if not path.exists():
        print(f"MISSING FILE: {path}")
        return None
    return pd.read_csv(path, sep="\t")

print("=" * 72)
print("AUDIT 1: SOURCE MOTIF COUNTS")
print("=" * 72)

source = read_tsv(features / "single_copy_sequence_features_corrected.tsv")
if source is not None:
    source["positive"] = (
        source["has_candidate_spacer_4_12"]
        .astype(str).str.lower().map({"true": 1, "false": 0})
    )
    counts = source.groupby("species")["positive"].agg(["sum", "count"])
    counts["prevalence_percent"] = 100 * counts["sum"] / counts["count"]
    print(counts.to_string())

    for name in [
        "Figure3_species_motif_prevalence.tsv",
        "Figure8_motif_prevalence.tsv",
    ]:
        p = out / name
        table = read_tsv(p)
        print(f"\n{name}:")
        if table is not None:
            print("Columns:", table.columns.tolist())
            print(table.to_string(index=False))

print("\n" + "=" * 72)
print("AUDIT 2: GEE SOURCE TABLE")
print("=" * 72)

gee = read_tsv(features / "adjusted_gee_all_30_pairwise.tsv")
if gee is not None:
    print("Rows:", len(gee))
    print("Models:", gee["model"].value_counts().to_dict())
    print("\nSpecies-vs-ZM4 rows:")
    mask = (
        (gee["model"] == "all_windows_adjusted")
        & (
            ((gee["species_A"] != "ZM4") & (gee["species_B"] == "ZM4"))
            | ((gee["species_A"] == "ZM4") & (gee["species_B"] != "ZM4"))
        )
    )
    cols = [
        "species_A", "species_B", "odds_ratio_A_vs_B",
        "OR_CI95_lower", "OR_CI95_upper",
        "Wald_p_value", "BH_FDR_across_30_tests"
    ]
    print(gee.loc[mask, cols].to_string(index=False))

print("\n" + "=" * 72)
print("AUDIT 3: FIGURE 8 EXPORTED VALUES")
print("=" * 72)

for name in ["Figure8_adjusted_GEE.tsv", "Figure6_species_vs_ZM4.tsv"]:
    table = read_tsv(out / name)
    print(f"\n{name}:")
    if table is not None:
        print("Rows:", len(table))
        print("Columns:", table.columns.tolist())
        print(table.to_string(index=False))

print("\n" + "=" * 72)
print("AUDIT 4: DINUCLEOTIDE NULL SUMMARY")
print("=" * 72)

null = read_tsv(features / "single_copy_ushuffle_species_summary.tsv")
if null is not None:
    cols = [
        c for c in [
            "species", "n_sequences", "observed_positive",
            "observed_prevalence_percent", "null_mean_prevalence_percent",
            "observed_minus_null_percentage_points",
            "empirical_upper_tail_p", "replicates", "BH_FDR"
        ] if c in null.columns
    ]
    print(null[cols].to_string(index=False))

print("\nAudit report generated in terminal only; no source files changed.")
