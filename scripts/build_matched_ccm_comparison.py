#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("14_expanded_study")
FEATURES = BASE / "features/single_copy_upstream/single_copy_sequence_features_corrected.tsv"
CCM_MAP = BASE / "orthology/audit/strict_ZM4_CCM_OrthoFinder_mapping.tsv"
OUTDIR = BASE / "results/matched_ccm"

OUTDIR.mkdir(parents=True, exist_ok=True)

TARGET_CCM = {
    "OG0000678": "phosphoglycerate kinase",
    "OG0000449": "triose-phosphate isomerase",
    "OG0000462": "glucose-6-phosphate isomerase",
    "OG0000450": "enolase",
}

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

features = read_tsv(FEATURES)
ccm_mapping = read_tsv(CCM_MAP)

# Validate the four prespecified CCM orthogroups against the mapping table.
mapped_targets = {
    r["orthogroup"]: r
    for r in ccm_mapping
    if r["orthogroup"] in TARGET_CCM
    and r["single_copy_in_all_six"].lower() == "true"
}

if set(mapped_targets) != set(TARGET_CCM):
    raise ValueError(
        "CCM mapping validation failed. "
        f"Expected {sorted(TARGET_CCM)}, found {sorted(mapped_targets)}"
    )

# Index features and check one record per orthogroup/species.
feature_map = {}
by_species = defaultdict(list)

for row in features:
    key = (row["species"], row["orthogroup"])
    if key in feature_map:
        raise ValueError(f"Duplicate feature record: {key}")
    feature_map[key] = row
    by_species[row["species"]].append(row)

species = sorted(by_species)
orthogroups = sorted({r["orthogroup"] for r in features})

if len(orthogroups) != 348:
    raise ValueError(f"Expected 348 orthogroups; found {len(orthogroups)}")

for sp in species:
    if len(by_species[sp]) != 348:
        raise ValueError(f"{sp}: expected 348 records")

for og in orthogroups:
    for sp in species:
        if (sp, og) not in feature_map:
            raise ValueError(f"Missing matched record: {sp}, {og}")

# Write the four-gene CCM table.
ccm_rows = []

for og, description in TARGET_CCM.items():
    for sp in species:
        r = feature_map[(sp, og)]
        ccm_rows.append({
            "orthogroup": og,
            "CCM_gene": description,
            "species": sp,
            "protein_id": r["protein_id"],
            "locus_tag": r["locus_tag"],
            "product": r["product"],
            "start_codon": r["start_codon"],
            "motif_positive": int(
                r["has_candidate_spacer_4_12"].lower() == "true"
            ),
            "best_candidate_motif": r["best_candidate_motif"],
            "best_candidate_spacer_nt": r["best_candidate_spacer_nt"],
            "GC_percent": float(r["GC_percent"]),
            "purine_percent": float(r["purine_percent"]),
            "overlap_flag": r["overlap_flag"],
        })

with (OUTDIR / "four_CCM_genes_by_species.tsv").open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(ccm_rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(ccm_rows)

# Summarize CCM versus non-CCM within each species.
summary = []

for sp in species:
    for group_name, group_ogs in [
        ("Four_single_copy_CCM", set(TARGET_CCM)),
        ("Other_single_copy_orthogroups",
         set(orthogroups) - set(TARGET_CCM)),
    ]:
        rows = [
            feature_map[(sp, og)]
            for og in sorted(group_ogs)
        ]

        n = len(rows)
        positives = sum(
            r["has_candidate_spacer_4_12"].lower() == "true"
            for r in rows
        )

        overlaps = [
            r for r in rows
            if r["overlap_flag"].lower() == "true"
        ]

        summary.append({
            "species": sp,
            "group": group_name,
            "n_orthogroups": n,
            "motif_positive": positives,
            "motif_positive_percent": round(100 * positives / n, 3),
            "mean_GC_percent": round(
                sum(float(r["GC_percent"]) for r in rows) / n, 3
            ),
            "overlapping_windows": len(overlaps),
            "overlapping_window_percent": round(
                100 * len(overlaps) / n, 3
            ),
        })

with (OUTDIR / "CCM_vs_other_single_copy_summary.tsv").open(
    "w", newline=""
) as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)

# Produce an orthogroup-level matched table across all six species.
matched = []

for og in orthogroups:
    group = "Four_single_copy_CCM" if og in TARGET_CCM else "Other_single_copy"
    row = {
        "orthogroup": og,
        "group": group,
        "CCM_description": TARGET_CCM.get(og, ""),
    }

    for sp in species:
        r = feature_map[(sp, og)]
        row[f"{sp}_motif_positive"] = int(
            r["has_candidate_spacer_4_12"].lower() == "true"
        )
        row[f"{sp}_GC_percent"] = float(r["GC_percent"])
        row[f"{sp}_overlap"] = r["overlap_flag"]
        row[f"{sp}_start_codon"] = r["start_codon"]

    matched.append(row)

with (OUTDIR / "matched_orthogroup_feature_matrix.tsv").open(
    "w", newline=""
) as f:
    writer = csv.DictWriter(
        f, fieldnames=list(matched[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(matched)

print("Matched CCM comparison completed.")
print("Species:", ", ".join(species))
print("Matched orthogroups:", len(orthogroups))
print("CCM orthogroups:", len(TARGET_CCM))
print("Other single-copy orthogroups:", len(orthogroups) - len(TARGET_CCM))
print("\nSpecies-level descriptive summary:")

for r in summary:
    print(
        f'{r["species"]:18s} {r["group"]:30s} '
        f'n={r["n_orthogroups"]:3d} '
        f'motif={r["motif_positive_percent"]:6.2f}% '
        f'mean_GC={r["mean_GC_percent"]:6.2f}% '
        f'overlap={r["overlapping_window_percent"]:6.2f}%'
    )

print("\nOutputs:")
print(OUTDIR / "four_CCM_genes_by_species.tsv")
print(OUTDIR / "CCM_vs_other_single_copy_summary.tsv")
print(OUTDIR / "matched_orthogroup_feature_matrix.tsv")
print("\nInterpretation: descriptive comparison only; four CCM genes")
print("are insufficient for a strong general pathway-level claim.")
