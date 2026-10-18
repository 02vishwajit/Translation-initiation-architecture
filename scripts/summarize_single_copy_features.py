#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

src = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_sequence_features.tsv"
)
outdir = src.parent

with src.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

# Four confirmed single-copy CCM orthogroups.
ccm_ogs = {
    "OG0000678",  # phosphoglycerate kinase
    "OG0000449",  # triose-phosphate isomerase
    "OG0000462",  # glucose-6-phosphate isomerase
    "OG0000450",  # enolase
}

def summarize(group):
    n = len(group)
    positives = sum(
        r["has_candidate_spacer_4_12"].lower() == "true"
        for r in group
    )
    return {
        "n": n,
        "motif_positive": positives,
        "motif_percent": round(100 * positives / n, 2) if n else "",
        "mean_GC_percent": round(
            sum(float(r["GC_percent"]) for r in group) / n, 3
        ) if n else "",
        "mean_purine_percent": round(
            sum(float(r["purine_percent"]) for r in group) / n, 3
        ) if n else "",
    }

summaries = []

# Species-wide summary, stratified by CDS overlap.
for sp in sorted({r["species"] for r in rows}):
    for overlap in ("False", "True"):
        group = [
            r for r in rows
            if r["species"] == sp
            and r["overlap_flag"].lower() == overlap.lower()
        ]
        if not group:
            continue
        s = summarize(group)
        summaries.append({
            "analysis": "species_by_overlap",
            "group": sp,
            "subset": "overlap" if overlap == "True" else "no_overlap",
            **s,
        })

# CCM versus non-CCM background, separately by species and overlap.
for sp in sorted({r["species"] for r in rows}):
    for label, predicate in (
        ("four_single_copy_CCM", lambda r: r["orthogroup"] in ccm_ogs),
        ("other_single_copy_orthogroups", lambda r: r["orthogroup"] not in ccm_ogs),
    ):
        for overlap in ("False", "True"):
            group = [
                r for r in rows
                if r["species"] == sp
                and predicate(r)
                and r["overlap_flag"].lower() == overlap.lower()
            ]
            if not group:
                continue
            s = summarize(group)
            summaries.append({
                "analysis": "CCM_vs_background",
                "group": sp,
                "subset": label + (
                    "_overlap" if overlap == "True" else "_no_overlap"
                ),
                **s,
            })

outfile = outdir / "single_copy_feature_summary.tsv"
with outfile.open("w", newline="") as f:
    fields = [
        "analysis", "group", "subset", "n", "motif_positive",
        "motif_percent", "mean_GC_percent", "mean_purine_percent",
    ]
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summaries)

print("Output:", outfile)
print("\nSpecies-wide summary by overlap status:")
for r in summaries:
    if r["analysis"] == "species_by_overlap":
        print(
            f'{r["group"]:18s} {r["subset"]:12s} '
            f'n={r["n"]:3d} motif={r["motif_percent"]:6.2f}% '
            f'GC={r["mean_GC_percent"]:6.2f}% '
            f'purine={r["mean_purine_percent"]:6.2f}%'
        )

print("\nFour single-copy CCM orthogroups:")
for r in summaries:
    if (
        r["analysis"] == "CCM_vs_background"
        and r["subset"].startswith("four_single_copy_CCM")
    ):
        print(
            f'{r["group"]:18s} {r["subset"]:35s} '
            f'n={r["n"]:2d} motif={r["motif_percent"]:6.2f}%'
        )
