#!/usr/bin/env python3

import csv
from pathlib import Path

base = Path("14_expanded_study/features/single_copy_upstream")

features_file = base / "single_copy_sequence_features.tsv"
null_file = base / "single_copy_dinucleotide_null.tsv"
output_file = base / "single_copy_dinucleotide_null_repaired.tsv"

with features_file.open(newline="") as f:
    feature_rows = list(csv.DictReader(f, delimiter="\t"))

with null_file.open(newline="") as f:
    null_rows = list(csv.DictReader(f, delimiter="\t"))

feature_map = {}
for row in feature_rows:
    key = (row["species"], row["orthogroup"])
    if key in feature_map:
        raise ValueError(f"Duplicate feature key: {key}")
    feature_map[key] = row

repaired = []
seen = set()

for row in null_rows:
    key = (row["species"], row["orthogroup"])
    if key not in feature_map:
        raise ValueError(f"No matching feature row: {key}")
    if key in seen:
        raise ValueError(f"Duplicate null row: {key}")
    seen.add(key)

    feat = feature_map[key]

    # Confirm the stored observed result agrees with the original feature scan.
    expected = int(feat["has_candidate_spacer_4_12"].lower() == "true")
    if int(row["observed_candidate"]) != expected:
        raise ValueError(
            f"Observed motif mismatch for {key}: "
            f"null={row['observed_candidate']}, feature={expected}"
        )

    row["sequence_id"] = feat["protein_id"]
    row["gc_fraction"] = f"{float(feat['GC_percent']) / 100:.6f}"
    row["overlaps_cds"] = feat["overlap_flag"]
    repaired.append(row)

if len(seen) != len(feature_map):
    missing = set(feature_map) - seen
    raise ValueError(f"Missing null rows: {len(missing)}")

with output_file.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(repaired[0].keys()), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(repaired)

print("Metadata repair and matching validation passed.")
print(f"Feature rows: {len(feature_rows)}")
print(f"Null rows: {len(null_rows)}")
print(f"Matched rows: {len(repaired)}")
print(f"Output: {output_file}")
print("Note: sequence_id is populated with protein_id.")
