#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

base = Path("14_expanded_study/orthology/orthofinder_input/OrthoFinder")
results = base / "Results_six_genome_v1"
orthogroups_file = results / "Orthogroups/Orthogroups.tsv"
single_copy_file = results / "Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

outdir = Path("14_expanded_study/orthology/single_copy_audit")
outdir.mkdir(parents=True, exist_ok=True)

with single_copy_file.open() as f:
    single_copy_ids = {
        line.strip() for line in f
        if line.strip() and not line.startswith("#")
    }

with orthogroups_file.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames
    species = [x for x in fields if x != "Orthogroup"]
    rows = list(reader)

selected = []
problems = []

for row in rows:
    og = row["Orthogroup"]
    if og not in single_copy_ids:
        continue

    counts = {
        sp: len([x.strip() for x in (row.get(sp) or "").split(",") if x.strip()])
        for sp in species
    }

    if all(counts[sp] == 1 for sp in species):
        selected.append(row)
    else:
        problems.append((og, counts))

outfile = outdir / "verified_single_copy_orthogroups.tsv"
with outfile.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(selected)

print("Species:", ", ".join(species))
print("OrthoFinder single-copy IDs:", len(single_copy_ids))
print("Verified one protein per species:", len(selected))
print("Unexpected membership patterns:", len(problems))
print("Output:", outfile)

if problems:
    print("\nFirst unexpected patterns:")
    for og, counts in problems[:10]:
        print(og, counts)
