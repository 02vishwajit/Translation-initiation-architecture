#!/usr/bin/env python3

import csv
import re
from pathlib import Path
from collections import Counter

BASE = Path("14_expanded_study")
INPUT = BASE / "orthology/orthofinder_input"
RESULTS = INPUT / "OrthoFinder/Results_six_genome_v1"
OG_FILE = RESULTS / "Orthogroups/Orthogroups.tsv"
SINGLE_FILE = RESULTS / "Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
OUT = BASE / "orthology/audit"
OUT.mkdir(parents=True, exist_ok=True)

if not OG_FILE.exists():
    raise SystemExit(f"ERROR: missing {OG_FILE}")

# Map protein identifiers to FASTA-header descriptions.
descriptions = {}
for faa in sorted(INPUT.glob("*.faa")):
    with faa.open() as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                accession = header.split()[0]
                descriptions[(faa.stem, accession)] = header

def proteins(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

with OG_FILE.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    species = [x for x in reader.fieldnames if x != "Orthogroup"]
    rows = list(reader)

single_copy_ids = set()
if SINGLE_FILE.exists():
    with SINGLE_FILE.open() as handle:
        single_copy_ids = {
            line.strip() for line in handle if line.strip()
        }

audit_rows = []
candidate_rows = []
all_species_groups = 0
single_copy_all_species = 0
assigned_gene_count = 0
groups_by_species_count = Counter()

# Broad candidate terms. These flag records for review; they do not
# establish pathway membership or functional activity.
pattern = re.compile(
    r"glycol|glucon|glucose|fructose|mannose|xylose|arabinose|"
    r"pentose|phosphoglucon|transketolase|transaldolase|"
    r"pyruvate|phosphofruct|fructose.?bisphosphate|"
    r"glyceraldehyde.?3.?phosphate|triose.?phosphate|"
    r"phosphoglycerate|enolase|phosphopyruvate|"
    r"lactate|ethanol|alcohol dehydrogenase|"
    r"pyruvate decarboxylase|acetaldehyde|acetyl.?coa|"
    r"phosphotransferase system|sugar transporter|"
    r"hexokinase|glucokinase|phosphoglucomutase|"
    r"phosphogluconate|ribulose.?phosphate|"
    r"ribose.?phosphate isomerase|"
    r"ferredoxin.?dependent oxidoreductase",
    re.IGNORECASE
)

for row in rows:
    og = row["Orthogroup"]
    by_species = {sp: proteins(row.get(sp, "")) for sp in species}
    present = sum(bool(v) for v in by_species.values())
    total = sum(len(v) for v in by_species.values())
    assigned_gene_count += total
    groups_by_species_count[present] += 1

    complete = present == len(species)
    single_copy = complete and all(len(v) == 1 for v in by_species.values())

    if complete:
        all_species_groups += 1
    if single_copy:
        single_copy_all_species += 1

    audit_rows.append({
        "orthogroup": og,
        "species_present": present,
        "species_total": len(species),
        "genes_total": total,
        "all_species_present": complete,
        "single_copy_all_species": single_copy,
        "in_orthofinder_single_copy_list": og in single_copy_ids,
        **{f"{sp}_gene_count": len(by_species[sp]) for sp in species}
    })

    hits = []
    for sp, ids in by_species.items():
        for accession in ids:
            desc = descriptions.get((sp, accession), accession)
            if pattern.search(desc):
                hits.append((sp, accession, desc))

    if hits:
        candidate_rows.append({
            "orthogroup": og,
            "species_present": present,
            "genes_total": total,
            "all_species_present": complete,
            "single_copy_all_species": single_copy,
            "candidate_protein_hits": len(hits),
            "candidate_hits": " || ".join(
                f"{sp}:{acc}:{desc}" for sp, acc, desc in hits
            )
        })

def write_tsv(path, records):
    if not records:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(records[0].keys()),
            delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(records)

write_tsv(OUT / "orthogroup_audit.tsv", audit_rows)
write_tsv(OUT / "carbon_metabolism_candidate_groups.tsv", candidate_rows)

summary = OUT / "orthofinder_audit_summary.txt"
with summary.open("w") as handle:
    handle.write(f"Species: {', '.join(species)}\n")
    handle.write(f"Orthogroups in table: {len(rows)}\n")
    handle.write(f"Genes represented in orthogroup table: {assigned_gene_count}\n")
    handle.write(f"Groups with all species present: {all_species_groups}\n")
    handle.write(f"Complete single-copy groups (computed from table): {single_copy_all_species}\n")
    handle.write(f"Groups in OrthoFinder single-copy list: {len(single_copy_ids)}\n")
    handle.write(f"Groups flagged by broad annotation keywords: {len(candidate_rows)}\n")
    handle.write("\nOrthogroups by number of species represented:\n")
    for n in sorted(groups_by_species_count):
        handle.write(f"{n}/{len(species)} species: {groups_by_species_count[n]}\n")

print(summary.read_text())
print(f"Audit table: {OUT / 'orthogroup_audit.tsv'}")
print(f"Candidate table: {OUT / 'carbon_metabolism_candidate_groups.tsv'}")
