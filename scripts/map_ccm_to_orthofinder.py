#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("14_expanded_study")
RESULTS = (
    BASE / "orthology/orthofinder_input/OrthoFinder/"
    "Results_six_genome_v1"
)
OG_FILE = RESULTS / "Orthogroups/Orthogroups.tsv"
CCM_FILE = Path("07_harmonized_analysis/ccm_strict_broad_gene_table.tsv")
ZM4_FILE = Path("07_harmonized_analysis/zm4_harmonized_ccm.tsv")
INPUT = BASE / "orthology/orthofinder_input"
OUTDIR = BASE / "orthology/audit"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "strict_ZM4_CCM_OrthoFinder_mapping.tsv"
UNMAPPED = OUTDIR / "strict_ZM4_CCM_unmapped.tsv"

for path in (OG_FILE, CCM_FILE, ZM4_FILE):
    if not path.exists():
        raise SystemExit(f"ERROR: required file not found: {path}")

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def proteins(cell):
    return [x.strip() for x in (cell or "").split(",") if x.strip()]

# Read FASTA headers to obtain product descriptions for mapped proteins.
descriptions = {}
for faa in INPUT.glob("*.faa"):
    with faa.open() as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                accession = header.split()[0]
                descriptions[(faa.stem, accession)] = header

# Map each protein accession to its orthogroup.
og_rows = read_tsv(OG_FILE)
species = [c for c in og_rows[0] if c != "Orthogroup"]
protein_to_og = {}
og_to_species = {}

for row in og_rows:
    og = row["Orthogroup"]
    og_to_species[og] = {
        sp: proteins(row.get(sp, "")) for sp in species
    }
    for sp in species:
        for acc in og_to_species[og][sp]:
            protein_to_og.setdefault((sp, acc), []).append(og)

# Preserve the existing harmonized scores and classifications.
harmonized = {
    r["protein_id"]: r
    for r in read_tsv(ZM4_FILE)
    if r.get("protein_id")
}

# Only genes already classified STRICT in the existing source table.
ccm_rows = [
    r for r in read_tsv(CCM_FILE)
    if r.get("organism") == "ZM4" and r.get("status") == "STRICT"
]

mapped_rows = []
unmapped_rows = []

for ccm in ccm_rows:
    acc = ccm.get("protein_id", "").strip()
    locus = ccm.get("locus_tag", "")
    h = harmonized.get(acc, {})

    groups = protein_to_og.get(("ZM4", acc), [])
    if not groups:
        unmapped_rows.append({
            "organism": "ZM4",
            "locus_tag": locus,
            "gene": ccm.get("gene", ""),
            "protein_id": acc,
            "product": ccm.get("product", ""),
            "status": "ZM4_PROTEIN_NOT_FOUND_IN_ORTHOGROUP_TABLE"
        })
        continue

    for og in groups:
        members = og_to_species[og]
        counts = {sp: len(members[sp]) for sp in species}
        all_species = all(counts[sp] > 0 for sp in species)
        single_copy = all(counts[sp] == 1 for sp in species)

        row = {
            "ZM4_locus_tag": locus,
            "ZM4_gene": ccm.get("gene", ""),
            "ZM4_protein_id": acc,
            "ZM4_product": ccm.get("product", ""),
            "original_CCM_status": ccm.get("status", ""),
            "original_evidence": ccm.get("evidence", ""),
            "original_rank": h.get("rank", ""),
            "original_integrated_score": h.get("integrated_score", ""),
            "original_harmonized_CCM_class": h.get("harmonized_ccm", ""),
            "orthogroup": og,
            "species_present": sum(counts[sp] > 0 for sp in species),
            "all_six_species_present": all_species,
            "single_copy_in_all_six": single_copy,
        }

        for sp in species:
            ids = members[sp]
            row[f"{sp}_protein_count"] = len(ids)
            row[f"{sp}_protein_ids"] = ";".join(ids)
            row[f"{sp}_products"] = " || ".join(
                descriptions.get((sp, pid), pid) for pid in ids
            )

        if len(groups) > 1:
            row["mapping_note"] = "Protein accession associated with multiple orthogroups; investigate"
        elif single_copy:
            row["mapping_note"] = "Complete single-copy orthogroup; still verify gene function"
        elif all_species:
            row["mapping_note"] = "All species present, but at least one species has multiple proteins"
        else:
            row["mapping_note"] = "At least one species absent from this orthogroup"

        mapped_rows.append(row)

def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

write_tsv(OUTFILE, mapped_rows)
write_tsv(UNMAPPED, unmapped_rows)

print(f"Strict ZM4 CCM entries in source: {len(ccm_rows)}")
print(f"Mapping rows produced: {len(mapped_rows)}")
print(f"Unmapped CCM entries: {len(unmapped_rows)}")
print(f"Mapping table: {OUTFILE}")
print(f"Unmapped table: {UNMAPPED}")

if mapped_rows:
    print("\nMapping summary:")
    print(f"  Unique ZM4 proteins mapped: {len(set(r['ZM4_protein_id'] for r in mapped_rows))}")
    print(f"  Complete six-species groups: {sum(r['all_six_species_present'] for r in mapped_rows)}")
    print(f"  Complete single-copy groups: {sum(r['single_copy_in_all_six'] for r in mapped_rows)}")
