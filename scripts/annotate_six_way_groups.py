from pathlib import Path
import csv
import re

root = Path("14_expanded_study")
matrix = root / "orthology/six_genome_matrix/strict_six_way_ortholog_matrix.tsv"
outdir = root / "orthology/six_genome_matrix"
strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

def read_fasta(path):
    annotations = {}
    if not path.exists():
        print(f"WARNING: missing FASTA: {path}")
        return annotations

    with path.open() as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            accession = header.split()[0]
            # NCBI protein headers commonly use:
            # accession [protein product] [organism]
            product = re.search(r"\[([^\]]+)\]", header)
            annotations[accession] = (
                product.group(1) if product else header
            )
    return annotations

annotation = {}
for strain in strains:
    faa = root / "genomes/standardized" / strain / f"{strain}.faa"
    annotation[strain] = read_fasta(faa)

with matrix.open(newline="") as f:
    groups = list(csv.DictReader(f, delimiter="\t"))

annotated_rows = []
for group in groups:
    row = {"ortholog_group": group["ortholog_group"]}
    for strain in strains:
        accession = group[strain]
        row[f"{strain}_protein"] = accession
        row[f"{strain}_product"] = annotation[strain].get(
            accession, "ANNOTATION_NOT_FOUND"
        )
    annotated_rows.append(row)

out = outdir / "strict_six_way_orthologs_annotated.tsv"
with out.open("w", newline="") as f:
    fields = ["ortholog_group"]
    for strain in strains:
        fields += [f"{strain}_protein", f"{strain}_product"]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(annotated_rows)

print(f"Annotated matrix: {out}")
print(f"Groups: {len(annotated_rows)}")
for strain in strains:
    missing = sum(
        row[f"{strain}_product"] == "ANNOTATION_NOT_FOUND"
        for row in annotated_rows
    )
    print(f"{strain}: {missing} protein annotations not found")

print("\nFirst 10 groups and product annotations:")
for row in annotated_rows[:10]:
    print("\n" + row["ortholog_group"])
    for strain in strains:
        print(f"  {strain}: {row[f'{strain}_product']}")
