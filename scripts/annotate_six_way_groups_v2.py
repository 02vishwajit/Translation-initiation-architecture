from pathlib import Path
from collections import defaultdict
import csv
import re

root = Path("14_expanded_study")
outdir = root / "orthology/six_genome_matrix"
matrix = outdir / "strict_six_way_ortholog_matrix.tsv"

strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

# accession -> annotation fields
annotations = {s: {} for s in strains}

# First: use GenBank CDS feature qualifiers where available.
try:
    from Bio import SeqIO

    for strain in strains:
        gbff = root / "genomes/standardized" / strain / f"{strain}.gbff"
        if not gbff.exists():
            print(f"WARNING: missing GenBank file: {gbff}")
            continue

        for record in SeqIO.parse(str(gbff), "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue

                q = feature.qualifiers
                protein_ids = q.get("protein_id", [])
                if not protein_ids:
                    continue

                product = q.get("product", [""])[0]
                gene = q.get("gene", [""])[0]
                locus = q.get("locus_tag", [""])[0]

                for protein_id in protein_ids:
                    annotations[strain][protein_id] = {
                        "product": product or "product_not_annotated",
                        "gene": gene or "",
                        "locus_tag": locus or "",
                        "source": "GenBank_CDS"
                    }

except ImportError:
    print("WARNING: Biopython unavailable; using FASTA-header fallback only.")

# Fallback: parse protein description BEFORE the organism in brackets.
for strain in strains:
    faa = root / "genomes/standardized" / strain / f"{strain}.faa"
    if not faa.exists():
        continue

    with faa.open() as f:
        for line in f:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            accession = header.split()[0]

            if accession in annotations[strain]:
                continue

            description = header[len(accession):].strip()
            # NCBI protein FASTA commonly ends with [organism].
            description = re.sub(r"\s*\[[^\]]*\]\s*$", "", description).strip()
            description = re.sub(r"\s+\(.*?\)\s*$", "", description).strip()

            annotations[strain][accession] = {
                "product": description or "product_not_annotated",
                "gene": "",
                "locus_tag": "",
                "source": "FASTA_header"
            }

with matrix.open(newline="") as f:
    groups = list(csv.DictReader(f, delimiter="\t"))

out = outdir / "strict_six_way_orthologs_annotated_v2.tsv"
fields = ["ortholog_group"]

for strain in strains:
    fields.extend([
        f"{strain}_protein",
        f"{strain}_gene",
        f"{strain}_product",
        f"{strain}_locus_tag",
        f"{strain}_annotation_source"
    ])

with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for group in groups:
        row = {"ortholog_group": group["ortholog_group"]}

        for strain in strains:
            accession = group[strain]
            ann = annotations[strain].get(accession, {})

            row[f"{strain}_protein"] = accession
            row[f"{strain}_gene"] = ann.get("gene", "")
            row[f"{strain}_product"] = ann.get(
                "product", "ANNOTATION_NOT_FOUND"
            )
            row[f"{strain}_locus_tag"] = ann.get("locus_tag", "")
            row[f"{strain}_annotation_source"] = ann.get("source", "")

        writer.writerow(row)

print(f"\nCorrected annotation table: {out}")
print(f"Ortholog groups: {len(groups)}")

for strain in strains:
    missing = sum(
        annotations[strain].get(group[strain], {}).get(
            "product", "ANNOTATION_NOT_FOUND"
        ) in ("", "ANNOTATION_NOT_FOUND", "product_not_annotated")
        for group in groups
    )
    print(f"{strain}: {missing}/{len(groups)} products unavailable")
