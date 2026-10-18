#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter
from Bio import SeqIO

BASE = Path("14_expanded_study")
GENOMES = BASE / "genomes/standardized"
OG_FILE = (
    BASE / "orthology/single_copy_audit/"
    "verified_single_copy_orthogroups.tsv"
)
OUT = BASE / "features/single_copy_upstream"
OUT.mkdir(parents=True, exist_ok=True)

WINDOW = 50

with OG_FILE.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    og_rows = list(reader)
    species = [x for x in reader.fieldnames if x != "Orthogroup"]

# Map each expected protein to its orthogroup.
wanted = {}
for row in og_rows:
    og = row["Orthogroup"]
    for sp in species:
        proteins = [
            x.strip() for x in (row.get(sp) or "").split(",")
            if x.strip()
        ]
        for pid in proteins:
            wanted[(sp, pid)] = og

records = []
seen = set()

for sp in species:
    gbff = GENOMES / sp / f"{sp}.gbff"

    if not gbff.exists():
        print("WARNING: missing GenBank file:", gbff)
        continue

    for rec in SeqIO.parse(str(gbff), "genbank"):
        cds_features = [
            feat for feat in rec.features if feat.type == "CDS"
        ]

        for feat in cds_features:
            protein_ids = feat.qualifiers.get("protein_id", [])
            if not protein_ids:
                continue

            pid = protein_ids[0].strip()
            key = (sp, pid)

            if key not in wanted:
                continue

            start = int(feat.location.start)
            end = int(feat.location.end)
            strand = feat.location.strand
            contig_len = len(rec.seq)

            # Extract upstream DNA in the gene's 5' -> 3' orientation.
            if strand == 1:
                lo, hi = start - WINDOW, start
                if lo < 0:
                    seq = ""
                    status = "CONTIG_BOUNDARY_INCOMPLETE"
                else:
                    seq = str(rec.seq[lo:hi]).upper()
                    status = "OK"

            elif strand == -1:
                lo, hi = end, end + WINDOW
                if hi > contig_len:
                    seq = ""
                    status = "CONTIG_BOUNDARY_INCOMPLETE"
                else:
                    seq = str(
                        rec.seq[lo:hi].reverse_complement()
                    ).upper()
                    status = "OK"

            else:
                lo, hi = start, end
                seq = ""
                status = "UNKNOWN_STRAND"

            # Flag any other CDS intersecting the upstream genomic interval.
            overlaps = []
            if status == "OK":
                for other in cds_features:
                    if other is feat:
                        continue
                    os = int(other.location.start)
                    oe = int(other.location.end)
                    if os < hi and oe > lo:
                        overlaps.append(
                            other.qualifiers.get(
                                "locus_tag", ["unknown"]
                            )[0]
                        )

            cds_seq = str(feat.extract(rec.seq)).upper()
            record = {
                "orthogroup": wanted[key],
                "species": sp,
                "protein_id": pid,
                "locus_tag": feat.qualifiers.get(
                    "locus_tag", [""]
                )[0],
                "gene": feat.qualifiers.get("gene", [""])[0],
                "product": feat.qualifiers.get("product", [""])[0],
                "contig": rec.id,
                "strand": strand,
                "start_codon": cds_seq[:3],
                "upstream_length_nt": len(seq),
                "upstream_sequence_5to3": seq,
                "status": status,
                "overlapping_CDS": ";".join(sorted(set(overlaps))),
            }
            records.append(record)
            seen.add(key)

# Explicitly record expected proteins that could not be matched to CDS features.
for (sp, pid), og in wanted.items():
    if (sp, pid) not in seen:
        records.append({
            "orthogroup": og,
            "species": sp,
            "protein_id": pid,
            "locus_tag": "",
            "gene": "",
            "product": "",
            "contig": "",
            "strand": "",
            "start_codon": "",
            "upstream_length_nt": 0,
            "upstream_sequence_5to3": "",
            "status": "PROTEIN_ID_NOT_FOUND_IN_GENBANK_CDS",
            "overlapping_CDS": "",
        })

tsv = OUT / "single_copy_upstream50.tsv"
with tsv.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(records[0].keys()), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(records)

fasta = OUT / "single_copy_upstream50.fasta"
with fasta.open("w") as f:
    for r in records:
        if r["status"] != "OK" or not r["upstream_sequence_5to3"]:
            continue
        f.write(
            f'>{r["orthogroup"]}|{r["species"]}|'
            f'{r["protein_id"]}|{r["locus_tag"]}\n'
        )
        f.write(r["upstream_sequence_5to3"] + "\n")

print("Orthogroups:", len(og_rows))
print("Species:", len(species))
print("Expected protein records:", len(wanted))
print("Records written:", len(records))
print("TSV:", tsv)
print("FASTA:", fasta)

print("\nExtraction statuses:")
for status, count in sorted(
    Counter(r["status"] for r in records).items()
):
    print(f"{status}: {count}")

print("\nRecords with overlapping CDS:",
      sum(bool(r["overlapping_CDS"]) for r in records))

print("\nRecords per species:")
for sp in species:
    subset = [r for r in records if r["species"] == sp]
    valid = sum(r["status"] == "OK" for r in subset)
    print(f"{sp}: {valid}/{len(subset)} valid windows")
