#!/usr/bin/env python3

import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

BASE = Path("14_expanded_study")
GENOMES = BASE / "genomes/standardized"
RESULTS = (
    BASE / "orthology/orthofinder_input/OrthoFinder/"
    "Results_six_genome_v1"
)
OG_FILE = RESULTS / "Orthogroups/Orthogroups.tsv"
OUT = BASE / "features/ccm_upstream"
OUT.mkdir(parents=True, exist_ok=True)

WINDOW = 50
TARGET_OGS = {
    "OG0000678",  # phosphoglycerate kinase
    "OG0000449",  # triose-phosphate isomerase
    "OG0000226",  # glyceraldehyde-3-phosphate dehydrogenase
    "OG0000462",  # glucose-6-phosphate isomerase
    "OG0000450",  # enolase
}

def parse_proteins(cell):
    return [x.strip() for x in (cell or "").split(",") if x.strip()]

# Read target orthogroups and protein membership.
with OG_FILE.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    species = [x for x in fields if x != "Orthogroup"]
    target_rows = [
        r for r in reader if r["Orthogroup"] in TARGET_OGS
    ]

wanted = {}
for row in target_rows:
    og = row["Orthogroup"]
    for sp in species:
        for pid in parse_proteins(row.get(sp, "")):
            wanted.setdefault((sp, pid), []).append(og)

# Match protein accessions to CDS features in each GenBank file.
records = []
seen = set()

for sp in species:
    gbff = GENOMES / sp / f"{sp}.gbff"
    if not gbff.exists():
        print(f"WARNING: missing GenBank file: {gbff}")
        continue

    for rec in SeqIO.parse(str(gbff), "genbank"):
        cds_features = [
            f for f in rec.features if f.type == "CDS"
        ]

        for feature in cds_features:
            protein_ids = feature.qualifiers.get("protein_id", [])
            if not protein_ids:
                continue

            pid = protein_ids[0].strip()
            key = (sp, pid)
            if key not in wanted:
                continue

            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            contig_len = len(rec.seq)

            if strand == 1:
                upstream_start = start - WINDOW
                upstream_end = start
                if upstream_start < 0:
                    upstream = ""
                    status = "CONTIG_BOUNDARY_INCOMPLETE"
                else:
                    upstream = str(rec.seq[upstream_start:upstream_end]).upper()
                    status = "OK"
            elif strand == -1:
                upstream_start = end
                upstream_end = end + WINDOW
                if upstream_end > contig_len:
                    upstream = ""
                    status = "CONTIG_BOUNDARY_INCOMPLETE"
                else:
                    upstream = str(
                        rec.seq[upstream_start:upstream_end].reverse_complement()
                    ).upper()
                    status = "OK"
            else:
                upstream = ""
                status = "UNKNOWN_STRAND"

            cds_seq = str(feature.extract(rec.seq)).upper()
            start_codon = cds_seq[:3] if len(cds_seq) >= 3 else ""
            locus = feature.qualifiers.get("locus_tag", [""])[0]
            gene = feature.qualifiers.get("gene", [""])[0]
            product = feature.qualifiers.get("product", [""])[0]

            # Identify other CDS features overlapping the upstream genomic interval.
            overlaps = []
            if status == "OK":
                for other in cds_features:
                    if other is feature:
                        continue
                    os = int(other.location.start)
                    oe = int(other.location.end)
                    if os < upstream_end and oe > upstream_start:
                        overlaps.append(
                            other.qualifiers.get("locus_tag", ["unknown"])[0]
                        )

            row_out = {
                "orthogroup": ";".join(wanted[key]),
                "species": sp,
                "protein_id": pid,
                "locus_tag": locus,
                "gene": gene,
                "product": product,
                "contig": rec.id,
                "strand": strand,
                "cds_start_0based": start,
                "cds_end_0based": end,
                "start_codon": start_codon,
                "upstream_window_requested_nt": WINDOW,
                "upstream_window_length_nt": len(upstream),
                "upstream_sequence_5to3": upstream,
                "upstream_status": status,
                "other_CDS_overlapping_window": ";".join(sorted(set(overlaps))),
            }
            records.append(row_out)
            seen.add(key)

# Add explicit records for proteins that were not mapped to a CDS.
for (sp, pid), ogs in wanted.items():
    if (sp, pid) not in seen:
        records.append({
            "orthogroup": ";".join(ogs),
            "species": sp,
            "protein_id": pid,
            "locus_tag": "",
            "gene": "",
            "product": "",
            "contig": "",
            "strand": "",
            "cds_start_0based": "",
            "cds_end_0based": "",
            "start_codon": "",
            "upstream_window_requested_nt": WINDOW,
            "upstream_window_length_nt": 0,
            "upstream_sequence_5to3": "",
            "upstream_status": "PROTEIN_ID_NOT_FOUND_IN_GENBANK_CDS",
            "other_CDS_overlapping_window": "",
        })

outfile = OUT / "five_CCM_orthogroups_upstream50.tsv"
with outfile.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=list(records[0].keys()),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(records)

fasta = OUT / "five_CCM_orthogroups_upstream50.fasta"
with fasta.open("w") as handle:
    for r in records:
        if r["upstream_status"] == "OK" and r["upstream_sequence_5to3"]:
            header = (
                f'>{r["orthogroup"].replace(";", "_")}|{r["species"]}|'
                f'{r["protein_id"]}|{r["locus_tag"]}'
            )
            handle.write(header + "\n")
            seq = r["upstream_sequence_5to3"]
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i+60] + "\n")

print("Target orthogroups:", len(TARGET_OGS))
print("Expected protein records:", len(wanted))
print("Extracted/missing records:", len(records))
print("TSV:", outfile)
print("FASTA:", fasta)

print("\nStatus counts:")
from collections import Counter
for status, n in sorted(Counter(r["upstream_status"] for r in records).items()):
    print(f"{status}: {n}")

print("\nRecords per orthogroup:")
for og in sorted(TARGET_OGS):
    subset = [r for r in records if og in r["orthogroup"].split(";")]
    print(
        og,
        "records =", len(subset),
        "valid windows =", sum(r["upstream_status"] == "OK" for r in subset),
    )
