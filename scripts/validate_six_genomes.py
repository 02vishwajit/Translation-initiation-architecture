from pathlib import Path
import csv
import re

base = Path("14_expanded_study/genomes/standardized")
strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

def read_fasta(path):
    records = []
    header = None
    chunks = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first header in {path}")
                chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks)))
    return records

rows = []

for strain in strains:
    folder = base / strain
    genome = folder / f"{strain}.fna"
    protein = folder / f"{strain}.faa"
    gff = folder / f"{strain}.gff"
    gbff = folder / f"{strain}.gbff"

    row = {"strain_id": strain}
    missing = [p.name for p in [genome, protein, gff, gbff] if not p.is_file()]
    row["missing_files"] = ";".join(missing)

    if genome.is_file():
        dna = read_fasta(genome)
        lengths = [len(seq) for _, seq in dna]
        sequence = "".join(seq.upper() for _, seq in dna)
        row.update({
            "genome_records": len(dna),
            "genome_total_bp": sum(lengths),
            "genome_min_record_bp": min(lengths, default=0),
            "genome_max_record_bp": max(lengths, default=0),
            "genome_GC_percent": round(
                100 * (sequence.count("G") + sequence.count("C")) /
                max(1, sum(c in "ACGT" for c in sequence)), 3
            ),
            "genome_non_ACGT_characters": sum(c not in "ACGTN" for c in sequence),
            "duplicate_genome_headers": len(dna) - len(set(h.split()[0] for h, _ in dna)),
        })

    if protein.is_file():
        proteins = read_fasta(protein)
        lengths = [len(seq) for _, seq in proteins]
        allowed = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
        row.update({
            "protein_records": len(proteins),
            "protein_total_aa": sum(lengths),
            "protein_min_aa": min(lengths, default=0),
            "protein_max_aa": max(lengths, default=0),
            "duplicate_protein_headers": len(proteins) -
                len(set(h.split()[0] for h, _ in proteins)),
            "protein_invalid_characters": sum(
                sum(c.upper() not in allowed for c in seq)
                for _, seq in proteins
            ),
            "protein_internal_stops": sum(
                seq[:-1].count("*") for _, seq in proteins
            ),
        })

    if gff.is_file():
        cds = 0
        genes = 0
        with gff.open() as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    genes += fields[2] == "gene"
                    cds += fields[2] == "CDS"
        row["gff_gene_features"] = genes
        row["gff_CDS_features"] = cds

    if gbff.is_file():
        text = gbff.read_text(errors="replace")
        row["gbff_CDS_features"] = len(re.findall(r"^\s{5}CDS\s+", text, re.M))
        row["gbff_record_count"] = len(re.findall(r"^LOCUS\s+", text, re.M))

    rows.append(row)

out = base / "six_genome_validation.tsv"
fields = sorted({k for r in rows for k in r})
with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Validation report:", out)
for r in rows:
    print(
        r["strain_id"],
        "| genome records:", r.get("genome_records", "NA"),
        "| total bp:", r.get("genome_total_bp", "NA"),
        "| proteins:", r.get("protein_records", "NA"),
        "| GFF CDS:", r.get("gff_CDS_features", "NA"),
        "| GBFF CDS:", r.get("gbff_CDS_features", "NA"),
        "| missing:", r.get("missing_files", "")
    )
