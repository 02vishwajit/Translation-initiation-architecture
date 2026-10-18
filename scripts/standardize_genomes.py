from pathlib import Path
import gzip
import shutil
import csv
import json

PROJECT = Path.home() / "alpha_translation_project"
RAW = PROJECT / "00_raw_data"
NEW = PROJECT / "14_expanded_study/genomes"
OUT = NEW / "standardized"
OUT.mkdir(parents=True, exist_ok=True)

# ID: source directory, assembly accession
sources = {
    "ZM4": (RAW / "zm4", "GCF_000007105.1"),
    "NOVO": (RAW / "novo", "GCF_000013325.1"),
    "RSP": (RAW / "rsphaeroides", "GCF_000012905.2"),
    "ZM_NCIMB11163": (
        NEW / "ZM_NCIMB11163/extracted/ncbi_dataset/data/GCF_000024245.1",
        "GCF_000024245.1"
    ),
    "TSA_JW_SL_YS485": (
        NEW / "TSA_JW_SL_YS485/extracted/ncbi_dataset/data/GCF_000307585.2",
        "GCF_000307585.2"
    ),
    "CT_ATCC27405": (
        NEW / "CT_ATCC27405/extracted/ncbi_dataset/data/GCF_000015865.1",
        "GCF_000015865.1"
    ),
}

def open_text(path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()

def find_file(folder, accession, kind):
    patterns = {
        "fna": [f"{accession}*genomic.fna", "genomic.fna"],
        "faa": ["protein.faa"],
        "gff": ["genomic.gff", "genomic.gff3"],
        "gbff": ["genomic.gbff"],
    }

    for pattern in patterns[kind]:
        matches = list(folder.glob(pattern))
        if matches:
            return matches[0]

    # Original files may have the assembly accession in their filename.
    for pattern in patterns[kind]:
        matches = list(folder.glob(pattern + ".gz"))
        if matches:
            return matches[0]

    return None

def fasta_stats(path):
    records = 0
    lengths = []
    current = 0

    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                if records:
                    lengths.append(current)
                records += 1
                current = 0
            else:
                current += len(line.strip())

    if records:
        lengths.append(current)

    return records, sum(lengths), min(lengths) if lengths else 0

summary = []

for strain, (folder, accession) in sources.items():
    if not folder.exists():
        raise FileNotFoundError(f"Missing source directory: {folder}")

    target = OUT / strain
    target.mkdir(parents=True, exist_ok=True)

    row = {
        "strain_id": strain,
        "assembly_accession": accession,
        "source_directory": str(folder),
    }

    for kind in ["fna", "faa", "gff", "gbff"]:
        src = find_file(folder, accession, kind)
        if src is None:
            print(f"WARNING: {strain}: missing {kind} source")
            row[kind + "_source"] = ""
            continue

        dest = target / f"{strain}.{kind}"
        with open_text(src) as source, dest.open("w") as output:
            shutil.copyfileobj(source, output)

        row[kind + "_source"] = str(src)
        row[kind + "_bytes"] = dest.stat().st_size

        if kind in ["fna", "faa"]:
            n, total, min_len = fasta_stats(dest)
            row[kind + "_records"] = n
            row[kind + "_total_length"] = total
            row[kind + "_minimum_record_length"] = min_len
            if n == 0 or total == 0:
                raise ValueError(f"Empty FASTA detected: {dest}")

    summary.append(row)
    print(f"Prepared {strain}")

fields = sorted({key for row in summary for key in row})
with (OUT / "genome_file_summary.tsv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"\nSummary written to {OUT / 'genome_file_summary.tsv'}")
