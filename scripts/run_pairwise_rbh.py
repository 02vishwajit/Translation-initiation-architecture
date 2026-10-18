from pathlib import Path
from itertools import combinations
import subprocess
import csv
import shutil

project = Path.home() / "alpha_translation_project"
base = project / "14_expanded_study"
proteomes = base / "genomes/standardized"
outdir = base / "orthology/pairwise_rbh"
dbdir = outdir / "blastdb"
logdir = base / "logs/blastp"

outdir.mkdir(parents=True, exist_ok=True)
dbdir.mkdir(parents=True, exist_ok=True)
logdir.mkdir(parents=True, exist_ok=True)

strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

for tool in ("makeblastdb", "blastp"):
    if shutil.which(tool) is None:
        raise SystemExit(f"Required executable not found: {tool}")

# Build protein databases.
for strain in strains:
    faa = proteomes / strain / f"{strain}.faa"
    if not faa.is_file() or faa.stat().st_size == 0:
        raise SystemExit(f"Missing/empty proteome: {faa}")

    prefix = dbdir / strain
    if not Path(str(prefix) + ".pin").exists():
        cmd = [
            "makeblastdb", "-in", str(faa), "-dbtype", "prot",
            "-parse_seqids", "-out", str(prefix)
        ]
        with (logdir / f"{strain}_makeblastdb.log").open("w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)

def best_hits(path):
    # BLAST outfmt:
    # qseqid sseqid pident length qlen slen evalue bitscore
    best = {}
    with path.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                continue
            q, s, pident, length, qlen, slen, evalue, bitscore = fields
            score = float(bitscore)
            row = {
                "query": q,
                "subject": s,
                "pident": float(pident),
                "alignment_length": int(length),
                "query_length": int(qlen),
                "subject_length": int(slen),
                "evalue": float(evalue),
                "bitscore": score,
            }
            if q not in best or score > best[q]["bitscore"]:
                best[q] = row
            elif score == best[q]["bitscore"] and s < best[q]["subject"]:
                # Deterministic tie-breaking; ties are not biological certainty.
                best[q] = row
    return best

def run_blast(query_strain, target_strain):
    path = outdir / f"{query_strain}__vs__{target_strain}.tsv"
    cmd = [
        "blastp",
        "-query", str(proteomes / query_strain / f"{query_strain}.faa"),
        "-db", str(dbdir / target_strain),
        "-evalue", "1e-5",
        "-seg", "yes",
        "-max_target_seqs", "10000",
        "-num_threads", "4",
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
        "-out", str(path),
    ]
    with (logdir / f"{query_strain}__vs__{target_strain}.log").open("w") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return path

pairs = list(combinations(strains, 2))
summary = []

for a, b in pairs:
    print(f"Comparing {a} <-> {b}", flush=True)
    ab_file = run_blast(a, b)
    ba_file = run_blast(b, a)

    ab = best_hits(ab_file)
    ba = best_hits(ba_file)

    rbh_rows = []
    for query, hit in ab.items():
        subject = hit["subject"]
        reverse = ba.get(subject)

        if reverse and reverse["subject"] == query:
            rbh_rows.append({
                "strain_a": a,
                "protein_a": query,
                "strain_b": b,
                "protein_b": subject,
                "identity_percent": hit["pident"],
                "alignment_length": hit["alignment_length"],
                "protein_a_length": hit["query_length"],
                "protein_b_length": hit["subject_length"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "coverage_a": hit["alignment_length"] / max(1, hit["query_length"]),
                "coverage_b": hit["alignment_length"] / max(1, hit["subject_length"]),
            })

    outfile = outdir / f"{a}__{b}__RBH.tsv"
    fields = [
        "strain_a", "protein_a", "strain_b", "protein_b",
        "identity_percent", "alignment_length",
        "protein_a_length", "protein_b_length", "evalue", "bitscore",
        "coverage_a", "coverage_b"
    ]
    with outfile.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rbh_rows)

    summary.append({
        "strain_a": a,
        "strain_b": b,
        "proteins_a_with_hits": len(ab),
        "proteins_b_with_hits": len(ba),
        "reciprocal_best_hit_pairs": len(rbh_rows),
        "rbh_file": str(outfile),
    })
    print(f"  RBH pairs: {len(rbh_rows)}", flush=True)

summary_file = outdir / "pairwise_rbh_summary.tsv"
with summary_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"\nSummary written: {summary_file}")
