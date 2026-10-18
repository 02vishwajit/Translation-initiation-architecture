from pathlib import Path
from itertools import combinations
import csv

base = Path("14_expanded_study/orthology/pairwise_rbh")

strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

def normalize_id(identifier):
    """Normalize BLAST identifiers like ref|WP_123.1| to WP_123.1."""
    identifier = identifier.strip()
    if "|" in identifier:
        parts = identifier.split("|")
        # BLAST database identifiers often use database|accession|description
        if len(parts) >= 3 and parts[1]:
            return parts[1]
    return identifier.split()[0]

def best_hits(path):
    best = {}

    with path.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                continue

            q, s, pident, length, qlen, slen, evalue, bitscore = fields
            q = normalize_id(q)
            s = normalize_id(s)
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

            old = best.get(q)
            if (
                old is None
                or score > old["bitscore"]
                or (score == old["bitscore"] and s < old["subject"])
            ):
                best[q] = row

    return best

summary = []

for a, b in combinations(strains, 2):
    ab_file = base / f"{a}__vs__{b}.tsv"
    ba_file = base / f"{b}__vs__{a}.tsv"

    if not ab_file.exists() or not ba_file.exists():
        print(f"WARNING: missing BLAST result for {a} / {b}")
        continue

    ab = best_hits(ab_file)
    ba = best_hits(ba_file)
    rows = []

    for q, hit in ab.items():
        s = hit["subject"]
        reverse = ba.get(s)

        if reverse and reverse["subject"] == q:
            rows.append({
                "strain_a": a,
                "protein_a": q,
                "strain_b": b,
                "protein_b": s,
                "identity_percent": hit["pident"],
                "alignment_length": hit["alignment_length"],
                "protein_a_length": hit["query_length"],
                "protein_b_length": hit["subject_length"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "coverage_a": hit["alignment_length"] / max(1, hit["query_length"]),
                "coverage_b": hit["alignment_length"] / max(1, hit["subject_length"]),
            })

    out = base / f"{a}__{b}__RBH.tsv"
    fields = [
        "strain_a", "protein_a", "strain_b", "protein_b",
        "identity_percent", "alignment_length",
        "protein_a_length", "protein_b_length", "evalue", "bitscore",
        "coverage_a", "coverage_b"
    ]

    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary.append({
        "strain_a": a,
        "strain_b": b,
        "proteins_a_with_hits": len(ab),
        "proteins_b_with_hits": len(ba),
        "reciprocal_best_hit_pairs": len(rows),
        "rbh_file": str(out),
    })

    print(f"{a} vs {b}: {len(rows)} RBH pairs")

summary_file = base / "pairwise_rbh_summary.tsv"
fields = [
    "strain_a", "strain_b", "proteins_a_with_hits",
    "proteins_b_with_hits", "reciprocal_best_hit_pairs", "rbh_file"
]

with summary_file.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"\nUpdated summary: {summary_file}")
