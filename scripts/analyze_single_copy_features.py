#!/usr/bin/env python3

import csv
import re
from pathlib import Path
from collections import Counter

src = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_upstream50.tsv"
)
outdir = src.parent
outfile = outdir / "single_copy_sequence_features.tsv"

MOTIFS = ("AGGAGG", "AGGAG", "GGAGG", "GAGG", "AGGA", "GGAG")
MIN_SPACER = 4
MAX_SPACER = 12
TARGET_SPACER = 7

def find_hits(seq):
    hits = []
    for motif in MOTIFS:
        for match in re.finditer(f"(?={motif})", seq):
            start = match.start()
            end = start + len(motif)
            hits.append({
                "motif": motif,
                "start": start,
                "end": end,
                "length": len(motif),
                "spacer": len(seq) - end,
            })
    return hits

def collapse_overlapping_hits(hits):
    """Group overlapping motif intervals and retain one representative
    per interval cluster. Prefer the longest motif, then a spacer near 7."""
    if not hits:
        return []

    hits = sorted(hits, key=lambda h: (h["start"], h["end"]))
    clusters = []
    cluster = [hits[0]]
    cluster_end = hits[0]["end"]

    for hit in hits[1:]:
        if hit["start"] < cluster_end:
            cluster.append(hit)
            cluster_end = max(cluster_end, hit["end"])
        else:
            clusters.append(cluster)
            cluster = [hit]
            cluster_end = hit["end"]
    clusters.append(cluster)

    representatives = []
    for cluster in clusters:
        representative = sorted(
            cluster,
            key=lambda h: (
                -h["length"],
                abs(h["spacer"] - TARGET_SPACER),
                h["start"],
            )
        )[0]
        representatives.append(representative)

    return representatives

rows = []
with src.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        seq = r["upstream_sequence_5to3"].upper()
        n = len(seq)

        gc_count = seq.count("G") + seq.count("C")
        purine_count = seq.count("A") + seq.count("G")

        hits = find_hits(seq)
        representatives = collapse_overlapping_hits(hits)

        plausible = [
            h for h in representatives
            if MIN_SPACER <= h["spacer"] <= MAX_SPACER
        ]
        plausible.sort(
            key=lambda h: (
                abs(h["spacer"] - TARGET_SPACER),
                -h["length"],
                h["start"],
            )
        )
        best = plausible[0] if plausible else None

        r.update({
            "GC_percent": round(100 * gc_count / n, 3) if n else "",
            "purine_percent": round(100 * purine_count / n, 3) if n else "",
            "raw_motif_match_count": len(hits),
            "overlap_collapsed_candidate_count": len(representatives),
            "candidate_motifs_spacer_4_12": ";".join(
                f'{h["motif"]}@{h["start"]+1}:spacer{h["spacer"]}'
                for h in plausible
            ),
            "has_candidate_spacer_4_12": bool(plausible),
            "best_candidate_motif": best["motif"] if best else "",
            "best_candidate_spacer_nt": best["spacer"] if best else "",
            "best_candidate_start_1based": best["start"] + 1 if best else "",
            "overlap_flag": bool(r["overlapping_CDS"]),
            "start_codon_class": (
                "ATG" if r["start_codon"] == "ATG"
                else "Alternative_start"
                if r["start_codon"] in {"GTG", "TTG"}
                else "Other_or_unexpected"
            ),
        })
        rows.append(r)

with outfile.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0].keys()), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("Output:", outfile)
print("Total rows:", len(rows))
print("Unique orthogroups:", len({r["orthogroup"] for r in rows}))
print("Sequence length distribution:",
      dict(Counter(r["upstream_length_nt"] for r in rows)))
print("Candidate motif at spacer 4–12:",
      sum(r["has_candidate_spacer_4_12"] for r in rows))
print("Overlapping CDS windows:",
      sum(r["overlap_flag"] for r in rows))

print("\nBy species:")
for sp in sorted({r["species"] for r in rows}):
    group = [r for r in rows if r["species"] == sp]
    n = sum(r["has_candidate_spacer_4_12"] for r in group)
    print(
        f"{sp}: motif={n}/{len(group)}, "
        f"mean_GC={sum(float(r['GC_percent']) for r in group)/len(group):.2f}%"
    )

print("\nBy orthogroup, sequences with a candidate motif:")
for og in sorted({r["orthogroup"] for r in rows}):
    group = [r for r in rows if r["orthogroup"] == og]
    n = sum(r["has_candidate_spacer_4_12"] for r in group)
    if n == 0 or n == len(group):
        continue
    print(f"{og}: {n}/{len(group)}")
