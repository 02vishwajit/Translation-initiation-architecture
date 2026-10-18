#!/usr/bin/env python3

import csv
import re
from pathlib import Path

src = Path(
    "14_expanded_study/features/ccm_upstream/"
    "five_CCM_orthogroups_upstream50.tsv"
)
out = src.parent / "five_CCM_sequence_features.tsv"

# Short candidate motifs only; these are screening patterns,
# not proof of functional ribosome-binding sites.
motifs = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA", "GGAG")

def pct(n, d):
    return 100.0 * n / d if d else None

results = []

with src.open(newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        seq = r["upstream_sequence_5to3"].upper()
        n = len(seq)

        gc = pct(seq.count("G") + seq.count("C"), n)
        purine = pct(seq.count("A") + seq.count("G"), n)

        matches = []
        for motif in motifs:
            for m in re.finditer(f"(?={motif})", seq):
                start = m.start()
                end = start + len(motif)
                # Sequence ends directly before the annotated start codon.
                spacer = n - end
                matches.append((spacer, motif, start + 1))

        matches.sort(key=lambda x: (x[0], -len(x[1]), x[2]))
        nearest = matches[0] if matches else None

        r.update({
            "GC_percent": round(gc, 3) if gc is not None else "",
            "purine_percent": round(purine, 3) if purine is not None else "",
            "exact_candidate_motif_count": len(matches),
            "candidate_motifs": ";".join(
                f"{motif}@{pos}:spacer{spacer}"
                for spacer, motif, pos in matches
            ),
            "nearest_candidate_motif": nearest[1] if nearest else "",
            "nearest_motif_start_1based": nearest[2] if nearest else "",
            "nearest_motif_spacer_nt": nearest[0] if nearest else "",
            "overlap_flag": bool(r["other_CDS_overlapping_window"]),
        })
        results.append(r)

with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0].keys()),
                            delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

print("Output:", out)
print("Records:", len(results))
print("Sequences with at least one exact candidate motif:",
      sum(int(r["exact_candidate_motif_count"]) > 0 for r in results))
print("Sequences with overlapping CDS:",
      sum(r["overlap_flag"] for r in results))
print("\nMotif counts by orthogroup:")
for og in sorted({r["orthogroup"] for r in results}):
    subset = [r for r in results if r["orthogroup"] == og]
    count = sum(int(r["exact_candidate_motif_count"]) > 0 for r in subset)
    print(f"{og}: {count}/{len(subset)} sequences with a candidate motif")
