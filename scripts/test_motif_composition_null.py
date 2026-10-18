#!/usr/bin/env python3

import csv
import random
import re
from pathlib import Path
from collections import defaultdict

src = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_sequence_features.tsv"
)
outdir = src.parent
N_SHUFFLES = 1000
SEED = 20261017

# Keep the motif definition and spacer range identical to the observed scan.
MOTIFS = ("AGGAGG", "AGGAG", "GGAGG", "GAGG", "AGGA", "GGAG")
MIN_SPACER = 4
MAX_SPACER = 12

def has_candidate(seq):
    n = len(seq)
    for motif in MOTIFS:
        for match in re.finditer(f"(?={motif})", seq):
            spacer = n - (match.start() + len(motif))
            if MIN_SPACER <= spacer <= MAX_SPACER:
                return True
    return False

def shuffled_sequence(seq, rng):
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)

rng = random.Random(SEED)

with src.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

results = []
species_stats = defaultdict(lambda: {
    "n": 0, "observed": 0, "null_positive_total": 0,
    "null_probability_sum": 0.0
})

for i, row in enumerate(rows, start=1):
    seq = row["upstream_sequence_5to3"].upper()
    observed = has_candidate(seq)

    null_positive = 0
    for _ in range(N_SHUFFLES):
        shuffled = shuffled_sequence(seq, rng)
        null_positive += has_candidate(shuffled)

    null_probability = null_positive / N_SHUFFLES

    # Plus-one correction for a simple empirical comparison.
    # This is not a multiple-testing-adjusted P value.
    empirical_p = (
        1 + sum(
            has_candidate(shuffled_sequence(seq, rng))
            for _ in range(0)
        )
    ) / (N_SHUFFLES + 1)

    results.append({
        "orthogroup": row["orthogroup"],
        "species": row["species"],
        "protein_id": row["protein_id"],
        "observed_candidate": int(observed),
        "shuffle_positive_count": null_positive,
        "shuffle_probability": round(null_probability, 4),
        "observed_minus_shuffle_probability":
            round(int(observed) - null_probability, 4),
        "overlap_flag": row["overlap_flag"],
        "GC_percent": row["GC_percent"],
        "purine_percent": row["purine_percent"],
    })

    s = species_stats[row["species"]]
    s["n"] += 1
    s["observed"] += int(observed)
    s["null_positive_total"] += null_positive
    s["null_probability_sum"] += null_probability

outfile = outdir / "single_copy_motif_composition_null.tsv"
with outfile.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0].keys()), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("Output:", outfile)
print("Sequences:", len(results))
print("Shuffles per sequence:", N_SHUFFLES)
print("\nSpecies-level observed vs composition-null expectation:")

for sp in sorted(species_stats):
    s = species_stats[sp]
    n = s["n"]
    observed_pct = 100 * s["observed"] / n
    expected_pct = 100 * s["null_probability_sum"] / n
    difference = observed_pct - expected_pct

    print(
        f"{sp:18s} observed={observed_pct:6.2f}% "
        f"null_expected={expected_pct:6.2f}% "
        f"difference={difference:+6.2f} percentage points"
    )

print("\nNote: mononucleotide-preserving null only;")
print("not a dinucleotide-preserving model or a functional RBS test.")
