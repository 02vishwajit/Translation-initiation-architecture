#!/usr/bin/env python3

import csv
import random
import re
from pathlib import Path
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq

SRC = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_sequence_features.tsv"
)
OUT = SRC.parent / "single_copy_dinucleotide_null.tsv"

N_SHUFFLES = 1000
SEED = 20261018
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


def dinucleotide_counts(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))


def dinucleotide_shuffle(seq, rng):
    """Randomized Eulerian-trail shuffle preserving dinucleotide counts."""
    if len(seq) < 3:
        return seq

    adjacency = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        adjacency[a].append(b)

    for node in adjacency:
        rng.shuffle(adjacency[node])

    # Hierholzer's algorithm for an Eulerian trail.
    stack = [seq[0]]
    path = []
    local = {node: list(edges) for node, edges in adjacency.items()}

    while stack:
        node = stack[-1]
        if local.get(node):
            stack.append(local[node].pop())
        else:
            path.append(stack.pop())

    shuffled = "".join(reversed(path))

    if len(shuffled) != len(seq):
        raise RuntimeError("Shuffle length mismatch")
    if Counter(shuffled) != Counter(seq):
        raise RuntimeError("Shuffle did not preserve mononucleotide counts")
    if dinucleotide_counts(shuffled) != dinucleotide_counts(seq):
        raise RuntimeError("Shuffle did not preserve dinucleotide counts")

    return shuffled


def main():
    rng = random.Random(SEED)

    with SRC.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    output = []
    summary = defaultdict(lambda: {
        "n": 0, "observed": 0, "null_positive": 0
    })

    for idx, row in enumerate(rows, start=1):
        seq = row["upstream_sequence_5to3"].upper()
        if not seq or any(base not in "ACGT" for base in seq):
            raise ValueError(
                f"Unexpected or ambiguous sequence at row {idx}: {seq}"
            )

        observed = int(has_candidate(seq))
        null_positive = 0

        for _ in range(N_SHUFFLES):
            shuffled = dinucleotide_shuffle(seq, rng)
            null_positive += int(has_candidate(shuffled))

        null_probability = null_positive / N_SHUFFLES
        summary[row["species"]]["n"] += 1
        summary[row["species"]]["observed"] += observed
        summary[row["species"]]["null_positive"] += null_positive

        output.append({
            "species": row["species"],
            "orthogroup": row.get("orthogroup", ""),
            "sequence_id": row.get("sequence_id", ""),
            "observed_candidate": observed,
            "shuffle_positive_count": null_positive,
            "n_shuffles": N_SHUFFLES,
            "null_probability": f"{null_probability:.6f}",
            "observed_minus_null": f"{observed-null_probability:.6f}",
            "gc_fraction": row.get("gc_fraction", ""),
            "overlaps_cds": row.get("overlaps_cds", ""),
        })

        if idx % 100 == 0:
            print(f"Processed {idx}/{len(rows)} sequences")

    fields = list(output[0].keys())
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)

    print("\nSpecies-level observed vs dinucleotide-null expectation:")
    for species, values in sorted(summary.items()):
        n = values["n"]
        observed_pct = 100 * values["observed"] / n
        expected_pct = 100 * values["null_positive"] / (n * N_SHUFFLES)
        print(
            f"{species:18s} observed={observed_pct:6.2f}% "
            f"null_expected={expected_pct:6.2f}% "
            f"difference={observed_pct-expected_pct:+6.2f} percentage points"
        )

    print(f"\nOutput: {OUT}")
    print(f"Sequences: {len(rows)}")
    print(f"Shuffles per sequence: {N_SHUFFLES}")
    print("Dinucleotide-count preservation was checked for every shuffle.")


if __name__ == "__main__":
    main()
