#!/usr/bin/env python3

import csv
import random
import re
from pathlib import Path
from collections import defaultdict, Counter

SRC = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_sequence_features.tsv"
)

OUT = SRC.parent / "single_copy_dinucleotide_null_v2.tsv"

N_SHUFFLES = 1000
SEED = 20261020

MOTIFS = (
    "AGGAGG",
    "AGGAG",
    "GGAGG",
    "GAGG",
    "AGGA",
    "GGAG",
)

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
    """
    EXACTLY reproduce the logic used in
    analyze_single_copy_features.py.
    """
    if not hits:
        return []

    hits = sorted(
        hits,
        key=lambda h: (h["start"], h["end"])
    )

    clusters = []
    cluster = [hits[0]]
    cluster_end = hits[0]["end"]

    for hit in hits[1:]:
        if hit["start"] < cluster_end:
            cluster.append(hit)
            cluster_end = max(
                cluster_end,
                hit["end"]
            )
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

        representatives.append(
            representative
        )

    return representatives


def candidate_status(seq):
    hits = find_hits(seq)

    representatives = collapse_overlapping_hits(
        hits
    )

    plausible = [
        h for h in representatives
        if MIN_SPACER <= h["spacer"] <= MAX_SPACER
    ]

    return int(bool(plausible))


def dinucleotide_counts(seq):
    return Counter(
        seq[i:i+2]
        for i in range(len(seq)-1)
    )


def dinucleotide_shuffle(seq, rng):
    """
    Randomized Eulerian traversal preserving
    exact dinucleotide counts.
    """

    if len(seq) < 3:
        return seq

    adjacency = defaultdict(list)

    for a, b in zip(seq[:-1], seq[1:]):
        adjacency[a].append(b)

    for node in adjacency:
        rng.shuffle(adjacency[node])

    local = {
        node: list(edges)
        for node, edges in adjacency.items()
    }

    stack = [seq[0]]
    path = []

    while stack:
        node = stack[-1]

        if local.get(node):
            stack.append(
                local[node].pop()
            )
        else:
            path.append(
                stack.pop()
            )

    shuffled = "".join(
        reversed(path)
    )

    # Mandatory integrity checks.
    if len(shuffled) != len(seq):
        raise RuntimeError(
            "Sequence-length mismatch"
        )

    if Counter(shuffled) != Counter(seq):
        raise RuntimeError(
            "Mononucleotide counts changed"
        )

    if (
        dinucleotide_counts(shuffled)
        != dinucleotide_counts(seq)
    ):
        raise RuntimeError(
            "Dinucleotide counts changed"
        )

    return shuffled


def as_bool(value):
    return str(value).strip().lower() == "true"


def main():

    rng = random.Random(SEED)

    with SRC.open(newline="") as f:
        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    print(
        f"Loaded {len(rows)} sequences."
    )

    # ------------------------------------------------
    # First validate observed classifications.
    # ------------------------------------------------

    mismatches = []

    for row in rows:

        seq = (
            row["upstream_sequence_5to3"]
            .upper()
        )

        recalculated = candidate_status(seq)

        stored = int(
            as_bool(
                row[
                    "has_candidate_spacer_4_12"
                ]
            )
        )

        if recalculated != stored:
            mismatches.append({
                "species": row["species"],
                "orthogroup": row["orthogroup"],
                "protein_id": row["protein_id"],
                "stored": stored,
                "recalculated": recalculated,
            })

    print(
        "Observed classification mismatches:",
        len(mismatches)
    )

    if mismatches:

        for x in mismatches[:20]:
            print(x)

        raise RuntimeError(
            "Observed motif detector does not "
            "match original feature analysis."
        )

    print(
        "Observed motif classifications "
        "validated successfully."
    )

    # ------------------------------------------------
    # Dinucleotide null.
    # ------------------------------------------------

    output = []

    species_summary = defaultdict(
        lambda: {
            "n": 0,
            "observed": 0,
            "shuffle_positive": 0,
        }
    )

    for idx, row in enumerate(rows, start=1):

        seq = (
            row["upstream_sequence_5to3"]
            .upper()
        )

        if (
            not seq
            or any(
                base not in "ACGT"
                for base in seq
            )
        ):
            raise ValueError(
                f"Invalid sequence "
                f"{row['species']} "
                f"{row['orthogroup']}"
            )

        observed = candidate_status(seq)

        positive = 0

        for _ in range(N_SHUFFLES):

            shuffled = dinucleotide_shuffle(
                seq,
                rng
            )

            positive += candidate_status(
                shuffled
            )

        null_probability = (
            positive / N_SHUFFLES
        )

        sp = row["species"]

        species_summary[sp]["n"] += 1

        species_summary[sp][
            "observed"
        ] += observed

        species_summary[sp][
            "shuffle_positive"
        ] += positive

        output.append({
            "species":
                sp,

            "orthogroup":
                row["orthogroup"],

            "protein_id":
                row["protein_id"],

            "locus_tag":
                row["locus_tag"],

            "gene":
                row["gene"],

            "observed_candidate":
                observed,

            "shuffle_positive_count":
                positive,

            "n_shuffles":
                N_SHUFFLES,

            "null_probability":
                f"{null_probability:.6f}",

            "observed_minus_null":
                f"{observed-null_probability:.6f}",

            "GC_percent":
                row["GC_percent"],

            "purine_percent":
                row["purine_percent"],

            "overlap_flag":
                row["overlap_flag"],

            "start_codon":
                row["start_codon"],
        })

        if idx % 100 == 0:
            print(
                f"Processed "
                f"{idx}/{len(rows)}"
            )

    fields = list(
        output[0].keys()
    )

    with OUT.open(
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(output)

    # ------------------------------------------------
    # Species summary.
    # ------------------------------------------------

    print(
        "\nSpecies-level observed vs "
        "dinucleotide-null expectation:"
    )

    for sp in sorted(species_summary):

        x = species_summary[sp]

        n = x["n"]

        observed_pct = (
            100 *
            x["observed"] /
            n
        )

        null_pct = (
            100 *
            x["shuffle_positive"] /
            (n * N_SHUFFLES)
        )

        print(
            f"{sp:18s} "
            f"observed={observed_pct:6.2f}% "
            f"null_expected={null_pct:6.2f}% "
            f"difference="
            f"{observed_pct-null_pct:+6.2f} "
            f"percentage points"
        )

    print(
        f"\nOutput: {OUT}"
    )

    print(
        f"Sequences: {len(rows)}"
    )

    print(
        f"Shuffles per sequence: "
        f"{N_SHUFFLES}"
    )

    print(
        "Observed and shuffled sequences "
        "used identical motif logic."
    )

    print(
        "Exact dinucleotide preservation "
        "checked for every shuffle."
    )


if __name__ == "__main__":
    main()
