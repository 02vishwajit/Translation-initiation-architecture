#!/usr/bin/env python3

import csv
import random
import re
from pathlib import Path
from collections import Counter, defaultdict

BASE = Path("14_expanded_study/features/single_copy_upstream")
SRC = BASE / "single_copy_sequence_features.tsv"

OUT_SEQ = BASE / "single_copy_dinucleotide_null_v3.tsv"
OUT_SPECIES = BASE / "single_copy_dinucleotide_species_v3.tsv"
OUT_REPLICATES = BASE / "single_copy_dinucleotide_replicates_v3.tsv"

N_SHUFFLES = 1000
SEED = 20261021

MOTIFS = ("AGGAGG", "AGGAG", "GGAGG", "GAGG", "AGGA", "GGAG")
MIN_SPACER = 4
MAX_SPACER = 12


def candidate_status(seq):
    """Accept any motif hit whose spacer is 4–12 nt."""
    for motif in MOTIFS:
        for match in re.finditer(f"(?={motif})", seq):
            end = match.start() + len(motif)
            spacer = len(seq) - end
            if MIN_SPACER <= spacer <= MAX_SPACER:
                return 1
    return 0


def dinucleotide_counts(seq):
    return Counter(seq[i:i + 2] for i in range(len(seq) - 1))


def dinucleotide_shuffle(seq, rng):
    """Randomized Eulerian traversal; verify exact dinucleotide preservation."""
    if len(seq) < 3:
        return seq

    adjacency = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        adjacency[a].append(b)

    for node in adjacency:
        rng.shuffle(adjacency[node])

    local = {node: list(edges) for node, edges in adjacency.items()}
    stack = [seq[0]]
    path = []

    while stack:
        node = stack[-1]
        if local.get(node):
            stack.append(local[node].pop())
        else:
            path.append(stack.pop())

    shuffled = "".join(reversed(path))

    if len(shuffled) != len(seq):
        raise RuntimeError("Sequence length changed during shuffling")
    if Counter(shuffled) != Counter(seq):
        raise RuntimeError("Mononucleotide counts changed")
    if dinucleotide_counts(shuffled) != dinucleotide_counts(seq):
        raise RuntimeError("Dinucleotide counts changed")

    return shuffled


def main():
    rng = random.Random(SEED)

    with SRC.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    if not rows:
        raise RuntimeError(f"No rows found in {SRC}")

    required = {"species", "orthogroup", "upstream_sequence_5to3"}
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Missing required columns: {sorted(missing)}")

    print(f"Loaded {len(rows)} sequences from {SRC}")

    species_rows = defaultdict(list)
    output_rows = []
    classification_changes = 0

    # Validate input and calculate corrected observed classifications.
    for row in rows:
        seq = row["upstream_sequence_5to3"].upper()

        if not seq or any(base not in "ACGT" for base in seq):
            raise ValueError(
                f"Invalid DNA sequence: {row.get('species')} "
                f"{row.get('orthogroup')}"
            )

        observed = candidate_status(seq)
        old_value = row.get("has_candidate_spacer_4_12", "")
        old_value = str(old_value).strip().lower() in {"true", "1", "yes"}

        if observed != int(old_value):
            classification_changes += 1

        new_row = dict(row)
        new_row["corrected_observed_candidate"] = observed
        new_row["corrected_motif_rule"] = "filter_spacer_before_classification"
        new_row["n_shuffles"] = N_SHUFFLES
        new_row["shuffle_positive_count"] = 0
        species_rows[row["species"]].append(new_row)
        output_rows.append(new_row)

    print(f"Observed classifications changed: {classification_changes}")

    # Species-level replicate counts: one count per complete shuffle replicate.
    replicate_rows = []
    null_sums = defaultdict(int)

    for replicate in range(1, N_SHUFFLES + 1):
        replicate_counts = defaultdict(int)

        for row in output_rows:
            seq = row["upstream_sequence_5to3"].upper()
            shuffled = dinucleotide_shuffle(seq, rng)
            positive = candidate_status(shuffled)

            row["shuffle_positive_count"] += positive
            replicate_counts[row["species"]] += positive

        for species, group in species_rows.items():
            count = replicate_counts[species]
            null_sums[species] += count
            replicate_rows.append({
                "replicate": replicate,
                "species": species,
                "n_sequences": len(group),
                "shuffle_positive": count,
                "shuffle_percent": 100 * count / len(group),
            })

        if replicate % 100 == 0:
            print(f"Completed shuffle replicate {replicate}/{N_SHUFFLES}")

    # Per-sequence output.
    for row in output_rows:
        row["null_probability"] = (
            row["shuffle_positive_count"] / N_SHUFFLES
        )
        row["observed_minus_null"] = (
            row["corrected_observed_candidate"] - row["null_probability"]
        )

    with OUT_SEQ.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    with OUT_REPLICATES.open("w", newline="") as f:
        fields = ["replicate", "species", "n_sequences",
                  "shuffle_positive", "shuffle_percent"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(replicate_rows)

    # Empirical upper-tail tests based on complete replicate-level species counts.
    species_results = []

    for species, group in sorted(species_rows.items()):
        observed = sum(r["corrected_observed_candidate"] for r in group)
        null_counts = [
            r["shuffle_positive"]
            for r in replicate_rows
            if r["species"] == species
        ]
        exceed = sum(count >= observed for count in null_counts)
        p_empirical = (exceed + 1) / (len(null_counts) + 1)
        expected = sum(null_counts) / len(null_counts)

        species_results.append({
            "species": species,
            "n_sequences": len(group),
            "observed_positive": observed,
            "observed_percent": 100 * observed / len(group),
            "null_expected_positive": expected,
            "null_expected_percent": 100 * expected / len(group),
            "observed_minus_expected": observed - expected,
            "empirical_p_upper_tail": p_empirical,
            "n_replicates": len(null_counts),
            "seed": SEED,
        })

    # Benjamini–Hochberg correction across the six species-level tests.
    order = sorted(
        range(len(species_results)),
        key=lambda i: species_results[i]["empirical_p_upper_tail"]
    )
    m = len(order)
    running = 1.0

    for rank_index in range(m - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        p = species_results[idx]["empirical_p_upper_tail"]
        running = min(running, p * m / rank)
        species_results[idx]["BH_FDR"] = min(1.0, running)

    with OUT_SPECIES.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(species_results[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(species_results)

    print("\nSpecies-level results (corrected motif rule)")
    for r in species_results:
        print(
            f'{r["species"]:18s} '
            f'observed={r["observed_percent"]:6.2f}% '
            f'null={r["null_expected_percent"]:6.2f}% '
            f'p={r["empirical_p_upper_tail"]:.4g} '
            f'FDR={r["BH_FDR"]:.4g}'
        )

    print("\nNew files written:")
    print(OUT_SEQ)
    print(OUT_SPECIES)
    print(OUT_REPLICATES)
    print("Existing v2 files were not modified.")


if __name__ == "__main__":
    main()
