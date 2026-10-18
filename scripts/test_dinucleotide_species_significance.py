#!/usr/bin/env python3

import csv
import random
from collections import defaultdict
from pathlib import Path

SRC = Path(
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_dinucleotide_null.tsv"
)
OUT = SRC.parent / "single_copy_dinucleotide_species_significance.tsv"

N_SIM = 200_000
SEED = 20261019
rng = random.Random(SEED)

groups = defaultdict(list)

with SRC.open(newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        groups[row["species"]].append(row)

results = []

for species, rows in sorted(groups.items()):
    observed = sum(int(r["observed_candidate"]) for r in rows)
    probabilities = [float(r["null_probability"]) for r in rows]

    simulated_counts = []
    for _ in range(N_SIM):
        total = sum(rng.random() < p for p in probabilities)
        simulated_counts.append(total)

    exceedances = sum(x >= observed for x in simulated_counts)
    p_empirical = (exceedances + 1) / (N_SIM + 1)

    null_mean = sum(probabilities)
    null_variance = sum(p * (1 - p) for p in probabilities)

    results.append({
        "species": species,
        "n_sequences": len(rows),
        "observed_positive": observed,
        "observed_percent": f"{100 * observed / len(rows):.3f}",
        "null_expected_positive": f"{null_mean:.3f}",
        "null_expected_percent": f"{100 * null_mean / len(rows):.3f}",
        "observed_minus_expected": f"{observed - null_mean:.3f}",
        "null_sd": f"{null_variance ** 0.5:.3f}",
        "empirical_p_upper_tail": f"{p_empirical:.8g}",
        "simulations": N_SIM,
    })

fields = list(results[0].keys())
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

print("Species-level Monte Carlo test")
print("Note: p-value resolution is limited by the simulation count.\n")
for r in results:
    print(
        f'{r["species"]:18s} '
        f'observed={r["observed_percent"]:>7s}% '
        f'null={r["null_expected_percent"]:>7s}% '
        f'p={r["empirical_p_upper_tail"]}'
    )

print(f"\nOutput: {OUT}")
