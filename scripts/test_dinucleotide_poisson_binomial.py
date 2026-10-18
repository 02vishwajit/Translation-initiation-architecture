#!/usr/bin/env python3

import csv
import math
from pathlib import Path
from collections import defaultdict

BASE = Path("14_expanded_study/features/single_copy_upstream")
SRC = BASE / "single_copy_dinucleotide_null_v2.tsv"
OUT = BASE / "single_copy_dinucleotide_poisson_binomial.tsv"


def poisson_binomial_upper_tail(probabilities, observed):
    """
    Exact dynamic-programming Poisson-binomial tail,
    conditional on the supplied per-sequence probabilities.
    """
    dp = [1.0]

    for p in probabilities:
        new = [0.0] * (len(dp) + 1)

        for k, probability in enumerate(dp):
            new[k] += probability * (1.0 - p)
            new[k + 1] += probability * p

        dp = new

    return min(1.0, max(0.0, sum(dp[observed:])))


def bh_adjust(pvalues):
    """Benjamini-Hochberg adjustment."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adjusted = [1.0] * m
    running = 1.0

    for rank_index in range(m - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(
            running,
            pvalues[original_index] * m / rank
        )
        adjusted[original_index] = min(1.0, running)

    return adjusted


groups = defaultdict(list)

with SRC.open(newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        groups[row["species"]].append(row)

results = []

for species, rows in sorted(groups.items()):
    observed = sum(int(r["observed_candidate"]) for r in rows)

    probabilities = [
        int(r["shuffle_positive_count"]) / int(r["n_shuffles"])
        for r in rows
    ]

    expected = sum(probabilities)
    variance = sum(p * (1 - p) for p in probabilities)
    sd = math.sqrt(variance)

    pvalue = poisson_binomial_upper_tail(
        probabilities, observed
    )

    results.append({
        "species": species,
        "n_sequences": len(rows),
        "observed_positive": observed,
        "observed_percent": 100 * observed / len(rows),
        "null_expected_positive": expected,
        "null_expected_percent": 100 * expected / len(rows),
        "observed_minus_expected": observed - expected,
        "null_sd": sd,
        "p_value_conditional": pvalue,
    })

adjusted = bh_adjust([r["p_value_conditional"] for r in results])

for row, qvalue in zip(results, adjusted):
    row["BH_FDR"] = qvalue

fields = list(results[0].keys())

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for row in results:
        writer.writerow({
            k: (
                f"{v:.10g}" if isinstance(v, float) else v
            )
            for k, v in row.items()
        })

print("Species-level Poisson-binomial tests")
print("P-values are conditional on estimated shuffle probabilities.\n")

for row in results:
    print(
        f'{row["species"]:18s} '
        f'observed={row["observed_percent"]:6.2f}% '
        f'null={row["null_expected_percent"]:6.2f}% '
        f'P={row["p_value_conditional"]:.4g} '
        f'FDR={row["BH_FDR"]:.4g}'
    )

print(f"\nOutput: {OUT}")
