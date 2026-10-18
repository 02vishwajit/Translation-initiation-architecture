import csv
from pathlib import Path
from collections import defaultdict
from itertools import combinations

from scipy.stats import binomtest, chi2

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_sequence_features_corrected.tsv"
OUT_Q = BASE / "matched_cochran_q.tsv"
OUT_PAIR = BASE / "matched_pairwise_mcnemar.tsv"

with INPUT.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

species = sorted({r["species"] for r in rows})
groups = defaultdict(dict)

for r in rows:
    og = r["orthogroup"]
    sp = r["species"]
    groups[og][sp] = {
        "positive": r["has_candidate_spacer_4_12"].strip().lower() == "true",
        "overlap": r["overlap_flag"].strip().lower() == "true",
    }

if len(species) != 6 or len(groups) != 348:
    raise SystemExit("Expected six species and 348 orthogroups.")
if any(set(v) != set(species) for v in groups.values()):
    raise SystemExit("Incomplete orthogroup/species matching.")

# Cochran's Q statistic for matched binary outcomes.
k = len(species)
row_totals = []
col_totals = {sp: 0 for sp in species}

for og in sorted(groups):
    values = [int(groups[og][sp]["positive"]) for sp in species]
    row_totals.append(sum(values))
    for sp, value in zip(species, values):
        col_totals[sp] += value

N = len(row_totals)
T = sum(row_totals)
sum_col_sq = sum(v * v for v in col_totals.values())
sum_row_sq = sum(v * v for v in row_totals)

denominator = k * T - sum_row_sq
if denominator <= 0:
    raise SystemExit("Cochran Q undefined: insufficient outcome variation.")

Q = (k - 1) * (k * sum_col_sq - T * T) / denominator
p_Q = chi2.sf(Q, df=k - 1)

with OUT_Q.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["n_orthogroups", "n_species", "Q", "df", "p_value"])
    writer.writerow([N, k, Q, k - 1, p_Q])

# Exact two-sided McNemar tests: exact binomial test among discordant pairs.
pair_results = []
for a, b in combinations(species, 2):
    a_only = 0
    b_only = 0
    both = 0
    neither = 0

    for members in groups.values():
        av = members[a]["positive"]
        bv = members[b]["positive"]
        if av and bv:
            both += 1
        elif av:
            a_only += 1
        elif bv:
            b_only += 1
        else:
            neither += 1

    discordant = a_only + b_only
    p = (
        binomtest(a_only, n=discordant, p=0.5, alternative="two-sided").pvalue
        if discordant else 1.0
    )

    pair_results.append({
        "species_A": a,
        "species_B": b,
        "both_positive": both,
        "A_only": a_only,
        "B_only": b_only,
        "neither": neither,
        "discordant_total": discordant,
        "prevalence_difference_percentage_points":
            100 * (a_only - b_only) / N,
        "mcnemar_exact_p": p,
    })

# Benjamini-Hochberg correction across all 15 pairwise tests.
order = sorted(
    range(len(pair_results)),
    key=lambda i: pair_results[i]["mcnemar_exact_p"]
)
m = len(pair_results)
adjusted = [1.0] * m
running = 1.0

for j in range(m - 1, -1, -1):
    i = order[j]
    rank = j + 1
    running = min(running, pair_results[i]["mcnemar_exact_p"] * m / rank)
    adjusted[i] = min(1.0, running)

for i, result in enumerate(pair_results):
    result["BH_FDR"] = adjusted[i]

with OUT_PAIR.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(pair_results[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(pair_results)

print("COCHRAN'S Q: MATCHED SIX-SPECIES TEST")
print(f"Orthogroups: {N}")
print(f"Q = {Q:.6f}")
print(f"df = {k-1}")
print(f"p = {p_Q:.6g}")

print("\nEXACT PAIRWISE McNEMAR TESTS")
print("Species_A\tSpecies_B\tA_only\tB_only\tDifference_pp\tExact_p\tBH_FDR")
for r in sorted(pair_results, key=lambda x: x["mcnemar_exact_p"]):
    print(
        f"{r['species_A']}\t{r['species_B']}\t"
        f"{r['A_only']}\t{r['B_only']}\t"
        f"{r['prevalence_difference_percentage_points']:.2f}\t"
        f"{r['mcnemar_exact_p']:.6g}\t{r['BH_FDR']:.6g}"
    )

print("\nSaved:")
print(OUT_Q)
print(OUT_PAIR)
