import csv
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
import ushuffle

BASE = Path("14_expanded_study/features/single_copy_upstream")
FASTA = BASE / "single_copy_upstream50.fasta"
TABLE = BASE / "single_copy_sequence_features_corrected.tsv"
OUT = BASE / "single_copy_ushuffle_nonoverlap_summary.tsv"

N_REPLICATES = 1000
SEED = 20261018
MOTIFS = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA", "GGAG")

def motif_positive(seq):
    for motif in MOTIFS:
        start = 0
        while True:
            pos = seq.find(motif, start)
            if pos < 0:
                break
            spacer = len(seq) - (pos + len(motif))
            if 4 <= spacer <= 12:
                return True
            start = pos + 1
    return False

def dinuc(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))

with open(TABLE, newline="") as f:
    table_rows = list(csv.DictReader(f, delimiter="\t"))

table = {
    (r["orthogroup"], r["species"]): r
    for r in table_rows
}

by_species = defaultdict(list)

for rec in SeqIO.parse(FASTA, "fasta"):
    parts = rec.id.split("|")
    key = (parts[0], parts[1])
    if key not in table:
        raise SystemExit(f"Missing feature row: {key}")

    row = table[key]
    overlap = row["overlap_flag"].strip().lower() == "true"
    if overlap:
        continue

    seq = str(rec.seq).upper()
    if len(seq) != 50 or set(seq) - set("ACGT"):
        raise SystemExit(f"Invalid sequence: {rec.id}")

    observed = motif_positive(seq)
    recorded = row["has_candidate_spacer_4_12"].strip().lower() == "true"
    if observed != recorded:
        raise SystemExit(f"Motif scanner mismatch: {rec.id}")

    by_species[parts[1]].append((rec.id, seq, observed))

if not by_species:
    raise SystemExit("No non-overlapping sequences found.")

ushuffle.set_seed(SEED)
summary = []

for species in sorted(by_species):
    items = by_species[species]
    n = len(items)
    observed = sum(x[2] for x in items)
    null_counts = []
    invalid = 0

    for rep in range(N_REPLICATES):
        positive = 0
        for record_id, seq, original_positive in items:
            shuffled = ushuffle.shuffle(seq.encode("ascii"), 2).decode("ascii")

            if (
                len(shuffled) != len(seq)
                or shuffled[0] != seq[0]
                or shuffled[-1] != seq[-1]
                or dinuc(shuffled) != dinuc(seq)
            ):
                invalid += 1
                continue

            positive += int(motif_positive(shuffled))

        if invalid:
            raise SystemExit(f"Invalid shuffle detected for {species}")
        null_counts.append(positive)

    null_mean = sum(null_counts) / N_REPLICATES
    exceedances = sum(x >= observed for x in null_counts)
    p_empirical = (exceedances + 1) / (N_REPLICATES + 1)
    sorted_counts = sorted(null_counts)

    summary.append({
        "species": species,
        "n_nonoverlap": n,
        "observed_positive": observed,
        "observed_percent": 100 * observed / n,
        "null_mean_positive": null_mean,
        "null_mean_percent": 100 * null_mean / n,
        "enrichment_percentage_points": 100 * (observed-null_mean) / n,
        "null_2_5_percentile_positive": sorted_counts[int(0.025*N_REPLICATES)],
        "null_97_5_percentile_positive": sorted_counts[int(0.975*N_REPLICATES)],
        "empirical_upper_tail_p": p_empirical,
        "replicates_at_least_observed": exceedances,
        "replicates": N_REPLICATES,
        "invalid_shuffles": invalid,
    })

# Benjamini-Hochberg correction across the six species.
order = sorted(
    range(len(summary)),
    key=lambda i: summary[i]["empirical_upper_tail_p"]
)
m = len(summary)
adjusted = [1.0] * m
running = 1.0

for j in range(m-1, -1, -1):
    i = order[j]
    running = min(running, summary[i]["empirical_upper_tail_p"] * m / (j+1))
    adjusted[i] = min(1.0, running)

for i, row in enumerate(summary):
    row["BH_FDR"] = adjusted[i]

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary[0]), delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(summary)

print("NON-OVERLAPPING WINDOWS: OBSERVED VS DINUCLEOTIDE NULL")
for r in summary:
    print(
        f"{r['species']}: n={r['n_nonoverlap']}, "
        f"observed={r['observed_percent']:.2f}%, "
        f"null={r['null_mean_percent']:.2f}%, "
        f"enrichment={r['enrichment_percentage_points']:.2f} pp, "
        f"p={r['empirical_upper_tail_p']:.4g}, "
        f"FDR={r['BH_FDR']:.4g}"
    )

print("\nSaved:", OUT)
