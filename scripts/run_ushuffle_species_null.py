import csv
import random
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
import ushuffle

BASE = Path("14_expanded_study/features/single_copy_upstream")
FASTA = BASE / "single_copy_upstream50.fasta"
FEATURES = BASE / "single_copy_sequence_features_corrected.tsv"

OUT_REPLICATES = BASE / "single_copy_ushuffle_species_replicates.tsv"
OUT_SUMMARY = BASE / "single_copy_ushuffle_species_summary.tsv"
OUT_VALIDATION = BASE / "single_copy_ushuffle_validation.tsv"

N_REPLICATES = 1000
SEED = 20261017
MOTIFS = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA", "GGAG")

random.seed(SEED)
ushuffle.set_seed(SEED)

def dinucleotide_counts(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))

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

records = list(SeqIO.parse(FASTA, "fasta"))
with open(FEATURES, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    feature_rows = list(reader)

feature_by_key = {
    (r["orthogroup"], r["species"]): r for r in feature_rows
}

if len(records) != 2088 or len(feature_rows) != 2088:
    raise SystemExit("STOP: Expected 2,088 FASTA and feature-table records.")

sequences = []
observed_by_species = Counter()
species_sequences = defaultdict(list)
validation = []

for record in records:
    parts = record.id.split("|")
    if len(parts) < 2:
        raise SystemExit(f"Unrecognized FASTA header: {record.id}")

    og, species = parts[0], parts[1]
    key = (og, species)
    if key not in feature_by_key:
        raise SystemExit(f"Missing feature-table row: {key}")

    seq = str(record.seq).upper()
    if len(seq) != 50 or set(seq) - set("ACGT"):
        raise SystemExit(f"Unexpected sequence length/alphabet: {record.id}")

    observed = motif_positive(seq)
    recorded = (
        feature_by_key[key]["has_candidate_spacer_4_12"]
        .strip().lower() == "true"
    )
    if observed != recorded:
        raise SystemExit(f"Motif classification mismatch: {record.id}")

    observed_by_species[species] += int(observed)
    species_sequences[species].append((record.id, seq, observed))
    sequences.append((record.id, species, seq, observed))

species_list = sorted(species_sequences)
if any(len(species_sequences[s]) != 348 for s in species_list):
    raise SystemExit("STOP: Each species must contain 348 sequences.")

replicate_rows = []
valid_rows = []
null_counts = {s: [] for s in species_list}

for rep in range(1, N_REPLICATES + 1):
    for species in species_list:
        positive_count = 0
        tested = 0
        failures = 0

        for record_id, seq, observed in species_sequences[species]:
            shuffled = ushuffle.shuffle(seq.encode("ascii"), 2).decode("ascii")
            tested += 1

            preserved = (
                len(shuffled) == len(seq)
                and shuffled[0] == seq[0]
                and shuffled[-1] == seq[-1]
                and dinucleotide_counts(shuffled) == dinucleotide_counts(seq)
            )

            if not preserved:
                failures += 1
                continue

            positive_count += int(motif_positive(shuffled))

        if failures:
            raise SystemExit(
                f"STOP: {failures} invalid shuffles in replicate {rep}, {species}"
            )

        null_counts[species].append(positive_count)
        replicate_rows.append({
            "replicate": rep,
            "species": species,
            "n_sequences": tested,
            "motif_positive": positive_count,
            "motif_prevalence_percent": 100 * positive_count / tested,
        })

    if rep % 50 == 0:
        print(f"Completed {rep}/{N_REPLICATES} randomized datasets.")

for species in species_list:
    valid_rows.append({
        "species": species,
        "sequences_tested_per_replicate": 348,
        "replicates": N_REPLICATES,
        "invalid_shuffles": 0,
        "dinucleotide_preservation": "PASS",
        "first_and_last_base_preservation": "PASS",
    })

with open(OUT_REPLICATES, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(replicate_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(replicate_rows)

with open(OUT_VALIDATION, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(valid_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(valid_rows)

summary_rows = []
for species in species_list:
    values = null_counts[species]
    observed = observed_by_species[species]
    exceedances = sum(x >= observed for x in values)
    mean_null = sum(values) / len(values)
    p_upper = (exceedances + 1) / (N_REPLICATES + 1)

    summary_rows.append({
        "species": species,
        "n_sequences": 348,
        "observed_positive": observed,
        "observed_prevalence_percent": 100 * observed / 348,
        "null_mean_positive": mean_null,
        "null_mean_prevalence_percent": 100 * mean_null / 348,
        "observed_minus_null_positive": observed - mean_null,
        "observed_minus_null_percentage_points":
            100 * (observed - mean_null) / 348,
        "null_2_5_percentile_positive":
            sorted(values)[int(0.025 * N_REPLICATES)],
        "null_97_5_percentile_positive":
            sorted(values)[int(0.975 * N_REPLICATES)],
        "empirical_upper_tail_p": p_upper,
        "replicates_at_least_observed": exceedances,
        "replicates": N_REPLICATES,
    })

# Benjamini-Hochberg adjustment across the six species.
order = sorted(
    range(len(summary_rows)),
    key=lambda i: summary_rows[i]["empirical_upper_tail_p"]
)
m = len(summary_rows)
adjusted = [1.0] * m
running = 1.0
for rank_index in range(m - 1, -1, -1):
    i = order[rank_index]
    rank = rank_index + 1
    p = summary_rows[i]["empirical_upper_tail_p"]
    running = min(running, p * m / rank)
    adjusted[i] = min(1.0, running)

for i, row in enumerate(summary_rows):
    row["BH_FDR"] = adjusted[i]

with open(OUT_SUMMARY, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(summary_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

print("\nOBSERVED VERSUS DINUCLEOTIDE-PRESERVING NULL")
for row in summary_rows:
    print(
        f"{row['species']}: observed "
        f"{row['observed_prevalence_percent']:.2f}%; null mean "
        f"{row['null_mean_prevalence_percent']:.2f}%; "
        f"difference "
        f"{row['observed_minus_null_percentage_points']:.2f} pp; "
        f"empirical p={row['empirical_upper_tail_p']:.4g}; "
        f"BH-FDR={row['BH_FDR']:.4g}"
    )

print("\nSaved:")
print(OUT_REPLICATES)
print(OUT_SUMMARY)
print(OUT_VALIDATION)
