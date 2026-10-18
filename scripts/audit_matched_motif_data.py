import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

base = Path("14_expanded_study/features/single_copy_upstream")
input_file = base / "single_copy_sequence_features_corrected.tsv"
output_file = base / "matched_motif_species_diagnostics.tsv"

with input_file.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

by_species = defaultdict(list)
by_og = defaultdict(dict)

for r in rows:
    species = r["species"]
    og = r["orthogroup"]
    positive = r["has_candidate_spacer_4_12"].strip().lower() == "true"
    overlap = r["overlap_flag"].strip().lower() == "true"

    try:
        gc = float(r["GC_percent"])
    except (ValueError, TypeError):
        raise SystemExit(f"Invalid GC value for {og}/{species}")

    codon = r["start_codon_class"].strip()
    item = {
        "og": og, "species": species, "positive": positive,
        "overlap": overlap, "gc": gc, "codon": codon
    }
    by_species[species].append(item)
    by_og[og][species] = item

species_list = sorted(by_species)
if len(by_og) != 348 or any(len(v) != 6 for v in by_og.values()):
    raise SystemExit("Expected 348 orthogroups with exactly six species each.")

print("SPECIES-LEVEL DIAGNOSTICS")
print("Species\tn\tPositive\tPrevalence_pct\tMean_GC_pct\tOverlap_pct")
diagnostics = []

for species in species_list:
    items = by_species[species]
    n = len(items)
    positives = sum(x["positive"] for x in items)
    row = {
        "species": species,
        "n": n,
        "positive": positives,
        "prevalence_pct": 100 * positives / n,
        "mean_GC_pct": mean(x["gc"] for x in items),
        "overlap_pct": 100 * sum(x["overlap"] for x in items) / n,
    }
    diagnostics.append(row)
    print(
        f"{species}\t{n}\t{positives}\t"
        f"{row['prevalence_pct']:.2f}\t{row['mean_GC_pct']:.2f}\t"
        f"{row['overlap_pct']:.2f}"
    )

print("\nSTART-CODON COUNTS BY SPECIES")
for species in species_list:
    counts = Counter(x["codon"] for x in by_species[species])
    print(species, dict(sorted(counts.items())))

print("\nORTHOGROUP MATCHING")
pattern_counts = Counter()
for og, members in by_og.items():
    # Number of species with a candidate motif in this orthogroup.
    pattern_counts[sum(x["positive"] for x in members.values())] += 1

for positives, count in sorted(pattern_counts.items()):
    print(f"{positives}/6 species positive: {count} orthogroups")

print("\nSpecies pairs: within-orthogroup matched differences")
print("Species_A\tSpecies_B\tBoth_observed\tA_only\tB_only\tNeither")
pair_rows = []

for i, a in enumerate(species_list):
    for b in species_list[i+1:]:
        both = a_only = b_only = neither = 0
        for members in by_og.values():
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
        pair_rows.append({
            "species_A": a, "species_B": b,
            "both_positive": both, "A_only": a_only,
            "B_only": b_only, "neither": neither
        })
        print(f"{a}\t{b}\t{both}\t{a_only}\t{b_only}\t{neither}")

with output_file.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(diagnostics[0]), delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(diagnostics)

print("\nSaved:", output_file)
print("Audit completed.")
