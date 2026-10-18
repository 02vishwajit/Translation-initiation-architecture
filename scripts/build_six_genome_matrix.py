from pathlib import Path
from itertools import combinations
import csv

BASE = Path("14_expanded_study/orthology/pairwise_rbh")
OUT = Path("14_expanded_study/orthology/six_genome_matrix")
OUT.mkdir(parents=True, exist_ok=True)

STRAINS = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

# Working filters; retain the original RBH files for sensitivity analyses.
MIN_IDENTITY = 50.0
MIN_COVERAGE = 0.70

def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

# Pairwise mapping: (strain, protein) -> protein in another strain.
partner = {}
edge_stats = {}
pairwise_rows = []

files = sorted(BASE.glob("*__*__RBH.tsv"))
if len(files) != 15:
    raise SystemExit(
        f"Expected 15 pairwise RBH files; found {len(files)} in {BASE}"
    )

for path in files:
    rows = read_tsv(path)
    if not rows:
        print(f"WARNING: empty RBH file: {path}")
        continue

    for r in rows:
        a, b = r["strain_a"], r["strain_b"]
        pa, pb = r["protein_a"], r["protein_b"]
        identity = float(r["identity_percent"])
        ca, cb = float(r["coverage_a"]), float(r["coverage_b"])

        accepted = (
            identity >= MIN_IDENTITY
            and ca >= MIN_COVERAGE
            and cb >= MIN_COVERAGE
        )

        pairwise_rows.append({
            "strain_a": a, "protein_a": pa,
            "strain_b": b, "protein_b": pb,
            "identity_percent": identity,
            "coverage_a": ca, "coverage_b": cb,
            "passes_working_filters": int(accepted),
            "source_file": path.name
        })

        if not accepted:
            continue

        key_ab = ((a, pa), b)
        key_ba = ((b, pb), a)

        # A protein should have only one RBH per other genome.
        # Flag conflicts instead of silently overwriting them.
        for key, value in [(key_ab, pb), (key_ba, pa)]:
            if key in partner and partner[key] != value:
                raise SystemExit(
                    f"Conflicting RBH mapping for {key}: "
                    f"{partner[key]} vs {value}"
                )
            partner[key] = value

        edge_stats[((a, pa), (b, pb))] = r
        edge_stats[((b, pb), (a, pa))] = r

with (OUT / "pairwise_rbh_filter_audit.tsv").open("w", newline="") as f:
    fields = [
        "strain_a", "protein_a", "strain_b", "protein_b",
        "identity_percent", "coverage_a", "coverage_b",
        "passes_working_filters", "source_file"
    ]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(pairwise_rows)

# Find groups connected across all six genomes.
# For each seed protein, follow its filtered RBH into every other strain,
# then require every pair of proteins in the group to have a filtered RBH.
groups = set()

all_seed_proteins = set()
for (strain, protein), other_strain in partner:
    all_seed_proteins.add((strain, protein))

for seed_strain, seed_protein in sorted(all_seed_proteins):
    group = {seed_strain: seed_protein}
    valid = True

    for target in STRAINS:
        if target == seed_strain:
            continue
        mapped = partner.get(((seed_strain, seed_protein), target))
        if mapped is None:
            valid = False
            break
        group[target] = mapped

    if not valid or len(group) != len(STRAINS):
        continue

    # Require pairwise consistency for all 15 genome pairs.
    for s1, s2 in combinations(STRAINS, 2):
        if ((s1, group[s1]), (s2, group[s2])) not in edge_stats:
            valid = False
            break

    if valid:
        groups.add(tuple(group[s] for s in STRAINS))

matrix_path = OUT / "strict_six_way_ortholog_matrix.tsv"
with matrix_path.open("w", newline="") as f:
    fields = ["ortholog_group"] + STRAINS
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for i, group in enumerate(sorted(groups), start=1):
        row = {"ortholog_group": f"SIXWAY_{i:05d}"}
        row.update(dict(zip(STRAINS, group)))
        w.writerow(row)

# Pairwise filtered counts and fraction retained.
summary_path = OUT / "pairwise_filter_summary.tsv"
with summary_path.open("w", newline="") as f:
    fields = [
        "strain_a", "strain_b", "raw_rbh_pairs",
        "filtered_rbh_pairs", "identity_cutoff",
        "coverage_cutoff_each_protein"
    ]
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for a, b in combinations(STRAINS, 2):
        subset = [
            r for r in pairwise_rows
            if {r["strain_a"], r["strain_b"]} == {a, b}
        ]
        passed = sum(r["passes_working_filters"] for r in subset)
        w.writerow({
            "strain_a": a,
            "strain_b": b,
            "raw_rbh_pairs": len(subset),
            "filtered_rbh_pairs": passed,
            "identity_cutoff": MIN_IDENTITY,
            "coverage_cutoff_each_protein": MIN_COVERAGE
        })

print(f"Pairwise filter audit: {OUT / 'pairwise_rbh_filter_audit.tsv'}")
print(f"Pairwise filter summary: {summary_path}")
print(f"Strict six-way matrix: {matrix_path}")
print(f"Strict six-way groups: {len(groups)}")
print("Filters: identity >= 50%; coverage >= 70% for both proteins")
