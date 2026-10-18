#!/usr/bin/env python3
"""Reproducible independent-hit scanner for upstream sequence candidates.

Input:
  features/single_copy_upstream/single_copy_upstream50.tsv

Output:
  features/single_copy_upstream/single_copy_sequence_features_corrected_rebuilt.tsv

Every occurrence of each specified motif is assessed independently.
This is a sequence-feature heuristic, not experimental validation of an RBS.
"""

import csv
import re
from collections import Counter
from pathlib import Path

BASE = Path("14_expanded_study/features/single_copy_upstream")
INPUT = BASE / "single_copy_upstream50.tsv"
OUTPUT = BASE / "single_copy_sequence_features_corrected_rebuilt.tsv"

MOTIFS = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA", "GGAG")
MIN_SPACER = 4
MAX_SPACER = 12
TARGET_SPACER = 7


def find_hits(sequence):
    """Return all motif occurrences, including overlapping occurrences."""
    hits = []
    sequence = sequence.upper()

    for motif in MOTIFS:
        for match in re.finditer(f"(?={re.escape(motif)})", sequence):
            start = match.start()
            end = start + len(motif)
            hits.append({
                "motif": motif,
                "start": start,
                "end": end,
                "length": len(motif),
                "spacer": len(sequence) - end,
            })

    return hits


def main():
    if not INPUT.is_file():
        raise SystemExit(f"Missing input file: {INPUT}")

    with INPUT.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = reader.fieldnames
        if not input_fields:
            raise SystemExit("Input table has no header.")
        rows = list(reader)

    if not rows:
        raise SystemExit("Input table is empty.")

    required = {
        "upstream_sequence_5to3",
        "orthogroup",
        "species",
        "overlapping_CDS",
        "start_codon",
    }
    missing = required.difference(input_fields)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    output_rows = []

    for row in rows:
        sequence = row["upstream_sequence_5to3"].upper()
        if not sequence:
            raise SystemExit(
                f"Empty sequence: {row.get('orthogroup')} "
                f"{row.get('species')}"
            )

        gc_count = sequence.count("G") + sequence.count("C")
        purine_count = sequence.count("A") + sequence.count("G")

        hits = find_hits(sequence)
        plausible = [
            hit for hit in hits
            if MIN_SPACER <= hit["spacer"] <= MAX_SPACER
        ]

        plausible.sort(key=lambda hit: (
            abs(hit["spacer"] - TARGET_SPACER),
            -hit["length"],
            hit["start"],
            hit["motif"],
        ))

        best = plausible[0] if plausible else None

        start_codon = row["start_codon"].upper()
        start_class = (
            "ATG" if start_codon == "ATG"
            else "Alternative_start"
            if start_codon in {"GTG", "TTG"}
            else "Other_or_unexpected"
        )

        # Preserve all original input fields, then set derived fields.
        row.update({
            "GC_percent": round(100 * gc_count / len(sequence), 3),
            "purine_percent": round(100 * purine_count / len(sequence), 3),
            "raw_motif_match_count": len(hits),
            "overlap_collapsed_candidate_count": len(hits),
            "candidate_motifs_spacer_4_12": ";".join(
                f'{hit["motif"]}@{hit["start"] + 1}:spacer{hit["spacer"]}'
                for hit in plausible
            ),
            "has_candidate_spacer_4_12": bool(plausible),
            "best_candidate_motif": best["motif"] if best else "",
            "best_candidate_spacer_nt": best["spacer"] if best else "",
            "best_candidate_start_1based": best["start"] + 1 if best else "",
            # A non-empty neighboring CDS locus tag indicates overlap.
            # Empty/NaN values indicate no overlap.
            "overlap_flag": bool(
                str(row["overlapping_CDS"]).strip()
                and str(row["overlapping_CDS"]).strip().lower()
                    not in {"nan", "none", "null"}
            ),
            "start_codon_class": start_class,
        })
        output_rows.append(row)

    # Keep original columns first; append any newly derived columns.
    fields = list(input_fields)
    for row in output_rows:
        for field in row:
            if field not in fields:
                fields.append(field)

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            lineterminator="\n", extrasaction="raise"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print("Input:", INPUT)
    print("Output:", OUTPUT)
    print("Rows:", len(output_rows))
    print("Unique orthogroups:",
          len({row["orthogroup"] for row in output_rows}))
    print("Sequence lengths:",
          dict(Counter(len(row["upstream_sequence_5to3"])
                       for row in output_rows)))
    print("Motif-positive:",
          sum(bool(row["has_candidate_spacer_4_12"])
              for row in output_rows))
    print("Overlapping CDS windows:",
          sum(bool(row["overlap_flag"]) for row in output_rows))

    print("\nMotif-positive by species:")
    species = sorted({row["species"] for row in output_rows})
    for name in species:
        group = [row for row in output_rows if row["species"] == name]
        positive = sum(bool(row["has_candidate_spacer_4_12"])
                       for row in group)
        print(f"{name}: {positive}/{len(group)} "
              f"({100 * positive / len(group):.2f}%)")


if __name__ == "__main__":
    main()
