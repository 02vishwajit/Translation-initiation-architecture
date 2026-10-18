
#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path

INPUT = Path(
    "14_expanded_study/features/ccm_upstream/"
    "five_CCM_sequence_features.tsv"
)

OUTPUT = Path(
    "14_expanded_study/features/ccm_upstream/"
    "five_CCM_motif_refined_v2.tsv"
)

MOTIFS = (
    "AGGAGG",
    "GGAGG",
    "AGGAG",
    "GAGG",
    "AGGA",
    "GGAG",
)

MIN_SPACER = 4
MAX_SPACER = 12
TARGET_SPACER = 7


def scan_motifs(sequence):
    """Find all overlapping motif occurrences independently."""
    sequence = sequence.upper()
    n = len(sequence)
    hits = []

    for motif in MOTIFS:
        for match in re.finditer(rf"(?={re.escape(motif)})", sequence):
            start = match.start()
            end = start + len(motif)

            hits.append({
                "motif": motif,
                "start_1based": start + 1,
                "end_1based": end,
                "spacer_nt": n - end,
                "length": len(motif),
            })

    return hits


def main():
    if not INPUT.is_file():
        sys.exit(f"ERROR: Input file not found: {INPUT}")

    results = []

    with INPUT.open(
        "r", newline="", encoding="utf-8-sig"
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            sys.exit("ERROR: Input file has no header.")

        required = {
            "orthogroup",
            "species",
            "upstream_sequence_5to3",
            "overlap_flag",
        }

        missing = required - set(reader.fieldnames)
        if missing:
            sys.exit(
                "ERROR: Missing columns: "
                + ", ".join(sorted(missing))
            )

        original_fields = list(reader.fieldnames)

        for row_number, row in enumerate(reader, start=2):
            sequence = (
                row.get("upstream_sequence_5to3") or ""
            ).strip().upper()

            if not sequence:
                sys.exit(
                    f"ERROR: Empty sequence at input line {row_number}"
                )

            hits = scan_motifs(sequence)

            # Filter each hit independently. Do not collapse overlapping
            # hits before applying the spacer criterion.
            valid = [
                h for h in hits
                if MIN_SPACER <= h["spacer_nt"] <= MAX_SPACER
            ]

            # Deterministic representative selection:
            # spacer nearest 7, longest motif, earliest position,
            # then alphabetical motif.
            ranked = sorted(
                valid,
                key=lambda h: (
                    abs(h["spacer_nt"] - TARGET_SPACER),
                    -h["length"],
                    h["start_1based"],
                    h["motif"],
                ),
            )

            best = ranked[0] if ranked else None

            # Deterministic order for the full candidate list.
            listed = sorted(
                valid,
                key=lambda h: (
                    h["start_1based"],
                    -h["length"],
                    h["motif"],
                ),
            )

            row["raw_motif_match_count"] = len(hits)
            row["valid_spacer_hit_count"] = len(valid)
            row["candidate_hits_spacer_4_12"] = ";".join(
                f'{h["motif"]}@{h["start_1based"]}:'
                f'spacer{h["spacer_nt"]}'
                for h in listed
            )
            row["has_candidate_spacer_4_12"] = bool(valid)

            row["best_candidate_motif"] = (
                best["motif"] if best else ""
            )
            row["best_candidate_spacer_nt"] = (
                best["spacer_nt"] if best else ""
            )
            row["best_candidate_start_1based"] = (
                best["start_1based"] if best else ""
            )

            # This is the source window's overlap annotation, not proof
            # that the selected motif itself overlaps a CDS.
            row["best_candidate_overlaps_CDS"] = (
                row.get("overlap_flag", "") if best else ""
            )

            results.append(row)

    if not results:
        sys.exit("ERROR: No input records were processed.")

    extra_fields = [
        "raw_motif_match_count",
        "valid_spacer_hit_count",
        "candidate_hits_spacer_4_12",
        "has_candidate_spacer_4_12",
        "best_candidate_motif",
        "best_candidate_spacer_nt",
        "best_candidate_start_1based",
        "best_candidate_overlaps_CDS",
    ]

    fields = list(original_fields)
    for field in extra_fields:
        if field not in fields:
            fields.append(field)

    # Write only to the new output path.
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(results)

    positives = sum(
        bool(row["has_candidate_spacer_4_12"])
        for row in results
    )

    print("Motif refinement completed.")
    print("Input:", INPUT)
    print("Output:", OUTPUT)
    print("Rows:", len(results))
    print("Positive rows:", positives)

    print("\nBy species:")
    species_names = sorted({
        row.get("species", "") for row in results
    })

    for species in species_names:
        group = [
            row for row in results
            if row.get("species", "") == species
        ]
        count = sum(
            bool(row["has_candidate_spacer_4_12"])
            for row in group
        )
        print(f"{species}: {count}/{len(group)}")

    print("\nAll valid hits are retained, including overlapping hits.")
    print("Candidates are heuristic, not experimentally validated.")


if __name__ == "__main__":
    main()
