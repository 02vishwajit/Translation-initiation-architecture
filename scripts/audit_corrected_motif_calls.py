from pathlib import Path
import re
import pandas as pd

base = Path("14_expanded_study/features/single_copy_upstream")
source = base / "single_copy_upstream50.tsv"
corrected = base / "single_copy_sequence_features_corrected.tsv"

motifs = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "AGGA", "GGAG")
minimum_spacer, maximum_spacer = 4, 12

def expected_positive(sequence):
    sequence = sequence.upper()
    for motif in motifs:
        for match in re.finditer(f"(?={motif})", sequence):
            spacer = len(sequence) - (match.start() + len(motif))
            if minimum_spacer <= spacer <= maximum_spacer:
                return True
    return False

src = pd.read_csv(source, sep="\t", dtype=str)
obs = pd.read_csv(corrected, sep="\t", dtype=str)

key = ["orthogroup", "species", "protein_id", "locus_tag"]

if src.duplicated(key).any() or obs.duplicated(key).any():
    raise SystemExit("Duplicate keys detected.")

merged = src.merge(
    obs[key + ["has_candidate_spacer_4_12"]],
    on=key,
    how="outer",
    indicator=True,
    validate="one_to_one",
)

print("Source rows:", len(src))
print("Corrected rows:", len(obs))
print("Join status:")
print(merged["_merge"].value_counts().to_string())

if not merged["_merge"].eq("both").all():
    raise SystemExit("Source/corrected key mismatch.")

# The stored upstream sequence must be unchanged.
check = src.merge(
    obs[key + ["upstream_sequence_5to3"]],
    on=key,
    suffixes=("_source", "_corrected"),
    validate="one_to_one",
)
sequence_equal = (
    check["upstream_sequence_5to3_source"]
    == check["upstream_sequence_5to3_corrected"]
)
print("Identical source/corrected sequences:",
      int(sequence_equal.sum()), "/", len(check))

# Attach the source sequence directly; merged contains only the key
# columns and the observed Boolean motif call.
expected = merged.merge(
    src[key + ["upstream_sequence_5to3"]],
    on=key,
    validate="one_to_one",
    suffixes=("", "_source"),
)
expected["expected_positive"] = expected[
    "upstream_sequence_5to3"
].map(expected_positive)

expected["observed_positive"] = (
    expected["has_candidate_spacer_4_12"]
    .str.lower().eq("true")
)

mismatch = expected[
    expected["expected_positive"] != expected["observed_positive"]
]

print("Independent motif-call mismatches:", len(mismatch))
if len(mismatch):
    print(mismatch[
        key + ["expected_positive", "observed_positive"]
    ].head(30).to_string(index=False))

print("\nExpected positives by species:")
print(expected.groupby("species")["expected_positive"].agg(
    ["sum", "count", "mean"]
).to_string())

out = Path("14_expanded_study/results/corrected_motif_call_audit.tsv")
expected[key + [
    "expected_positive", "observed_positive"
]].to_csv(out, sep="\t", index=False)
print("\nAudit report:", out)

if not sequence_equal.all():
    raise SystemExit("Source and corrected sequences differ.")
if len(mismatch):
    raise SystemExit("Corrected motif calls differ from independent scanner.")
