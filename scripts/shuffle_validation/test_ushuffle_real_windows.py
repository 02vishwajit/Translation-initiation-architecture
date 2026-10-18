from collections import Counter
from Bio import SeqIO
import ushuffle
import random

FASTA = (
    "14_expanded_study/features/single_copy_upstream/"
    "single_copy_upstream50.fasta"
)
SEED = 20261024
SHUFFLES_PER_SEQUENCE = 10

def dinuc(s):
    return Counter(s[i:i+2] for i in range(len(s)-1))

random.seed(SEED)
ushuffle.set_seed(SEED)

records = list(SeqIO.parse(FASTA, "fasta"))
print("FASTA:", FASTA)
print("Sequences loaded:", len(records))

tested = 0
invalid = 0
unchanged = 0
errors = []

for record in records:
    seq = str(record.seq).upper()

    if len(seq) != 50 or set(seq) - set("ACGT"):
        continue

    target = dinuc(seq)

    for _ in range(SHUFFLES_PER_SEQUENCE):
        try:
            result = ushuffle.shuffle(seq.encode(), 2)
            if isinstance(result, bytes):
                result = result.decode().upper()
            else:
                result = str(result).upper()

            tested += 1

            if len(result) != len(seq) or dinuc(result) != target:
                invalid += 1

            if result == seq:
                unchanged += 1

        except Exception as exc:
            errors.append((record.id, repr(exc)))
            break

print("\n--- Validation summary ---")
print("Shuffle attempts:", tested)
print("Invalid outputs:", invalid)
print("Unchanged outputs:", unchanged)
print("Errors:", errors[:5])
print("Dinucleotide preservation passed:", invalid == 0 and not errors)
