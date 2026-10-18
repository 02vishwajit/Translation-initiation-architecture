from collections import Counter
from itertools import permutations
import ushuffle
import random

SEQ = "TTCGAGCTGG"
N = 100_000
SEED = 20261023

def dinuc(s):
    return Counter(s[i:i+2] for i in range(len(s)-1))

target = dinuc(SEQ)

# Exact reference set for this short test sequence.
valid = {
    "".join(p)
    for p in permutations(SEQ)
    if dinuc("".join(p)) == target
}

print("Exact valid arrangements:", len(valid))
print("ushuffle module:", ushuffle.__file__)

# Inspect the callable signatures/docstrings.
for name in ("shuffle", "Shuffler", "set_seed"):
    obj = getattr(ushuffle, name)
    print(f"\n{name}:")
    print("  Doc:", getattr(obj, "__doc__", None))

# Try the documented function interface.
random.seed(SEED)
ushuffle.set_seed(SEED)

counts = Counter()
invalid = 0
errors = []

for i in range(N):
    try:
        result = ushuffle.shuffle(SEQ.encode(), 2)
        if isinstance(result, bytes):
            result = result.decode()
        else:
            result = str(result)

        if len(result) != len(SEQ) or dinuc(result) != target:
            invalid += 1
        else:
            counts[result] += 1
    except Exception as exc:
        errors.append(repr(exc))
        break

print("\n--- Results ---")
print("Requested shuffles:", N)
print("Distinct valid arrangements sampled:", len(counts))
print("Invalid outputs:", invalid)
print("Errors:", errors[:3])
if counts:
    print("Minimum frequency:", min(counts.values()))
    print("Maximum frequency:", max(counts.values()))
    print("Original sampled:", SEQ in counts)
    print("States not sampled:", len(valid - set(counts)))
    print("\nTen most frequent outputs:")
    for seq, n in counts.most_common(10):
        print(f"{seq}\t{n}\t{n/N:.5%}")
