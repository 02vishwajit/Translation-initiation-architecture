from collections import Counter
from itertools import permutations
import importlib.util
import ushuffle

SEQ = "TTCGAGCTGG"
N = 100_000
SEEDS = [101, 202, 303, 404, 505]

def dinuc(s):
    return Counter(s[i:i+2] for i in range(len(s)-1))

target = dinuc(SEQ)

valid = {
    "".join(p)
    for p in permutations(SEQ)
    if dinuc("".join(p)) == target
}

print("Exact state space:", len(valid))
print("Expected count per state:", N / len(valid))
print("SciPy available:", importlib.util.find_spec("scipy") is not None)

if importlib.util.find_spec("scipy") is not None:
    from scipy.stats import chisquare

for seed in SEEDS:
    ushuffle.set_seed(seed)
    counts = Counter()
    invalid = 0

    for _ in range(N):
        result = ushuffle.shuffle(SEQ.encode(), 2)
        if isinstance(result, bytes):
            result = result.decode()
        else:
            result = str(result)

        if len(result) != len(SEQ) or dinuc(result) != target:
            invalid += 1
        else:
            counts[result] += 1

    frequencies = [counts[s] for s in sorted(valid)]
    print(f"\nSeed: {seed}")
    print("Distinct states:", len(counts))
    print("Invalid outputs:", invalid)
    print("Minimum frequency:", min(frequencies))
    print("Maximum frequency:", max(frequencies))

    if importlib.util.find_spec("scipy") is not None:
        result = chisquare(frequencies)
        print("Chi-square statistic:", round(result.statistic, 3))
        print("Chi-square p-value:", f"{result.pvalue:.6g}")
