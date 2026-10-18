#!/usr/bin/env python3

from collections import Counter, defaultdict
from itertools import product
import random

ALPHABET = "ACGT"
SEED = 20261022
N_TEST = 20000


def dinuc(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))


def shuffle(seq, rng):
    adj = defaultdict(list)
    for a, b in zip(seq[:-1], seq[1:]):
        adj[a].append(b)

    for node in adj:
        rng.shuffle(adj[node])

    local = {k: list(v) for k, v in adj.items()}
    stack = [seq[0]]
    path = []

    while stack:
        node = stack[-1]
        if local.get(node):
            stack.append(local[node].pop())
        else:
            path.append(stack.pop())

    out = "".join(reversed(path))

    assert len(out) == len(seq)
    assert Counter(out) == Counter(seq)
    assert dinuc(out) == dinuc(seq)
    return out


def enumerate_valid(seq):
    target = dinuc(seq)
    valid = []

    for letters in product(ALPHABET, repeat=len(seq)-1):
        candidate = seq[0] + "".join(letters)
        if dinuc(candidate) == target:
            valid.append(candidate)

    return sorted(set(valid))


def main():
    rng = random.Random(SEED)

    # Short examples allow exact enumeration.
    examples = ["AABA", "ATAT", "AAGAA", "ACGAC"]
    # Use DNA-only examples.
    examples = ["ATAT", "AAGAA", "ACGAC"]

    for seq in examples:
        valid = enumerate_valid(seq)
        if not valid:
            continue

        counts = Counter(shuffle(seq, rng) for _ in range(N_TEST))
        expected = N_TEST / len(valid)

        chi_like = sum(
            (counts.get(s, 0) - expected) ** 2 / expected
            for s in valid
        )

        print(f"\nInput: {seq}")
        print(f"Distinct valid sequences: {len(valid)}")
        print(f"Unique sequences sampled: {len(counts)}")
        print(f"Expected count if uniform: {expected:.2f}")
        print(f"Frequency statistic (descriptive): {chi_like:.3f}")

        print("Sequence\tObserved_count\tExpected_if_uniform")
        for s in valid:
            print(f"{s}\t{counts.get(s, 0)}\t{expected:.2f}")


if __name__ == "__main__":
    main()
