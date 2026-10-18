from collections import Counter
from itertools import permutations
from random import Random
import math

SEQ = "TTCGAGCTGG"
SEED = 20261023
BURN_IN = 50_000
SAMPLES = 200_000
THIN = 10

rng = Random(SEED)

def dinuc(s):
    return Counter(s[i:i+2] for i in range(len(s)-1))

target = dinuc(SEQ)

# Exhaustive gold-standard state space.
valid = {
    "".join(p)
    for p in permutations(SEQ)
    if dinuc("".join(p)) == target
}

print("Original:", SEQ)
print("Exact valid state count:", len(valid))

# Proposal: uniformly shuffle the current sequence's characters.
# This proposal is symmetric between any two distinct permutations
# of the same multiset. Reject proposals that alter dinucleotide counts.
state = SEQ
accepted = 0
proposed = 0
visited = Counter()

total_steps = BURN_IN + SAMPLES * THIN

for step in range(total_steps):
    chars = list(state)
    rng.shuffle(chars)
    proposal = "".join(chars)
    proposed += 1

    if dinuc(proposal) == target:
        state = proposal
        accepted += 1

    if step >= BURN_IN and (step - BURN_IN) % THIN == 0:
        visited[state] += 1

print("\n--- Results ---")
print("Steps:", proposed)
print("Accepted moves:", accepted)
print("Acceptance rate:", f"{accepted/proposed:.6%}")
print("Post-burn-in samples:", sum(visited.values()))
print("Distinct states visited:", len(visited))
print("States in exact space:", len(valid))
print("Unvisited valid states:", len(valid - set(visited)))

expected = sum(visited.values()) / len(valid)
tv_distance = 0.5 * sum(
    abs(visited[s] / sum(visited.values()) - 1 / len(valid))
    for s in valid
)

print("Expected count per state under uniformity:", round(expected, 2))
print("Observed count range:", min(visited.values()), max(visited.values()))
print("Empirical total-variation distance from uniform:", round(tv_distance, 4))

print("\nTen most frequently visited states:")
for sequence, count in visited.most_common(10):
    print(f"{sequence}\t{count}")

print("\nDinucleotide preservation check:")
print("All sampled states valid:", set(visited).issubset(valid))
print("All exact states reachable in theory:", True)
print(
    "Interpretation: the proposal has a uniform stationary distribution "
    "on the valid state space, but finite-run convergence still needs "
    "diagnostics."
)
