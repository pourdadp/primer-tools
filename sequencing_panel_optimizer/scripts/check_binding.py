
import json

def reverse_complement(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(comp[base] for base in reversed(seq))

def are_bases_compatible(b1, b2):
    iupac_map = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
        'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'},
        'W': {'A','T'}, 'K': {'G','T'}, 'M': {'A','C'},
        'B': {'C','G','T'}, 'D': {'A','G','T'}, 'H': {'A','C','T'},
        'V': {'A','C','G'}, 'N': {'A','C','G','T'},
        'I': {'A','C','G','T'}
    }
    set1 = iupac_map.get(b1.upper(), {b1.upper()})
    set2 = iupac_map.get(b2.upper(), {b2.upper()})
    return bool(set1 & set2)

def semi_global_align(seq, primer, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    n = len(seq)
    m = len(primer)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    max_score = -float('inf')
    best_end_i = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if are_bases_compatible(seq[i-1], primer[j-1]):
                match = match_score
            else:
                match = mismatch_penalty
            diag = dp[i-1][j-1] + match
            up = dp[i-1][j] + gap_penalty
            left = dp[i][j-1] + gap_penalty
            dp[i][j] = max(diag, up, left, 0)
            if j == m and dp[i][j] > max_score:
                max_score = dp[i][j]
                best_end_i = i

    if max_score < 2:
        return None, None, None, None, None

    start_i = best_end_i
    start_j = m
    mismatches = 0
    gaps = 0

    while start_j > 0 and dp[start_i][start_j] > 0:
        if start_i > 0 and start_j > 0 and dp[start_i][start_j] == dp[start_i-1][start_j-1] + (match_score if are_bases_compatible(seq[start_i-1], primer[start_j-1]) else mismatch_penalty):
            if not are_bases_compatible(seq[start_i-1], primer[start_j-1]):
                mismatches += 1
            start_i -= 1
            start_j -= 1
        elif start_i > 0 and dp[start_i][start_j] == dp[start_i-1][start_j] + gap_penalty:
            gaps += 1
            start_i -= 1
        elif start_j > 0 and dp[start_i][start_j] == dp[start_i][start_j-1] + gap_penalty:
            gaps += 1
            start_j -= 1
        else:
            break

    return start_i, best_end_i, mismatches, gaps, max_score

def find_binding(sequence, primer, max_mismatch, max_gaps=0):
    positions = []
    start, end, mism, gaps, score = semi_global_align(sequence, primer)
    if start is not None and mism <= max_mismatch and gaps <= max_gaps:
        positions.append((start, end, mism, gaps, "forward", score))

    rc_primer = reverse_complement(primer)
    start, end, mism, gaps, score = semi_global_align(sequence, rc_primer)
    if start is not None and mism <= max_mismatch and gaps <= max_gaps:
        positions.append((start, end, mism, gaps, "reverse", score))

    return positions

with open("results/input_data.json", encoding='utf-8') as f:
    data = json.load(f)

binding_results = {}
rejected_primers = []
rejection_reasons = {}

for primer in data["primers"]:
    name = primer["name"]
    seq = primer["sequence"]
    positions = find_binding(data["sequence"], seq, data["max_mispairing"])

    forward_pos = [p for p in positions if p[4] == "forward"]
    reverse_pos = [p for p in positions if p[4] == "reverse"]

    reject_reason = None

    if forward_pos and reverse_pos:
        reject_reason = f"Binds in both directions (F:{len(forward_pos)} sites, R:{len(reverse_pos)} sites)"
    elif len(forward_pos) > 1:
        reject_reason = f"Binds to {len(forward_pos)} sites in Forward direction"
    elif len(reverse_pos) > 1:
        reject_reason = f"Binds to {len(reverse_pos)} sites in Reverse direction"
    elif not forward_pos and not reverse_pos:
        reject_reason = f"No binding site found (max_mismatch={data['max_mispairing']})"

    if reject_reason:
        rejected_primers.append(name)
        rejection_reasons[name] = reject_reason
        binding_results[name] = {
            "sequence": seq,
            "status": "REJECTED",
            "reason": reject_reason,
            "forward_positions": forward_pos,
            "reverse_positions": reverse_pos
        }
        print(f"WARNING: Primer '{name}' REJECTED: {reject_reason}")
    else:
        if forward_pos:
            binding_results[name] = {
                "sequence": seq,
                "status": "ACCEPTED",
                "direction": "forward",
                "positions": forward_pos
            }
        else:
            binding_results[name] = {
                "sequence": seq,
                "status": "ACCEPTED",
                "direction": "reverse",
                "positions": reverse_pos
            }
        print(f"INFO: Primer '{name}' ACCEPTED as {binding_results[name]['direction']}")

with open("results/binding_results.json", "w", encoding='utf-8') as out:
    json.dump(binding_results, out, indent=2)

accepted = sum(1 for v in binding_results.values() if v["status"] == "ACCEPTED")
print(f"\nSUMMARY: {accepted} primer(s) accepted, {len(rejected_primers)} rejected.")
if rejected_primers:
    print(f"Rejected: {', '.join(rejected_primers)}")
