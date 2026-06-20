
import json
import time
from tqdm import tqdm

def greedy_pair_coverage(valid_pairs, seq_len):
    intervals = []
    for pair_name, positions in valid_pairs.items():
        best_pos = max(positions, key=lambda x: x["product_length"])
        intervals.append({
            "pair": pair_name,
            "start": best_pos["f_start"],
            "end": best_pos["r_end"],
            "product_len": best_pos["product_length"],
            "tm_diff": best_pos["tm_diff"],
            "f_primer": pair_name.split('-')[0],
            "r_primer": pair_name.split('-')[1]
        })

    intervals.sort(key=lambda x: x["start"])
    selected = []
    covered = [False] * seq_len
    current_coverage = 0

    with tqdm(total=seq_len, desc="Covering target with amplicons", unit="bp") as pbar:
        pbar.update(0)
        while True:
            best_interval = None
            best_gain = 0
            for interval in intervals:
                if interval in selected:
                    continue
                start = interval["start"]
                end = interval["end"]
                new_coverage = 0
                for i in range(max(0, start), min(seq_len, end)):
                    if not covered[i]:
                        new_coverage += 1
                if new_coverage > best_gain:
                    best_gain = new_coverage
                    best_interval = interval
            if best_gain == 0 or best_interval is None:
                break
            selected.append(best_interval)
            start = best_interval["start"]
            end = best_interval["end"]
            for i in range(max(0, start), min(seq_len, end)):
                if not covered[i]:
                    covered[i] = True
                    current_coverage += 1
                    pbar.update(1)
            if current_coverage == seq_len:
                break

    return selected, current_coverage

def main():
    with open("results/valid_pairs.json", encoding='utf-8') as f:
        valid_pairs = json.load(f)
    with open("results/input_data.json", encoding='utf-8') as f:
        data = json.load(f)

    seq_len = len(data["sequence"])

    if not valid_pairs:
        print("WARNING: No valid amplicon pairs found.")
        result = {"selected_pairs": [], "coverage": 0, "total_length": seq_len, "percentage": 0, "num_pairs": 0, "time_seconds": 0}
        with open("results/best_coverage.json", "w", encoding='utf-8') as out:
            json.dump(result, out, indent=2)
        return

    print("\nStarting greedy optimization for overlapping amplicons...")
    start_time = time.time()
    selected, coverage = greedy_pair_coverage(valid_pairs, seq_len)
    elapsed = time.time() - start_time

    print(f"\nOptimization completed in {elapsed:.2f} seconds.")
    print(f"Covered {coverage} out of {seq_len} bases ({coverage/seq_len*100:.1f}%)")
    print(f"Number of amplicons used: {len(selected)}")

    result = {
        "selected_pairs": selected,
        "coverage": coverage,
        "total_length": seq_len,
        "percentage": coverage/seq_len*100,
        "num_pairs": len(selected),
        "time_seconds": elapsed
    }
    with open("results/best_coverage.json", "w", encoding='utf-8') as out:
        json.dump(result, out, indent=2)

    if selected:
        print("\nSelected amplicons (in order of selection):")
        for i, pair in enumerate(selected, 1):
            print(f"  {i:2d}. {pair['pair']} -> covers [{pair['start']}-{pair['end']}] (length: {pair['product_len']} bp, Tm diff: {pair['tm_diff']:.1f} deg C)")

if __name__ == "__main__":
    main()
