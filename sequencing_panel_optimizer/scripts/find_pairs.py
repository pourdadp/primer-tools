
import json

with open("results/binding_results.json", encoding='utf-8') as f:
    binding = json.load(f)
with open("results/input_data.json", encoding='utf-8') as f:
    data = json.load(f)
with open("results/primer_tm.json", encoding='utf-8') as f:
    primer_tm = json.load(f)
with open("results/pair_info.json", encoding='utf-8') as f:
    pair_info = json.load(f)

max_len = data["max_product_length"]
min_len = data["min_product_length"]
max_tm_diff = data["max_tm_diff"]

valid_pairs = {}

forward_primers = [name for name, info in binding.items() if info.get("status") == "ACCEPTED" and info.get("direction") == "forward"]
reverse_primers = [name for name, info in binding.items() if info.get("status") == "ACCEPTED" and info.get("direction") == "reverse"]

for f_primer in forward_primers:
    for r_primer in reverse_primers:
        pair_key = f"{f_primer}-{r_primer}"
        if pair_info.get(pair_key, {}).get("dimer_possible", False):
            continue
        if abs(primer_tm[f_primer] - primer_tm[r_primer]) > max_tm_diff:
            continue

        f_info = binding[f_primer]
        r_info = binding[r_primer]
        for f_pos in f_info["positions"]:
            for r_pos in r_info["positions"]:
                start_f, end_f, _, _, _, _ = f_pos
                start_r, end_r, _, _, _, _ = r_pos
                product_len = start_r - end_f
                if min_len <= product_len <= max_len:
                    valid_pairs.setdefault(pair_key, []).append({
                        "f_start": start_f,
                        "f_end": end_f,
                        "r_start": start_r,
                        "r_end": end_r,
                        "product_length": product_len,
                        "tm_f": primer_tm[f_primer],
                        "tm_r": primer_tm[r_primer],
                        "tm_diff": abs(primer_tm[f_primer] - primer_tm[r_primer])
                    })

with open("results/valid_pairs.json", "w", encoding='utf-8') as out:
    json.dump(valid_pairs, out, indent=2)
