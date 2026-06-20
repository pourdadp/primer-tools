
import json

def melting_temp(seq):
    return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))

def check_dimer(seq1, seq2):
    if len(seq1) >= 3 and len(seq2) >= 3:
        return seq1[-3:] == seq2[-3:][::-1]
    return False

with open("results/binding_results.json", encoding='utf-8') as f:
    binding = json.load(f)

primer_tm = {}
for name, info in binding.items():
    if info["status"] == "ACCEPTED":
        primer_tm[name] = melting_temp(info["sequence"])

with open("results/primer_tm.json", "w", encoding='utf-8') as f:
    json.dump(primer_tm, f)

primers = list(binding.keys())
pair_info = {}

for i in range(len(primers)):
    for j in range(i+1, len(primers)):
        p1 = primers[i]
        p2 = primers[j]
        seq1 = binding[p1]["sequence"]
        seq2 = binding[p2]["sequence"]
        tm1 = melting_temp(seq1)
        tm2 = melting_temp(seq2)
        delta_tm = abs(tm1 - tm2)
        dimer = check_dimer(seq1, seq2)
        pair_info[f"{p1}-{p2}"] = {
            "tm1": tm1,
            "tm2": tm2,
            "delta_tm": delta_tm,
            "dimer_possible": dimer
        }

with open("results/pair_info.json", "w", encoding='utf-8') as out:
    json.dump(pair_info, out, indent=2)
