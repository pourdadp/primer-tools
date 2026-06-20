
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches

with open("results/binding_results.json", encoding='utf-8') as f:
    binding = json.load(f)
with open("results/input_data.json", encoding='utf-8') as f:
    data = json.load(f)

seq = data["sequence"]
fig, ax = plt.subplots(figsize=(15, 3))
ax.set_xlim(0, len(seq))
ax.set_ylim(0, 1)
ax.plot([0, len(seq)], [0.5, 0.5], 'k-', lw=2)

colors = plt.cm.tab10.colors
idx = 0
for name, info in binding.items():
    if info.get("status") != "ACCEPTED":
        continue
    color = colors[idx % len(colors)]
    idx += 1
    direction = info.get("direction", "unknown")
    for pos in info["positions"]:
        start, end, mism, gaps, _, _ = pos
        rect = patches.Rectangle((start, 0.3), end-start, 0.4,
                                 linewidth=1, edgecolor=color, facecolor=color, alpha=0.5)
        ax.add_patch(rect)
        label = f"{name}\n({direction},{mism}mis"
        if gaps:
            label += f",{gaps}gaps"
        label += ")"
        ax.text((start+end)/2, 0.8, label, ha='center', fontsize=8, color=color)

ax.set_yticks([])
ax.set_xlabel("Position in target sequence")
plt.title("Primer binding sites (auto-detected directions)")
plt.tight_layout()
plt.savefig("results/primer_binding.png")
