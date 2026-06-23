import json
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_binding_report(binding_results, seq_len, output_dir):
    """Generate a formatted HTML report with a printable table of binding results."""
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Binding Results Report</title>
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 20px; background: #fff; }}
            h2 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; font-size: 11pt; }}
            th, td {{ border: 1px solid #333; padding: 8px 10px; text-align: center; }}
            th {{ background-color: #2c3e50; color: white; font-weight: bold; }}
            tr:nth-child(even) {{ background-color: #f9f9f9; }}
            .accepted {{ background-color: #d4edda !important; }}
            .rejected {{ background-color: #f8d7da !important; }}
            .legend {{ margin: 15px 0; padding: 10px; background: #f8f9fa; border: 1px solid #ddd; border-radius: 5px; }}
            .legend-item {{ display: inline-block; margin-right: 20px; }}
            .legend-color {{ display: inline-block; width: 20px; height: 20px; margin-right: 5px; vertical-align: middle; }}
            .legend-color.green {{ background-color: #d4edda; border: 1px solid #155724; }}
            .legend-color.red {{ background-color: #f8d7da; border: 1px solid #721c24; }}
            .sequence {{ font-family: monospace; font-size: 9pt; word-break: break-all; }}
            .description {{ font-size: 9pt; color: #555; text-align: left; padding: 4px 8px; background: #f1f1f1; border-radius: 3px; }}
            .field-desc {{ display: inline-block; margin: 2px 5px; font-size: 8pt; color: #666; }}
            @media print {{
                body {{ margin: 10px; font-size: 9pt; }}
                th, td {{ padding: 4px 6px; }}
                .no-print {{ display: none !important; }}
                th {{ background-color: #2c3e50 !important; color: white !important; }}
                .accepted {{ background-color: #d4edda !important; }}
                .rejected {{ background-color: #f8d7da !important; }}
            }}
        </style>
    </head>
    <body>
        <button onclick="window.print()" class="no-print" style="padding:8px 16px; background:#2c3e50; color:white; border:none; border-radius:4px; cursor:pointer; float:right;">🖨️ Print</button>
        <h2>🧬 Primer Binding Results</h2>
        <p><strong>Sequence Length:</strong> {} bp</p>
        <p><strong>Total Primers:</strong> {} &nbsp;|&nbsp;
           <strong>Accepted:</strong> <span style="color:#155724;">{}</span> &nbsp;|&nbsp;
           <strong>Rejected:</strong> <span style="color:#721c24;">{}</span></p>

        <div class="legend">
            <span class="legend-item"><span class="legend-color green"></span> Accepted</span>
            <span class="legend-item"><span class="legend-color red"></span> Rejected</span>
        </div>

        <table>
            <thead>
                <tr>
                    <th>Primer</th>
                    <th>Status</th>
                    <th>Direction</th>
                    <th>Position (start, end)</th>
                    <th>Mismatches</th>
                    <th>Gaps</th>
                    <th>Score</th>
                    <th>Sequence</th>
                </tr>
            </thead>
            <tbody>
    """

    for name, info in binding_results.items():
        status = info.get("status", "UNKNOWN")
        cls = "accepted" if status == "ACCEPTED" else "rejected"
        seq = info.get("sequence", "-")
        direction = info.get("direction", "-")
        reason = info.get("reason", "")

        if status == "ACCEPTED":
            positions = info.get("positions", [])
            if positions:
                for pos in positions:
                    start, end, mism, gaps, dir_val, score = pos
                    html_content += f"""
                <tr class="{cls}">
                    <td><strong>{name}</strong></td>
                    <td>{status}</td>
                    <td>{dir_val}</td>
                    <td>{start} → {end}</td>
                    <td>{mism}</td>
                    <td>{gaps}</td>
                    <td>{score}</td>
                    <td class="sequence">{seq}</td>
                </tr>
                    """
            else:
                html_content += f"""
                <tr class="{cls}">
                    <td><strong>{name}</strong></td>
                    <td>{status}</td>
                    <td>{direction}</td>
                    <td colspan="4">No binding positions</td>
                    <td class="sequence">{seq}</td>
                </tr>
                """
        else:
            html_content += f"""
                <tr class="{cls}">
                    <td><strong>{name}</strong></td>
                    <td>{status}</td>
                    <td>-</td>
                    <td colspan="4" style="text-align:left; font-size:9pt; color:#721c24;">⚠️ {reason}</td>
                    <td class="sequence">{seq}</td>
                </tr>
            """

    html_content += """
            </tbody>
        </table>

        <div style="margin-top:30px; font-size:9pt; color:#555; border-top:1px solid #ddd; padding-top:10px;">
            <h4>📖 Legend: Position Fields</h4>
            <table style="width:auto; font-size:9pt; border-collapse:collapse;">
                <thead>
                    <tr><th>Field</th><th>Description</th></tr>
                </thead>
                <tbody>
                    <tr><td><strong>start</strong></td><td>Start position of binding (0‑based)</td></tr>
                    <tr><td><strong>end</strong></td><td>End position of binding (exclusive)</td></tr>
                    <tr><td><strong>Mismatches</strong></td><td>Number of mismatches (substitutions)</td></tr>
                    <tr><td><strong>Gaps</strong></td><td>Number of indels (insertions/deletions)</td></tr>
                    <tr><td><strong>Direction</strong></td><td>Forward or Reverse orientation</td></tr>
                    <tr><td><strong>Score</strong></td><td>Semi‑global alignment score (higher = better)</td></tr>
                </tbody>
            </table>
        </div>
    </body>
    </html>
    """

    total = len(binding_results)
    accepted = sum(1 for v in binding_results.values() if v.get("status") == "ACCEPTED")
    rejected = total - accepted

    # Format the HTML with the actual values (order: seq_len, total, accepted, rejected)
    html_content = html_content.format(seq_len, total, accepted, rejected)

    report_path = os.path.join(output_dir, "binding_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    return report_path

# ==================== Original visualization function ====================

def main():
    with open("results/binding_results.json", encoding='utf-8') as f:
        binding = json.load(f)
    with open("results/input_data.json", encoding='utf-8') as f:
        data = json.load(f)

    seq = data["sequence"]
    seq_len = len(seq)

    # Generate binding report (HTML table)
    report_path = generate_binding_report(binding, seq_len, "results")
    print(f"Binding report generated: {report_path}")

    # ===== Draw binding map (original PNG) =====
    fig, ax = plt.subplots(figsize=(15, 3))
    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, 1)
    ax.plot([0, seq_len], [0.5, 0.5], 'k-', lw=2)

    colors = plt.cm.tab10.colors
    idx = 0
    for name, info in binding.items():
        if info.get("status") != "ACCEPTED":
            continue
        color = colors[idx % len(colors)]
        idx += 1
        direction = info.get("direction", "unknown")
        for pos in info.get("positions", []):
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
    print("Binding map saved: results/primer_binding.png")

if __name__ == "__main__":
    main()
