# quickngs/msa.py
"""
Multiple Sequence Alignment Module for QuickNGS
Uses MUSCLE (local) or Biopython's built‑in ClustalW wrapper.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

import subprocess
import os
import tempfile
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

# ---------- MUSCLE Alignment (Local, requires MUSCLE installed) ----------
def run_muscle(sequences, labels=None):
    """
    Run MUSCLE multiple sequence alignment.
    sequences: list of DNA/protein strings
    labels: optional list of names (default: seq_1, seq_2, ...)
    Returns: dict with alignment (FASTA), consensus, and HTML table.
    """
    if not labels:
        labels = [f"seq_{i+1}" for i in range(len(sequences))]
    
    # Write sequences to temp FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        for label, seq in zip(labels, sequences):
            f.write(f">{label}\n{seq}\n")
        input_file = f.name
    
    output_file = input_file + '.aln'
    
    try:
        # Run MUSCLE
        subprocess.run(['muscle', '-align', input_file, '-output', output_file],
                      check=True, capture_output=True, timeout=120)
        
        # Read alignment
        alignment = AlignIO.read(output_file, 'fasta')
        
        # Generate consensus
        consensus = ''
        for col in range(alignment.get_alignment_length()):
            bases = alignment[:, col]
            # Most common base (or gap)
            base_counts = {}
            for b in bases:
                if b != '-':
                    base_counts[b] = base_counts.get(b, 0) + 1
            if base_counts:
                consensus += max(base_counts, key=base_counts.get)
            else:
                consensus += '-'
        
        # Generate HTML with coloring
        html = generate_msa_html(alignment, labels)
        
        # Cleanup
        os.unlink(input_file)
        os.unlink(output_file)
        
        return {
            'alignment': str(alignment),
            'consensus': consensus,
            'html': html,
            'length': alignment.get_alignment_length(),
            'seq_count': len(sequences)
        }
    except subprocess.TimeoutExpired:
        return {'error': 'MUSCLE alignment timed out (>2 minutes).'}
    except FileNotFoundError:
        return {'error': 'MUSCLE is not installed. Please install: sudo apt install muscle'}
    except Exception as e:
        return {'error': str(e)}
    finally:
        if os.path.exists(input_file):
            os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

# ---------- Biopython Pairwise Fallback (No MUSCLE needed) ----------
def run_pairwise_alignment(seq1, seq2, label1='seq_1', label2='seq_2'):
    """
    Simple pairwise alignment using Biopython's built‑in pairwise2.
    Always available—no external tools needed.
    """
    from Bio import pairwise2
    
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
    if not alignments:
        return {'error': 'Could not align sequences.'}
    
    best = alignments[0]
    aligned1, aligned2, score, begin, end = best
    
    seq1_record = SeqRecord(Seq(aligned1), id=label1)
    seq2_record = SeqRecord(Seq(aligned2), id=label2)
    alignment = MultipleSeqAlignment([seq1_record, seq2_record])
    
    # Generate HTML
    html = generate_msa_html(alignment, [label1, label2])
    
    return {
        'alignment': str(alignment),
        'html': html,
        'length': len(aligned1),
        'seq_count': 2
    }

# ---------- HTML Generation with Coloring ----------
def generate_msa_html(alignment, labels):
    """Generate colored HTML table for MSA visualization."""
    # Color scheme for amino acids / nucleotides
    colors = {
        'A': '#ccffcc', 'G': '#ccffcc',  # Green for small nonpolar
        'T': '#ffe6e6', 'C': '#ffe6e6',  # Red for pyrimidines
        'U': '#ffe6e6',
        'D': '#ffcccc', 'E': '#ffcccc',  # Red for acidic
        'R': '#ccccff', 'K': '#ccccff', 'H': '#ccccff',  # Blue for basic
        'F': '#ffffcc', 'Y': '#ffffcc', 'W': '#ffffcc',  # Yellow for aromatic
        'I': '#ccffcc', 'L': '#ccffcc', 'V': '#ccffcc', 'M': '#ccffcc', 'P': '#ccffcc',
        'N': '#ffe6cc', 'Q': '#ffe6cc', 'S': '#ffe6cc', 'T': '#ffe6cc',
        '-': '#f0f0f0'  # Gray for gaps
    }
    
    html = '<div class="table-responsive"><table class="table table-sm table-bordered msa-table" style="font-family: monospace; font-size: 12px;">'
    
    # Header row with position markers
    html += '<tr><th></th>'
    length = alignment.get_alignment_length()
    for i in range(1, length + 1):
        if i % 10 == 0:
            html += f'<th style="padding:0 2px;font-size:10px;text-align:center;">{i}</th>'
        else:
            html += '<th style="padding:0 2px;"></th>'
    html += '</tr>'
    
    # Sequence rows
    for i, record in enumerate(alignment):
        html += f'<tr><td style="font-weight:bold;white-space:nowrap;">{labels[i] if i < len(labels) else record.id}</td>'
        for base in str(record.seq):
            color = colors.get(base.upper(), '#ffffff')
            html += f'<td style="background:{color};padding:0 2px;text-align:center;">{base}</td>'
        html += '</tr>'
    
    html += '</table></div>'
    return html
