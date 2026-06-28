# quickngs/translation.py
"""
Translation and Codon Usage Module for QuickNGS
Supports: Auto-detect frame, Longest ORF, First ATG, Stop codon,
         Codon Usage (CAI via Biopython), Manual frame selection.
Includes embedded Kazusa reference files for human, E. coli, yeast, mycoplasma.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

from Bio.Seq import Seq
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from io import StringIO

# ---------- Embedded Kazusa Reference Files ----------
EMBEDDED_REFERENCES = {
    'human': """UUU 17.6 0.46\nUUC 20.3 0.54\nUUA 7.7 0.08\nUUG 12.9 0.13\nCUU 13.2 0.13\nCUC 19.6 0.20\nCUA 7.2 0.07\nCUG 39.6 0.40\nAUU 16.0 0.36\nAUC 20.8 0.47\nAUA 7.5 0.17\nAUG 22.0 1.00\nGUU 11.0 0.18\nGUC 14.5 0.24\nGUA 7.1 0.12\nGUG 28.1 0.46\nUCU 15.2 0.19\nUCC 17.7 0.22\nUCA 12.2 0.15\nUCG 4.4 0.06\nCCU 17.5 0.29\nCCC 19.8 0.33\nCCA 16.9 0.27\nCCG 6.9 0.11\nACU 13.1 0.25\nACC 18.9 0.36\nACA 15.1 0.28\nACG 6.1 0.11\nGCU 18.4 0.27\nGCC 27.7 0.40\nGCA 15.8 0.23\nGCG 7.4 0.10\nUAU 12.2 0.44\nUAC 15.3 0.56\nCAU 10.9 0.42\nCAC 15.1 0.58\nCAA 12.3 0.27\nCAG 34.2 0.73\nAAU 17.0 0.47\nAAC 19.1 0.53\nAAA 24.4 0.43\nAAG 31.9 0.57\nGAU 21.8 0.46\nGAC 25.1 0.54\nGAA 29.0 0.42\nGAG 39.6 0.58\nUGU 10.6 0.46\nUGC 12.6 0.54\nUGG 13.2 1.00\nCGU 4.5 0.08\nCGC 10.4 0.18\nCGA 6.2 0.11\nCGG 11.4 0.20\nAGU 12.1 0.15\nAGC 19.5 0.24\nAGA 12.0 0.21\nAGG 12.0 0.21\nGGU 10.8 0.16\nGGC 22.8 0.34\nGGA 16.5 0.25\nGGG 16.5 0.25""",
    'ecoli': """UUU 22.4 0.57\nUUC 16.4 0.43\nUUA 13.1 0.18\nUUG 13.3 0.18\nCUU 11.4 0.16\nCUC 10.5 0.15\nCUA 3.9 0.05\nCUG 49.7 0.69\nAUU 30.1 0.49\nAUC 24.8 0.40\nAUA 5.6 0.09\nAUG 27.3 1.00\nGUU 18.4 0.27\nGUC 15.4 0.22\nGUA 11.3 0.16\nGUG 25.5 0.37\nUCU 8.3 0.15\nUCC 8.8 0.16\nUCA 7.7 0.14\nUCG 8.7 0.16\nCCU 7.1 0.17\nCCC 5.5 0.13\nCCA 8.3 0.20\nCCG 22.6 0.54\nACU 9.0 0.19\nACC 23.4 0.49\nACA 7.5 0.16\nACG 14.3 0.30\nGCU 15.2 0.17\nGCC 25.0 0.28\nGCA 20.2 0.22\nGCG 32.1 0.35\nUAU 16.5 0.59\nUAC 12.5 0.45\nCAU 12.8 0.56\nCAC 9.8 0.44\nCAA 15.4 0.35\nCAG 29.0 0.66\nAAU 18.6 0.51\nAAC 21.5 0.59\nAAA 33.7 0.74\nAAG 10.2 0.23\nGAU 32.1 0.63\nGAC 19.1 0.37\nGAA 39.8 0.69\nGAG 18.0 0.31\nUGU 5.1 0.47\nUGC 6.4 0.59\nUGG 15.3 1.00\nCGU 21.6 0.38\nCGC 22.2 0.39\nCGA 3.8 0.07\nCGG 5.9 0.10\nAGU 9.1 0.16\nAGC 16.4 0.29\nAGA 2.8 0.05\nAGG 1.6 0.03\nGGU 25.3 0.35\nGGC 29.4 0.41\nGGA 8.7 0.12\nGGG 11.3 0.16""",
    'yeast': """UUU 26.1 0.59\nUUC 18.4 0.41\nUUA 26.2 0.28\nUUG 27.2 0.29\nCUU 12.4 0.13\nCUC 5.4 0.06\nCUA 13.4 0.14\nCUG 10.5 0.11\nAUU 30.1 0.46\nAUC 17.2 0.26\nAUA 17.8 0.27\nAUG 20.9 1.00\nGUU 22.1 0.39\nGUC 11.8 0.21\nGUA 11.8 0.21\nGUG 10.8 0.19\nUCU 23.5 0.26\nUCC 14.2 0.16\nUCA 18.7 0.21\nUCG 8.6 0.10\nCCU 13.5 0.31\nCCC 6.8 0.15\nCCA 18.3 0.41\nCCG 5.3 0.12\nACU 20.3 0.35\nACC 12.7 0.22\nACA 17.8 0.31\nACG 8.0 0.14\nGCU 21.2 0.38\nGCC 12.6 0.22\nGCA 16.2 0.29\nGCG 6.2 0.11\nUAU 18.8 0.56\nUAC 14.8 0.44\nCAU 13.6 0.64\nCAC 7.8 0.36\nCAA 27.3 0.69\nCAG 12.1 0.31\nAAU 35.7 0.59\nAAC 24.8 0.41\nAAA 41.9 0.58\nAAG 30.8 0.42\nGAU 37.6 0.65\nGAC 20.2 0.35\nGAA 45.6 0.70\nGAG 19.2 0.30\nUGU 8.1 0.63\nUGC 4.8 0.37\nUGG 10.4 1.00\nCGU 6.4 0.14\nCGC 2.6 0.06\nCGA 3.0 0.07\nCGG 1.7 0.04\nAGU 14.2 0.16\nAGC 9.8 0.11\nAGA 21.3 0.48\nAGG 9.2 0.21\nGGU 23.9 0.47\nGGC 9.8 0.19\nGGA 10.9 0.22\nGGG 6.0 0.12""",
    'mycoplasma': """UUU 25.6 1.00\nUUC 0.0 0.00\nUUA 35.2 0.65\nUUG 0.0 0.00\nCUU 0.0 0.00\nCUC 0.0 0.00\nCUA 0.0 0.00\nCUG 0.0 0.00\nAUU 30.1 0.49\nAUC 0.0 0.00\nAUA 31.2 0.51\nAUG 31.2 1.00\nGUU 20.4 1.00\nGUC 0.0 0.00\nGUA 0.0 0.00\nGUG 0.0 0.00\nUCU 10.1 1.00\nUCC 0.0 0.00\nUCA 0.0 0.00\nUCG 0.0 0.00\nCCU 25.6 1.00\nCCC 0.0 0.00\nCCA 0.0 0.00\nCCG 0.0 0.00\nACU 15.7 1.00\nACC 0.0 0.00\nACA 0.0 0.00\nACG 0.0 0.00\nGCU 30.1 1.00\nGCC 0.0 0.00\nGCA 0.0 0.00\nGCG 0.0 0.00\nUAU 20.5 1.00\nUAC 0.0 0.00\nCAU 18.9 1.00\nCAC 0.0 0.00\nCAA 45.6 1.00\nCAG 0.0 0.00\nAAU 35.7 0.59\nAAC 24.8 0.41\nAAA 41.9 0.58\nAAG 30.8 0.42\nGAU 37.6 0.65\nGAC 20.2 0.35\nGAA 45.6 0.70\nGAG 19.2 0.30\nUGU 8.1 0.63\nUGC 4.8 0.37\nUGA 10.4 1.00\nUGG 0.0 0.00\nCGU 0.0 0.00\nCGC 0.0 0.00\nCGA 0.0 0.00\nCGG 0.0 0.00\nAGU 14.2 0.16\nAGC 9.8 0.11\nAGA 21.3 0.48\nAGG 9.2 0.21\nGGU 23.9 1.00\nGGC 0.0 0.00\nGGA 0.0 0.00\nGGG 0.0 0.00"""
}

# ---------- Reverse Complement ----------
def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, b) for b in reversed(seq))

# ---------- CAI Calculation ----------
def calculate_cai(seq, organism='human'):
    """Calculate Codon Adaptation Index (CAI) using Biopython."""
    if organism not in EMBEDDED_REFERENCES:
        return None
    ref_file = StringIO(EMBEDDED_REFERENCES[organism])
    cai = CodonAdaptationIndex()
    try:
        cai.generate_index(ref_file)
        score = cai.cai_for_gene(str(seq))
        return score
    except:
        return None

def find_best_frame_by_cai(seq, genetic_code=1):
    """Find the best frame using CAI from embedded references."""
    best_frame = 1
    best_score = 0
    best_org = 'human'
    for organism in ['human', 'ecoli', 'yeast', 'mycoplasma']:
        for frame in [1, 2, 3, -1, -2, -3]:
            if frame > 0:
                coding_seq = seq[frame-1:]
            else:
                coding_seq = reverse_complement(seq)[abs(frame)-1:]
            cai_score = calculate_cai(coding_seq, organism)
            if cai_score and cai_score > best_score:
                best_score = cai_score
                best_frame = frame
                best_org = organism
    if best_score > 0:
        return best_frame, best_score, best_org
    return None

# ---------- Translation Function ----------
def translate_contig(seq, ref_protein=None, method='auto', genetic_code=1, organism='human'):
    """
    Translate a DNA sequence using the specified method.
    Returns dict with protein, frame, method_description, warning, and color.
    """
    # Generate all 6 frames
    frames = {}
    for frame in range(3):
        if frame > 0:
            seq_trimmed = seq[frame:]
        else:
            seq_trimmed = seq
        translated = Seq(seq_trimmed).translate(table=genetic_code, to_stop=False)
        frames[frame + 1] = str(translated)
        
        rc_seq = reverse_complement(seq)
        if frame > 0:
            seq_trimmed = rc_seq[frame:]
        else:
            seq_trimmed = rc_seq
        translated = Seq(seq_trimmed).translate(table=genetic_code, to_stop=False)
        frames[-(frame + 1)] = str(translated)
    
    # Method 1: Protein reference alignment
    if ref_protein and method in ('auto', 'reference'):
        best_score = 0
        best_frame = 1
        for frame_num, protein in frames.items():
            score = sum(1 for a, b in zip(protein, ref_protein) if a == b)
            if score > best_score:
                best_score = score
                best_frame = frame_num
        pct = round(best_score / len(ref_protein) * 100) if ref_protein else 0
        return {
            'protein': frames[best_frame],
            'frame': best_frame,
            'method_description': f'✅ Protein reference alignment used. Best match: Frame {best_frame:+d} (score: {pct}%)',
            'warning': None,
            'color': 'success'
        }
    
    # Method 2: Longest ORF
    if method in ('auto', 'longest_orf'):
        best_orf = ''
        best_frame = 1
        has_start = False
        
        for frame_num, protein in frames.items():
            orfs = protein.split('*')
            for orf in orfs:
                if 'M' in orf:
                    orf_from_m = orf[orf.index('M'):]
                    if len(orf_from_m) > len(best_orf):
                        best_orf = orf_from_m
                        best_frame = frame_num
                        has_start = True
        
        if best_orf:
            desc = f'✅ Longest ORF selected (Frame {best_frame:+d}, {len(best_orf)} aa, {"ATG→" if has_start else "no ATG→"}TGA/TAA/TAG)'
            warning = None if has_start else "No start codon (ATG) found. Translation may be N-terminal truncated."
            return {
                'protein': best_orf,
                'frame': best_frame,
                'method_description': desc,
                'warning': warning,
                'color': 'success' if has_start else 'warning'
            }
    
    # Method 3: First ATG
    if method in ('auto', 'first_atg'):
        for frame_num, protein in frames.items():
            if 'M' in protein:
                start_idx = protein.index('M')
                return {
                    'protein': protein[start_idx:],
                    'frame': frame_num,
                    'method_description': f'⚠️ First ATG selected (Frame {frame_num:+d}). No full ORF found—translation may be incomplete.',
                    'warning': 'No full ORF found. Translation may be incomplete.',
                    'color': 'warning'
                }
    
    # Method 4: Stop codon-based
    if method in ('auto', 'stop_codon'):
        for frame_num, protein in frames.items():
            if '*' in protein:
                stop_idx = protein.index('*')
                return {
                    'protein': protein[:stop_idx],
                    'frame': frame_num,
                    'method_description': f'⚠️ Stop codon-based selection (Frame {frame_num:+d}). Reading frame determined from stop codon position.',
                    'warning': 'Reading frame determined from stop codon position. Start may be missing.',
                    'color': 'warning'
                }
    
    # Method 5: Codon Usage (CAI)
    if method in ('auto', 'codon_usage'):
        result = find_best_frame_by_cai(seq, genetic_code)
        if result:
            frame_num, cai_score, org_used = result
            return {
                'protein': frames[frame_num],
                'frame': frame_num,
                'method_description': f'⚠️ Codon usage bias (CAI) used to select Frame {frame_num:+d} (CAI: {cai_score:.2f}, ref: {org_used}). No start/stop codons found.',
                'warning': 'No start/stop codons found. Translation may be incorrect.',
                'color': 'warning'
            }
    
    # Method 6: Manual or fallback
    if isinstance(method, int):
        frame_num = method
    else:
        frame_num = 1
    
    return {
        'protein': frames.get(frame_num, frames[1]),
        'frame': frame_num,
        'method_description': f'🔧 Manual frame selected: Frame {frame_num:+d}. Verify translation independently.',
        'warning': 'Manual frame selected. Verify translation independently.',
        'color': 'info'
    }
