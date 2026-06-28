# quickngs/assembly.py
"""
DNA Assembly Module for QuickNGS
Algorithms: Greedy OLC, Semi‑Global Alignment, De Bruijn Graph
Supports AB1, FASTA, seq formats via Biopython
Handles multiple files, orientation control, quality trimming,
intelligent clustering, optional cluster merging, and guided assembly.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

from Bio import SeqIO
from io import StringIO
import os
import subprocess

# ---------- Tool availability check ----------
def is_tool_available(tool_name):
    try:
        subprocess.run([tool_name], capture_output=True, check=False)
        return True
    except FileNotFoundError:
        return False

# ---------- Quality Trimming ----------
def trim_by_quality(seq, qualities, threshold=20, window=5):
    if not qualities or len(qualities) != len(seq):
        return seq, qualities
    start = 0
    for i in range(0, len(qualities) - window + 1):
        avg = sum(qualities[i:i+window]) / window
        if avg >= threshold:
            start = i
            break
    end = len(qualities)
    for i in range(len(qualities) - window, -1, -1):
        avg = sum(qualities[i:i+window]) / window
        if avg >= threshold:
            end = i + window
            break
    return seq[start:end], qualities[start:end]

# ---------- File parsing ----------
def parse_uploaded_file(filepath, trim_low_quality=False, quality_threshold=20, window_size=5):
    filename = os.path.basename(filepath).lower()
    sequences = []
    trim_info = None
    try:
        if filename.endswith('.ab1'):
            record = SeqIO.read(filepath, 'abi')
            seq = str(record.seq)
            quals = record.letter_annotations.get('phred_quality', None)
            original_len = len(seq)
            if trim_low_quality and quals:
                seq, quals = trim_by_quality(seq, quals, quality_threshold, window_size)
                trim_info = {'original': original_len, 'trimmed': len(seq), 'removed': original_len - len(seq)}
            sequences.append(seq)
        elif filename.endswith(('.fasta', '.fa', '.fna')):
            for record in SeqIO.parse(filepath, 'fasta'):
                sequences.append(str(record.seq))
        else:
            with open(filepath, 'r') as f:
                lines = f.read().strip().split('\n')
                if lines and lines[0].startswith('>'):
                    for record in SeqIO.parse(StringIO('\n'.join(lines)), 'fasta'):
                        sequences.append(str(record.seq))
                else:
                    sequences = [line.strip() for line in lines if line.strip()]
    except Exception as e:
        raise ValueError(f"Could not read {filename}: {str(e)}")
    return sequences, trim_info

def parse_uploaded_files(filepaths, trim_low_quality=False, quality_threshold=20, window_size=5):
    all_reads = []
    all_trim_info = []
    for fp in filepaths:
        reads, trim_info = parse_uploaded_file(fp, trim_low_quality, quality_threshold, window_size)
        all_reads.extend(reads)
        if trim_info:
            all_trim_info.append(trim_info)
    return all_reads, all_trim_info

# ---------- Utilities ----------
def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(b, b) for b in reversed(seq))

# ---------- Overlap detection ----------
def fuzzy_overlap_ok(suffix, prefix, max_mismatch):
    run = 0
    for i in range(len(suffix)):
        if suffix[i] != prefix[i]:
            run += 1
            if run > max_mismatch:
                return False
        else:
            run = 0
    return True

def find_overlap(a, b, min_len, fuzzy_threshold=50, max_mismatch=3):
    best = 0
    max_len = min(len(a), len(b))
    for i in range(min_len, max_len + 1):
        a_suffix = a[-i:]
        b_prefix = b[:i]
        if i >= fuzzy_threshold:
            if fuzzy_overlap_ok(a_suffix, b_prefix, max_mismatch):
                best = i
        else:
            if a_suffix == b_prefix:
                best = i
    return best

# ---------- Greedy OLC ----------
def greedy_assemble(reads, min_overlap=3, fuzzy_threshold=50, max_mismatch=3, orientation='auto'):
    if not reads:
        return {'contig': '', 'steps': [], 'placements': [], 'remaining': [], 'total': 0}
    sorted_reads = sorted(reads, key=lambda x: -len(x))
    best_contig, best_steps, best_placements, best_remaining = '', [], [], sorted_reads[:]

    for seed_idx in range(len(sorted_reads)):
        temp_reads = sorted_reads[:]
        seed = temp_reads.pop(seed_idx)
        placements = [{'seq': seed, 'start': 0, 'end': len(seed) - 1, 'mismatches': set()}]
        steps = [{'seq': seed, 'type': 'seed', 'side': None, 'overlap': 0}]
        extended = True
        while temp_reads and extended:
            contig = resolve_and_build(placements)
            best_action = None
            for i, read in enumerate(temp_reads):
                if contig.find(read) != -1 or (orientation != 'forward' and contig.find(reverse_complement(read)) != -1):
                    best_action = {'type': 'contained', 'idx': i, 'read': read}
                    break
            if best_action:
                temp_reads.pop(best_action['idx'])
                steps.append({'seq': best_action['read'], 'type': 'contained'})
                continue
            best_overlap, best_idx, best_ori, best_side, best_seq = 0, -1, 'f', 'right', None
            for i, read in enumerate(temp_reads):
                ov = find_overlap(contig, read, min_overlap, fuzzy_threshold, max_mismatch)
                if ov > best_overlap:
                    best_overlap, best_idx, best_ori, best_side, best_seq = ov, i, 'f', 'right', read
                ov = find_overlap(read, contig, min_overlap, fuzzy_threshold, max_mismatch)
                if ov > best_overlap:
                    best_overlap, best_idx, best_ori, best_side, best_seq = ov, i, 'f', 'left', read
                if orientation != 'forward':
                    rc = reverse_complement(read)
                    ov = find_overlap(contig, rc, min_overlap, fuzzy_threshold, max_mismatch)
                    if ov > best_overlap:
                        best_overlap, best_idx, best_ori, best_side, best_seq = ov, i, 'r', 'right', rc
                    ov = find_overlap(rc, contig, min_overlap, fuzzy_threshold, max_mismatch)
                    if ov > best_overlap:
                        best_overlap, best_idx, best_ori, best_side, best_seq = ov, i, 'r', 'left', rc
            if best_idx == -1:
                extended = False
                break
            chosen_read = temp_reads.pop(best_idx)
            oriented_seq = best_seq if best_ori == 'f' else reverse_complement(chosen_read)
            steps.append({'seq': chosen_read, 'type': 'add', 'orientation': best_ori, 'side': best_side, 'overlap': best_overlap})
            mismatch_set = set()
            if best_overlap >= fuzzy_threshold:
                if best_side == 'right':
                    overlap_start = len(contig) - best_overlap
                    for j in range(best_overlap):
                        if contig[overlap_start + j] != oriented_seq[j]:
                            mismatch_set.add(overlap_start + j)
                else:
                    for j in range(best_overlap):
                        if oriented_seq[len(oriented_seq) - best_overlap + j] != contig[j]:
                            mismatch_set.add(j)
            if best_side == 'right':
                placements.append({'seq': chosen_read, 'start': len(contig) - best_overlap,
                                   'end': len(contig) - best_overlap + len(oriented_seq) - 1,
                                   'mismatches': mismatch_set})
            else:
                shift = len(oriented_seq) - best_overlap
                for p in placements:
                    p['start'] += shift
                    p['end'] += shift
                placements.append({'seq': chosen_read, 'start': 0, 'end': len(oriented_seq) - 1, 'mismatches': set()})
        contig = resolve_and_build(placements)
        if len(contig) > len(best_contig):
            best_contig, best_steps, best_placements, best_remaining = contig, steps, placements, temp_reads[:]
            break
    return {'contig': best_contig, 'steps': best_steps, 'placements': best_placements,
            'remaining': best_remaining, 'total': len(sorted_reads)}

def resolve_and_build(placements):
    if not placements:
        return ''
    max_pos = max(p['end'] for p in placements) + 1
    candidates = [[] for _ in range(max_pos)]
    for p in placements:
        seq = p['seq']
        for i in range(p['start'], p['end'] + 1):
            idx = i - p['start']
            candidates[i].append((seq[idx], min(idx, len(seq) - 1 - idx)))
    contig = []
    for c in candidates:
        best_base = max(c, key=lambda x: x[1])[0] if c else 'N'
        contig.append(best_base)
    return ''.join(contig)

# ---------- Semi‑Global Alignment ----------
MATCH, MISMATCH, GAP = 1, -1, -2

def semi_global_align(seq1, seq2):
    rows, cols = len(seq1) + 1, len(seq2) + 1
    score = [[0] * cols for _ in range(rows)]
    max_score, max_i, max_j = 0, 0, 0
    for i in range(1, rows):
        for j in range(1, cols):
            match = score[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            delete = score[i-1][j] + GAP
            insert = score[i][j-1] + GAP
            best = max(match, delete, insert, 0)
            score[i][j] = best
            if best > max_score:
                max_score, max_i, max_j = best, i, j
    i, j = max_i, max_j
    matches = contig_start = read_start = 0
    while i > 0 and j > 0 and score[i][j] > 0:
        if seq1[i-1] == seq2[j-1]:
            matches += 1
        if score[i][j] == score[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH):
            i -= 1; j -= 1
        elif score[i][j] == score[i-1][j] + GAP:
            i -= 1
        elif score[i][j] == score[i][j-1] + GAP:
            j -= 1
        else:
            break
        contig_start, read_start = i, j
    return {'score': max_score, 'matches': matches, 'contig_start': contig_start, 'contig_end': max_i,
            'read_start': read_start, 'read_end': max_j}

def merge_read_to_contig(contig, read, min_overlap, orientation='auto'):
    aln = semi_global_align(contig, read)
    best_aln = aln
    best_orientation = 'forward'
    if orientation != 'forward':
        rc = reverse_complement(read)
        aln_rc = semi_global_align(contig, rc)
        if aln_rc['matches'] > aln['matches']:
            best_aln = aln_rc
            best_orientation = 'reverse'
    if best_aln['matches'] >= min_overlap:
        read_seq = read if best_orientation == 'forward' else reverse_complement(read)
        left_overhang = read_seq[:best_aln['read_start']]
        right_overhang = read_seq[best_aln['read_end']:]
        left_contig = contig[:best_aln['contig_start']]
        right_contig = contig[best_aln['contig_end']:]
        return left_overhang + left_contig + contig[best_aln['contig_start']:best_aln['contig_end']] + right_contig + right_overhang
    return None

def assemble_with_alignment(reads, min_overlap=3, orientation='auto'):
    if not reads:
        return ''
    sorted_reads = sorted(reads, key=lambda x: -len(x))
    contig = sorted_reads[0]
    used = {0}
    changed = True
    while changed:
        changed = False
        for i, read in enumerate(sorted_reads):
            if i in used:
                continue
            new_contig = merge_read_to_contig(contig, read, min_overlap, orientation)
            if new_contig:
                contig = new_contig
                used.add(i)
                changed = True
                break
    return contig

# ---------- De Bruijn Graph ----------
def debruijn_assemble(reads, k=3):
    if not reads:
        return ''
    graph = {}
    for read in
