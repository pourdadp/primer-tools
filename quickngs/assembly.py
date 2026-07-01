# quickngs/assembly.py
"""
DNA Assembly Module for QuickNGS
Algorithms: Greedy OLC, Semi‑Global Alignment, De Bruijn Graph
Supports AB1, FASTA, seq formats via Biopython
Handles multiple files, orientation control, quality trimming,
intelligent clustering, optional cluster merging, guided assembly,
preserves original filenames, and reports trimming details.
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

# ---------- File parsing (now returns labels) ----------
def parse_uploaded_file(filepath, trim_low_quality=False, quality_threshold=20, window_size=5):
    filename = os.path.basename(filepath)
    sequences = []
    trim_info = None
    labels = []
    try:
        if filename.lower().endswith('.ab1'):
            record = SeqIO.read(filepath, 'abi')
            seq = str(record.seq)
            quals = record.letter_annotations.get('phred_quality', None)
            original_len = len(seq)
            if trim_low_quality and quals:
                seq, quals = trim_by_quality(seq, quals, quality_threshold, window_size)
                trim_info = {'original': original_len, 'trimmed': len(seq), 'removed': original_len - len(seq)}
            sequences.append(seq)
            labels.append(filename)
        elif filename.lower().endswith(('.fasta', '.fa', '.fna')):
            records = list(SeqIO.parse(filepath, 'fasta'))
            for i, record in enumerate(records):
                seq = str(record.seq)
                sequences.append(seq)
                if len(records) > 1:
                    labels.append(f"{filename} (seq {i+1})")
                else:
                    labels.append(filename)
        else:  # plain text / .seq / .txt
            with open(filepath, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]
                if lines and lines[0].startswith('>'):
                    records = list(SeqIO.parse(StringIO('\n'.join(lines)), 'fasta'))
                    for i, record in enumerate(records):
                        seq = str(record.seq)
                        sequences.append(seq)
                        if len(records) > 1:
                            labels.append(f"{filename} (seq {i+1})")
                        else:
                            labels.append(filename)
                else:
                    for i, line in enumerate(lines):
                        sequences.append(line)
                        if len(lines) > 1:
                            labels.append(f"{filename} (line {i+1})")
                        else:
                            labels.append(filename)
    except Exception as e:
        raise ValueError(f"Could not read {filename}: {str(e)}")
    return sequences, trim_info, labels

def parse_uploaded_files(filepaths, trim_low_quality=False, quality_threshold=20, window_size=5):
    all_reads = []
    all_trim_info = []
    all_labels = []
    trim_details = []
    for fp in filepaths:
        reads, trim_info, labels = parse_uploaded_file(fp, trim_low_quality, quality_threshold, window_size)
        all_reads.extend(reads)
        all_labels.extend(labels)
        if trim_info:
            trim_details.append({
                'filename': os.path.basename(fp),
                'original': trim_info['original'],
                'trimmed': trim_info['trimmed'],
                'removed': trim_info['removed']
            })
    return all_reads, all_trim_info, all_labels, trim_details

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

# ---------- De Bruijn Graph (Fixed: cycle detection, max length, auto k) ----------
def debruijn_assemble(reads, k=3):
    if not reads:
        return ''
    avg_read_len = sum(len(r) for r in reads) / len(reads)
    for attempt in range(3):
        graph = {}
        for read in reads:
            for i in range(len(read) - k + 1):
                kmer = read[i:i+k]
                prefix = kmer[:-1]
                suffix = kmer[1:]
                graph.setdefault(prefix, []).append(suffix)
        if not graph:
            return ''
        start = max(graph.keys(), key=lambda x: len(graph[x]))
        contig = start
        current = start
        visited_edges = set()
        max_contig_length = 100000
        while current in graph and graph[current]:
            if len(contig) > max_contig_length:
                break
            next_node = graph[current].pop(0)
            edge = (current, next_node)
            if edge in visited_edges:
                break
            visited_edges.add(edge)
            contig += next_node[-1]
            current = next_node
        if len(contig) >= avg_read_len * 1.5 or k >= 7:
            return contig
        k += 2
    return contig

# ---------- Intelligent Clustering ----------
def cluster_reads(reads, min_overlap=3, max_mismatch=3):
    n = len(reads)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            ov = find_overlap(reads[i], reads[j], min_overlap, 50, max_mismatch)
            if ov >= min_overlap:
                graph[i].add(j)
                graph[j].add(i)
            else:
                ov = find_overlap(reads[j], reads[i], min_overlap, 50, max_mismatch)
                if ov >= min_overlap:
                    graph[i].add(j)
                    graph[j].add(i)
    visited = set()
    clusters = []
    for i in range(n):
        if i not in visited:
            component = []
            stack = [i]
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    component.append(node)
                    stack.extend(graph[node] - visited)
            clusters.append(component)
    return clusters

# ---------- Cluster Merging (Scaffolding) ----------
def merge_clusters(contigs, min_overlap=10, max_mismatch=0):
    if len(contigs) <= 1:
        return contigs, {i: [i+1] for i in range(len(contigs))}
    n = len(contigs)
    merged = list(contigs)
    merge_map = {i: [i] for i in range(n)}
    changed = True
    while changed:
        changed = False
        for i in range(len(merged)):
            for j in range(i+1, len(merged)):
                if not merged[i] or not merged[j]:
                    continue
                ov = find_overlap(merged[i], merged[j], min_overlap, 50, max_mismatch)
                if ov >= min_overlap:
                    new_contig = merged[i] + merged[j][ov:]
                    merged[i] = new_contig
                    merge_map[i].extend(merge_map.pop(j, []))
                    merged[j] = ''
                    changed = True
                    break
                ov = find_overlap(merged[j], merged[i], min_overlap, 50, max_mismatch)
                if ov >= min_overlap:
                    new_contig = merged[j] + merged[i][ov:]
                    merged[j] = new_contig
                    merge_map[j].extend(merge_map.pop(i, []))
                    merged[i] = ''
                    changed = True
                    break
                rc_i = reverse_complement(merged[i])
                ov = find_overlap(rc_i, merged[j], min_overlap, 50, max_mismatch)
                if ov >= min_overlap:
                    new_contig = rc_i + merged[j][ov:]
                    merged[i] = new_contig
                    merge_map[i].extend(merge_map.pop(j, []))
                    merged[j] = ''
                    changed = True
                    break
                rc_j = reverse_complement(merged[j])
                ov = find_overlap(rc_j, merged[i], min_overlap, 50, max_mismatch)
                if ov >= min_overlap:
                    new_contig = rc_j + merged[i][ov:]
                    merged[j] = new_contig
                    merge_map[j].extend(merge_map.pop(i, []))
                    merged[i] = ''
                    changed = True
                    break
            if changed:
                break
        merged = [c for c in merged if c]
    final_merge_map = {}
    new_idx = 0
    for old_idx in range(len(contigs)):
        for m_idx, m_list in merge_map.items():
            if old_idx in m_list:
                if m_idx not in final_merge_map:
                    final_merge_map[m_idx] = []
                final_merge_map[m_idx].append(old_idx + 1)
                break
    return merged, final_merge_map

# ---------- Guided Assembly (with bcftools pipeline) ----------
def guided_assemble_fastq(r1_fastq, ref_fasta, results_folder):
    for tool in ['bwa', 'samtools', 'bcftools']:
        if not is_tool_available(tool):
            raise RuntimeError(f"{tool} is not installed.")

    if not os.path.exists(ref_fasta):
        raise RuntimeError(f"Reference file not found: {ref_fasta}")
    if os.path.getsize(ref_fasta) == 0:
        raise RuntimeError("Reference file is empty.")
    try:
        with open(ref_fasta, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                raise RuntimeError(
                    "The reference file does not appear to be a valid FASTA file. "
                    "FASTA files must start with '>' followed by a header. "
                    "Please upload a plain‑text FASTA (.fasta, .fa) file."
                )
            f.seek(0)
            ref_records = list(SeqIO.parse(f, "fasta"))
            if not ref_records:
                raise RuntimeError(
                    "Reference file contains no sequences. "
                    "Please check that the file is a valid FASTA with at least one sequence."
                )
    except UnicodeDecodeError:
        raise RuntimeError(
            "The reference file appears to be a binary file (e.g., AB1 trace). "
            "Please upload a plain‑text FASTA (.fasta, .fa) file as the reference."
        )
    except Exception as e:
        raise RuntimeError(f"Could not read reference file: {str(e)}")

    if not os.path.exists(ref_fasta + ".bwt"):
        subprocess.run(["bwa", "index", ref_fasta], check=True, capture_output=True)

    sample = os.path.basename(r1_fastq).replace('_R1.fastq', '')
    sam_file = os.path.join(results_folder, f"{sample}.sam")
    bam_file = os.path.join(results_folder, f"{sample}.bam")
    sorted_bam = os.path.join(results_folder, f"{sample}.sorted.bam")

    # Alignment (single-end)
    with open(sam_file, 'w') as out:
        subprocess.run(['bwa', 'mem', '-M', '-R', f'@RG\\tID:{sample}\\tSM:{sample}',
                        ref_fasta, r1_fastq], stdout=out, stderr=subprocess.PIPE, check=True)
    with open(bam_file, 'w') as out:
        subprocess.run(['samtools', 'view', '-bS', sam_file], stdout=out, check=True)
    subprocess.run(['samtools', 'sort', '-o', sorted_bam, bam_file], check=True)
    subprocess.run(['samtools', 'index', sorted_bam], check=True)

    # Generate consensus using a single pipeline (mpileup → call → consensus)
    consensus_file = os.path.join(results_folder, f"{sample}_consensus.fa")
    with open(consensus_file, 'w') as out:
        # mpileup with bcftools (no legacy -u option)
        mpileup = subprocess.Popen(
            ['bcftools', 'mpileup', '-Ou', '-f', ref_fasta, sorted_bam],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # call variants
        call = subprocess.Popen(
            ['bcftools', 'call', '-c'],
            stdin=mpileup.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # generate consensus
        consensus = subprocess.run(
            ['bcftools', 'consensus', '-f', ref_fasta],
            stdin=call.stdout, stdout=out, stderr=subprocess.PIPE, check=False
        )
        # Check for errors
        mpileup.stdout.close()
        call.stdout.close()
        mpileup.wait()
        call.wait()
        if consensus.returncode != 0:
            raise RuntimeError(f"bcftools consensus failed: {consensus.stderr.decode().strip()}")

    with open(consensus_file) as f:
        contig = f.read().strip()

    # Coverage statistics
    avg_depth, cov_20x = 0, 0
    try:
        depth_out = subprocess.run(['samtools', 'depth', sorted_bam], capture_output=True, text=True)
        if depth_out.returncode == 0:
            depths = [int(l.split()[2]) for l in depth_out.stdout.strip().split('\n') if l]
            if depths:
                avg_depth = sum(depths) / len(depths)
                cov_20x = sum(1 for d in depths if d >= 20) / len(depths) * 100
    except:
        pass

    return {
        'contig': contig, 'variants': [], 'avg_depth': avg_depth,
        'coverage_20x': cov_20x,
        'placements': [{'seq': contig, 'start': 0, 'end': len(contig)-1, 'mismatches': []}]
    }

# ---------- Master Assembly Function ----------
def assemble_reads_from_files(filepaths, mode='greedy', min_overlap=3, max_mismatch=3, k=3,
                              orientation='auto', ref_fasta=None, results_folder=None,
                              trim_low_quality=False, quality_threshold=20, window_size=5):
    all_trim_info = []
    all_labels = []
    trim_details = []

    if ref_fasta:
        reads, _, labels, _ = parse_uploaded_files(filepaths, trim_low_quality, quality_threshold, window_size)
        all_labels = labels
        r1_path = os.path.join(results_folder, 'combined_R1.fastq')
        with open(r1_path, 'w') as f1:
            for i, seq in enumerate(reads):
                qual = ''.join(chr(40+33) for _ in seq)
                f1.write(f"@read{i}\n{seq}\n+\n{qual}\n")
        return guided_assemble_fastq(r1_path, ref_fasta, results_folder)

    reads, all_trim_info, all_labels, trim_details = parse_uploaded_files(
        filepaths, trim_low_quality, quality_threshold, window_size)

    idx_to_label = {i: label for i, label in enumerate(all_labels)}
    clusters = cluster_reads(reads, min_overlap, max_mismatch)
    all_used_indices = set()

    if len(clusters) > 1:
        contigs = []
        cluster_read_indices = []
        for indices in clusters:
            clust_reads = [reads[i] for i in indices]
            if mode == 'greedy':
                assembly = greedy_assemble(clust_reads, min_overlap, 50, max_mismatch, orientation)
                contig = assembly['contig']
            elif mode == 'alignment':
                contig = assemble_with_alignment(clust_reads, min_overlap, orientation)
            elif mode == 'debruijn':
                contig = debruijn_assemble(clust_reads, k)
            else:
                contig = ''
            if contig:
                contigs.append(contig)
                cluster_read_indices.append(indices)
                all_used_indices.update(indices)

        merged_contigs, merge_map = merge_clusters(contigs, min_overlap=min_overlap*2, max_mismatch=max_mismatch)

        all_contigs = []
        for i, contig in enumerate(merged_contigs):
            contributing = merge_map.get(i, [i+1])
            all_reads_used = []
            for cluster_id in contributing:
                orig_idx = cluster_id - 1
                if orig_idx < len(cluster_read_indices):
                    all_reads_used.extend(cluster_read_indices[orig_idx])
            read_names = [idx_to_label[idx] for idx in sorted(set(all_reads_used))]
            all_contigs.append({
                'cluster_id': i + 1,
                'contig': contig,
                'reads_used': len(all_reads_used),
                'read_names': read_names,
                'length': len(contig),
                'merged_from': contributing if len(contributing) > 1 else None
            })

        unused_indices = [i for i in range(len(reads)) if i not in all_used_indices]
        unused_labels = [idx_to_label[i] for i in unused_indices]

        return {
            'contig': all_contigs[0]['contig'] if all_contigs else '',
            'all_contigs': all_contigs,
            'total_clusters': len(all_contigs),
            'initial_clusters': len(clusters),
            'reads_count': len(reads),
            'mode': mode,
            'multi_cluster': True,
            'trim_info': all_trim_info,
            'trim_details': trim_details,
            'unused_reads': unused_labels
        }
    else:
        if mode == 'greedy':
            res = greedy_assemble(reads, min_overlap, 50, max_mismatch, orientation)
        elif mode == 'alignment':
            contig = assemble_with_alignment(reads, min_overlap, orientation)
            res = {'contig': contig, 'steps': [], 'placements': [], 'remaining': [], 'total': len(reads)}
        elif mode == 'debruijn':
            contig = debruijn_assemble(reads, k)
            res = {'contig': contig, 'k': k, 'reads_count': len(reads)}
        else:
            res = {'error': 'Unknown mode'}
        res['reads_count'] = len(reads)
        res['multi_cluster'] = False
        res['all_contigs'] = [{'cluster_id': 1, 'contig': res['contig'], 'reads_used': len(reads),
                                'read_names': all_labels,
                                'length': len(res['contig']), 'merged_from': None}]
        res['total_clusters'] = 1
        res['trim_info'] = all_trim_info
        res['trim_details'] = trim_details
        res['unused_reads'] = []
        return res