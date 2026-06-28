# quickngs/assembly.py
# ... (تمام توابع قبلی بدون تغییر)

# ---------- Master Assembly Function ----------
def assemble_reads_from_files(filepaths, mode='greedy', min_overlap=3, max_mismatch=3, k=3,
                              orientation='auto', ref_fasta=None, results_folder=None,
                              trim_low_quality=False, quality_threshold=20, window_size=5):
    if ref_fasta:
        reads = parse_uploaded_files(filepaths, trim_low_quality, quality_threshold, window_size)
        r1_path = os.path.join(results_folder, 'combined_R1.fastq')
        r2_path = os.path.join(results_folder, 'combined_R2.fastq')
        with open(r1_path, 'w') as f1, open(r2_path, 'w') as f2:
            for i, seq in enumerate(reads):
                qual = ''.join(chr(40+33) for _ in seq)
                f1.write(f"@read{i}\n{seq}\n+\n{qual}\n")
                f2.write(f"@read{i}_R2\n\n+\n\n")
        return guided_assemble_fastq(r1_path, r2_path, ref_fasta, results_folder)

    reads = parse_uploaded_files(filepaths, trim_low_quality, quality_threshold, window_size)
    clusters = cluster_reads(reads, min_overlap, max_mismatch)

    if len(clusters) > 1:
        contigs = []
        cluster_read_indices = []
        for indices in clusters:
            clust_reads = [reads[i] for i in indices]               # <-- اصلاح‌شده
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

        merged_contigs, merge_map = merge_clusters(contigs, min_overlap=min_overlap*2, max_mismatch=max_mismatch)

        all_contigs = []
        for i, contig in enumerate(merged_contigs):
            contributing = merge_map.get(i, [i+1])
            all_reads_used = []
            for cluster_id in contributing:
                orig_idx = cluster_id - 1
                if orig_idx < len(cluster_read_indices):
                    all_reads_used.extend(cluster_read_indices[orig_idx])
            all_contigs.append({
                'cluster_id': i + 1,
                'contig': contig,
                'reads_used': len(all_reads_used),
                'read_names': [f"read_{r+1}" for r in sorted(set(all_reads_used))],
                'length': len(contig),
                'merged_from': contributing if len(contributing) > 1 else None
            })

        return {
            'contig': all_contigs[0]['contig'] if all_contigs else '',
            'all_contigs': all_contigs,
            'total_clusters': len(all_contigs),
            'initial_clusters': len(clusters),
            'reads_count': len(reads),
            'mode': mode,
            'multi_cluster': True
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
                                'read_names': [f"read_{i+1}" for i in range(len(reads))],
                                'length': len(res['contig']), 'merged_from': None}]
        res['total_clusters'] = 1
        return res
