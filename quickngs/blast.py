# quickngs/blast.py
"""
BLAST Module for QuickNGS
Supports Local BLAST+ and NCBI Web BLAST.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

import subprocess
import os
import tempfile
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------- Local BLAST (requires BLAST+ installed) ----------
def run_local_blast(sequence, db_path, program='blastn', evalue=0.001, max_hits=50):
    """
    Run local BLAST+ against a pre-indexed database.
    Requires: makeblastdb and blastn/blastp/blastx installed.
    """
    # Write query to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">query\n{sequence}\n")
        query_file = f.name
    
    try:
        cmd = [
            program,
            '-db', db_path,
            '-query', query_file,
            '-evalue', str(evalue),
            '-max_target_seqs', str(max_hits),
            '-outfmt', '5'  # XML output
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            return {'error': f"BLAST failed: {result.stderr}"}
        
        # Parse XML
        from io import StringIO
        blast_records = NCBIXML.parse(StringIO(result.stdout))
        hits = []
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    hits.append({
                        'title': alignment.title,
                        'length': alignment.length,
                        'score': hsp.score,
                        'evalue': hsp.expect,
                        'identity': round(hsp.identities / hsp.align_length * 100, 1),
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'sbjct_start': hsp.sbjct_start,
                        'sbjct_end': hsp.sbjct_end,
                        'align_length': hsp.align_length
                    })
        return hits[:max_hits]
    except subprocess.TimeoutExpired:
        return {'error': 'BLAST timed out (>5 minutes).'}
    except Exception as e:
        return {'error': str(e)}
    finally:
        os.unlink(query_file)

# ---------- NCBI Web BLAST (requires internet) ----------
def run_web_blast(sequence, database='nt', program='blastn', evalue=0.001, max_hits=50):
    """
    Submit BLAST search to NCBI via Biopython.
    This may take 30-120 seconds depending on queue.
    """
    try:
        record = SeqRecord(Seq(sequence), id="query")
        result_handle = NCBIWWW.qblast(program, database, record.format('fasta'), 
                                        expect=evalue, hitlist_size=max_hits)
        
        blast_records = NCBIXML.parse(result_handle)
        hits = []
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    hits.append({
                        'title': alignment.title,
                        'accession': alignment.accession,
                        'length': alignment.length,
                        'score': hsp.score,
                        'evalue': hsp.expect,
                        'identity': round(hsp.identities / hsp.align_length * 100, 1),
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'sbjct_start': hsp.sbjct_start,
                        'sbjct_end': hsp.sbjct_end,
                        'align_length': hsp.align_length,
                        'link': f"https://www.ncbi.nlm.nih.gov/nuccore/{alignment.accession}"
                    })
        return hits[:max_hits]
    except Exception as e:
        return {'error': f"Web BLAST failed: {str(e)}"}

# ---------- Master BLAST function ----------
def run_blast(sequence, mode='web', db='nt', program='blastn', evalue=0.001, max_hits=50):
    """
    Run BLAST search.
    mode: 'web' (NCBI) or 'local' (requires BLAST+ and indexed database)
    """
    if mode == 'local':
        return run_local_blast(sequence, db, program, evalue, max_hits)
    else:
        return run_web_blast(sequence, db, program, evalue, max_hits)
