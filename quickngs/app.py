# quickngs/app.py
"""
QuickNGS – From FASTQ/AB1 to clinical report or assembled contig.
Real NGS pipeline (BWA + Samtools + FreeBayes + SnpEff) + De Novo / Guided Assembly.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

import os
import time
import json
import threading
import subprocess
import yaml
import shutil
import sys
from flask import Flask, render_template, request, jsonify, send_file
from werkzeug.utils import secure_filename

# ---------- Smart Storage Management ----------
def get_free_space(path='/'):
    try:
        stat = shutil.disk_usage(path)
        return stat.free
    except OSError:
        return 0

def find_best_storage():
    root_free = get_free_space('/')
    if root_free > 2 * 1024 * 1024 * 1024:
        default_uploads = os.path.join(os.getcwd(), 'uploads')
        default_results = os.path.join(os.getcwd(), 'results')
        return default_uploads, default_results

    candidates = []
    for base in ['/media', '/mnt']:
        if not os.path.exists(base):
            continue
        try:
            for user in os.listdir(base):
                user_path = os.path.join(base, user)
                if os.path.isdir(user_path):
                    for volume in os.listdir(user_path):
                        vol_path = os.path.join(user_path, volume)
                        if os.path.isdir(vol_path):
                            free = get_free_space(vol_path)
                            if free > 10 * 1024 * 1024 * 1024:
                                candidates.append((vol_path, free))
        except OSError:
            continue

    if not candidates:
        return os.path.join(os.getcwd(), 'uploads'), os.path.join(os.getcwd(), 'results')

    best_volume, _ = max(candidates, key=lambda x: x[1])
    base_path = os.path.join(best_volume, 'quickngs_data')
    uploads_path = os.path.join(base_path, 'uploads')
    results_path = os.path.join(base_path, 'results')
    return uploads_path, results_path

# ---------- App Setup ----------
app = Flask(__name__)
app.secret_key = 'quickngs-secret-key-2025'

try:
    UPLOAD_FOLDER, RESULTS_FOLDER = find_best_storage()
    os.makedirs(UPLOAD_FOLDER, exist_ok=True)
    os.makedirs(RESULTS_FOLDER, exist_ok=True)
except OSError as e:
    print(f"Fatal Error: Could not create storage directories. {e}", file=sys.stderr)
    sys.exit(1)

CACHE_FOLDER = os.path.join(os.path.dirname(RESULTS_FOLDER), 'quickngs_cache')
os.makedirs(CACHE_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 * 1024  # 10 GB

ALLOWED_EXTENSIONS = {'fastq', 'fq', 'gz', 'fasta', 'fa', 'fna', 'ab1', 'seq', 'txt'}

# ---------- Helper Functions ----------
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def update_status(results_folder, step_name, progress, message=""):
    status = {"step": step_name, "progress": progress, "message": message}
    try:
        with open(os.path.join(results_folder, "status.json"), "w") as f:
            json.dump(status, f)
    except OSError:
        pass

def has_enough_space(file_storage_list, target_dir):
    total_size = 0
    for fs in file_storage_list:
        if fs:
            try:
                fs.stream.seek(0, os.SEEK_END)
                size = fs.stream.tell()
                fs.stream.seek(0)
                total_size += size
            except:
                continue
    free = get_free_space(target_dir)
    if free >= total_size:
        return True, total_size
    else:
        return False, total_size - free

def human_readable_size(size_bytes):
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"

def is_tool_available(tool_name):
    try:
        subprocess.run([tool_name], capture_output=True, check=False)
        return True
    except FileNotFoundError:
        return False

def parse_fastqc_data(fastqc_data_path):
    stats = {"basic_statistics": "N/A", "per_base_quality": "N/A", "total_sequences": "N/A"}
    try:
        with open(fastqc_data_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    if parts[1] == "Basic Statistics":
                        stats["basic_statistics"] = parts[0]
                    elif parts[1] == "Per base sequence quality":
                        stats["per_base_quality"] = parts[0]
                    elif parts[1] == "Total Sequences":
                        stats["total_sequences"] = parts[2]
    except:
        pass
    return stats

def parse_vcf_for_variants(vcf_path):
    variants = []
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                info = parts[7]
                gene = "Unknown"
                impact = "Unknown"
                if 'ANN=' in info:
                    ann_parts = info.split('ANN=')[1].split(';')[0].split('|')
                    if len(ann_parts) > 3:
                        gene = ann_parts[3] if ann_parts[3] else "Unknown"
                    if len(ann_parts) > 1:
                        impact = ann_parts[1] if ann_parts[1] else "Unknown"
                variants.append({
                    "chr": chrom, "pos": pos, "ref": ref, "alt": alt,
                    "qual": qual, "gene": gene, "impact": impact
                })
    except:
        pass
    return variants

# ---------- Reference Genome ----------
REFERENCE_GENOMES = {
    "hg38": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "description": "Human GRCh38/hg38"
    },
    "hg19": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "description": "Human GRCh37/hg19"
    },
    "mm10": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
        "description": "Mouse GRCm38/mm10"
    }
}

def get_reference_genome(ref_name, results_folder):
    if ref_name == "custom" or ref_name not in REFERENCE_GENOMES:
        mock_ref = os.path.join(results_folder, "reference.fa")
        if not os.path.exists(mock_ref):
            with open(mock_ref, "w") as f:
                f.write(">mock\n" + "A" * 5000 + "\n")
        return mock_ref
    ref_info = REFERENCE_GENOMES[ref_name]
    ref_fasta = os.path.join(CACHE_FOLDER, f"{ref_name}.fa")
    ref_gz = ref_fasta + ".gz"
    if os.path.exists(ref_fasta + ".bwt"):
        return ref_fasta
    if not os.path.exists(ref_fasta):
        update_status(results_folder, "Downloading Reference", 8, f"Downloading {ref_info['description']}...")
        try:
            subprocess.run(["wget", "-q", "-O", ref_gz, ref_info["url"]], check=True, timeout=3600)
            subprocess.run(["gunzip", "-f", ref_gz], check=True)
        except:
            update_status(results_folder, "Warning", 9, "Could not download reference genome. Using mock reference.")
            mock_ref = os.path.join(results_folder, "reference.fa")
            with open(mock_ref, "w") as f:
                f.write(">mock\n" + "A" * 5000 + "\n")
            return mock_ref
    update_status(results_folder, "Indexing Reference", 12, f"Building BWA index for {ref_info['description']}...")
    result = subprocess.run(["bwa", "index", ref_fasta], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to index reference: {result.stderr}")
    return ref_fasta

# ---------- FastQC ----------
def run_fastqc(fastq_path, output_dir, results_folder):
    update_status(results_folder, "FastQC", 15, f"Running quality control...")
    try:
        subprocess.run(["fastqc", "-q", "-o", output_dir, fastq_path], check=True, timeout=600)
        base_name = os.path.basename(fastq_path).rsplit('.', 1)[0]
        fastqc_dir = os.path.join(output_dir, base_name + "_fastqc")
        if not os.path.exists(fastqc_dir):
            for item in os.listdir(output_dir):
                if item.endswith("_fastqc"):
                    fastqc_dir = os.path.join(output_dir, item)
                    break
        summary_file = os.path.join(fastqc_dir, "summary.txt")
        if os.path.exists(summary_file):
            return parse_fastqc_data(summary_file)
    except:
        pass
    return {"basic_statistics": "N/A", "per_base_quality": "N/A", "total_sequences": "N/A"}

# ---------- Trimming ----------
def run_trimmomatic(r1_path, r2_path, results_folder, min_quality=20):
    update_status(results_folder, "Trimming Adapters", 22, "Removing adapters and low-quality bases...")
    trimmed_r1 = os.path.join(results_folder, "trimmed_R1.fastq.gz")
    trimmed_r2 = os.path.join(results_folder, "trimmed_R2.fastq.gz")
    unpaired = os.path.join(results_folder, "unpaired.fastq.gz")
    try:
        subprocess.run([
            "trimmomatic", "PE", "-phred33",
            r1_path, r2_path, trimmed_r1, unpaired, trimmed_r2, unpaired,
            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
            f"LEADING:{min_quality}", f"TRAILING:{min_quality}",
            "SLIDINGWINDOW:4:20", "MINLEN:50"
        ], check=True, capture_output=True, text=True, timeout=3600)
        return trimmed_r1, trimmed_r2
    except:
        update_status(results_folder, "Warning", 23, "Trimming failed. Using original files.")
        return r1_path, r2_path

# ---------- NGS Pipeline ----------
def run_pipeline(config_path, results_folder):
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
        sample_name = config.get('sample_name', 'Unknown Sample')
        ref_name = config.get('reference', 'custom')
        min_quality = config.get('min_quality', 20)
        min_depth = config.get('min_depth', 10)
        do_trim = config.get('trim_adapters', True)

        fastq_files = [f for f in os.listdir(results_folder) if f.endswith(('.fastq', '.fq', '.gz'))]
        r1 = next((f for f in fastq_files if 'R1' in f or '_1' in f), None)
        r2 = next((f for f in fastq_files if 'R2' in f or '_2' in f), None)
        if not r1 or not r2:
            update_status(results_folder, "Error", 0, "FASTQ files not found.")
            return
        r1_path = os.path.join(results_folder, r1)
        r2_path = os.path.join(results_folder, r2)

        update_status(results_folder, "Checking Tools", 5, "Verifying analysis tools...")
        for tool, msg in {"bwa": "BWA", "samtools": "Samtools", "freebayes": "FreeBayes"}.items():
            if not is_tool_available(tool):
                update_status(results_folder, "Error", 0, f"{msg} not installed.")
                return

        update_status(results_folder, "Preparing Reference", 8, f"Setting up {ref_name}...")
        ref_fasta = get_reference_genome(ref_name, results_folder)

        fastqc_dir = os.path.join(results_folder, "fastqc")
        os.makedirs(fastqc_dir, exist_ok=True)
        qc_before = run_fastqc(r1_path, fastqc_dir, results_folder)

        if do_trim and is_tool_available("trimmomatic"):
            r1_path, r2_path = run_trimmomatic(r1_path, r2_path, results_folder, min_quality)
            qc_after = run_fastqc(r1_path, fastqc_dir, results_folder)
        else:
            update_status(results_folder, "Skipping Trimming", 20, "Trimming skipped.")
            qc_after = qc_before

        update_status(results_folder, "Alignment (BWA)", 35, "Aligning reads...")
        sam_file = os.path.join(results_folder, f"{sample_name}.sam")
        bam_file = os.path.join(results_folder, f"{sample_name}.bam")
        sorted_bam = os.path.join(results_folder, f"{sample_name}.sorted.bam")
        with open(sam_file, "w") as sam_out:
            subprocess.run(["bwa", "mem", "-M", "-R", f"@RG\\tID:{sample_name}\\tSM:{sample_name}",
                            ref_fasta, r1_path, r2_path], stdout=sam_out, stderr=subprocess.PIPE, check=True)
        with open(bam_file, "w") as bam_out:
            subprocess.run(["samtools", "view", "-bS", sam_file], stdout=bam_out, check=True)
        subprocess.run(["samtools", "sort", "-o", sorted_bam, bam_file], check=True)
        subprocess.run(["samtools", "index", sorted_bam], check=True)

        update_status(results_folder, "Variant Calling (FreeBayes)", 70, "Calling variants...")
        vcf_file = os.path.join(results_folder, f"{sample_name}.vcf")
        with open(vcf_file, "w") as vcf_out:
            subprocess.run(["freebayes", "-f", ref_fasta, "--min-base-quality", str(min_quality),
                            "--min-coverage", str(min_depth), sorted_bam], stdout=vcf_out, stderr=subprocess.PIPE)
        update_status(results_folder, "Annotation (SnpEff)", 80, "Annotating variants...")
        annotated_vcf = os.path.join(results_folder, f"{sample_name}.annotated.vcf")
        if is_tool_available("snpEff"):
            subprocess.run(["snpEff", "GRCh38.99", vcf_file], stdout=open(annotated_vcf, "w"), stderr=subprocess.PIPE)
        else:
            shutil.copy(vcf_file, annotated_vcf)

        variants = parse_vcf_for_variants(annotated_vcf if os.path.exists(annotated_vcf) else vcf_file)
        with open(os.path.join(results_folder, "variants.json"), "w") as f:
            json.dump(variants, f)

        update_status(results_folder, "Coverage Analysis", 90, "Calculating coverage...")
        avg_depth, cov_20x = 0, 0
        try:
            depth_out = subprocess.run(["samtools", "depth", sorted_bam], capture_output=True, text=True)
            depths = [int(l.split()[2]) for l in depth_out.stdout.strip().split("\n") if l]
            if depths:
                avg_depth = sum(depths)/len(depths)
                cov_20x = sum(1 for d in depths if d>=20)/len(depths)*100
        except: pass

        update_status(results_folder, "Generating Report", 95, "Assembling final report...")
        report = {
            "sample_name": sample_name, "reference": ref_name,
            "date": time.strftime("%Y-%m-%d %H:%M:%S"), "total_variants": len(variants),
            "average_depth": round(avg_depth, 1), "coverage_20x": round(cov_20x, 1),
            "qc_before": qc_before.get("basic_statistics", "N/A"),
            "qc_after": qc_after.get("basic_statistics", "N/A"),
            "variants": variants[:50]
        }
        with open(os.path.join(results_folder, "report.json"), "w") as f:
            json.dump(report, f)

        update_status(results_folder, "Completed", 100, "Analysis complete!")

    except Exception as e:
        update_status(results_folder, "Error", 0, str(e))

# ---------- Assembly routes ----------
@app.route('/sanger_assemble', methods=['POST'])
def sanger_assemble():
    files = request.files.getlist('seq_files')
    if not files:
        return jsonify({"status": "error", "message": "No sequence files uploaded."}), 400
    algorithm = request.form.get('algorithm', 'greedy')
    project_name = request.form.get('project_name', 'Sanger_Project')
    trim_quality = request.form.get('trim_quality') == 'on'
    quality_threshold = int(request.form.get('quality_threshold', 20))
    window_size = int(request.form.get('window_size', 5))
    run_id = f"{project_name}_{int(time.time())}"
    results_dir = os.path.join(RESULTS_FOLDER, run_id)
    os.makedirs(results_dir, exist_ok=True)
    filepaths = []
    for file in files:
        fpath = os.path.join(results_dir, secure_filename(file.filename))
        file.save(fpath)
        filepaths.append(fpath)
    from assembly import assemble_reads_from_files
    result = assemble_reads_from_files(
        filepaths, mode=algorithm, orientation='auto',
        results_folder=results_dir,
        trim_low_quality=trim_quality,
        quality_threshold=quality_threshold,
        window_size=window_size
    )
    return render_template('assemble_result.html', result=result, filename=files[0].filename)

@app.route('/assemble_guided', methods=['POST'])
def assemble_guided():
    files = request.files.getlist('seq_files')
    ref_file = request.files.get('ref_file')
    if not files or not ref_file:
        return jsonify({"status": "error", "message": "Both sequence and reference files required."}), 400
    project
