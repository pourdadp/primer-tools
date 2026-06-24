
"""
QuickNGS - From FASTQ to report in one click.
Mobile-first Flask application.
"""

import os
import time
import random
import threading
import yaml
from flask import Flask, render_template, request, jsonify, session
from werkzeug.utils import secure_filename

# ---------- App Setup ----------
app = Flask(__name__)
app.secret_key = 'quickngs-secret-key-2025'

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500 MB max

ALLOWED_EXTENSIONS = {'fastq', 'fq', 'gz', 'fasta', 'fa', 'fna'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# ---------- Mock Pipeline (نمایشی) ----------
def run_pipeline(config_path, results_folder):
    """
    این تابع بعداً با ابزارهای واقعی (BWA, FreeBayes و...) جایگزین می‌شود.
    الان یک اجرای نمایشی با تأخیر انجام می‌دهد.
    """
    status = {"step": "QC", "progress": 10}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)
    time.sleep(1)

    status = {"step": "Trimming", "progress": 25}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)
    time.sleep(1)

    status = {"step": "Alignment (BWA)", "progress": 50}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)
    time.sleep(2)

    status = {"step": "Variant Calling (FreeBayes)", "progress": 75}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)
    time.sleep(2)

    status = {"step": "Annotation (SnpEff)", "progress": 90}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)
    time.sleep(1)

    # تولید خروجی ساختگی
    variants = [
        {"chr": "chr1", "pos": 123456, "ref": "A", "alt": "G", "gene": "BRCA1"},
        {"chr": "chr1", "pos": 234567, "ref": "C", "alt": "T", "gene": "TP53"},
    ]
    with open(os.path.join(results_folder, "variants.json"), "w") as f:
        json.dump(variants, f)

    status = {"step": "Report", "progress": 100}
    with open(os.path.join(results_folder, "status.json"), "w") as f:
        json.dump(status, f)

# ---------- Routes ----------
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/start', methods=['POST'])
def start_analysis():
    sample_name = request.form.get('sample_name', 'Sample')
    ref_genome = request.form.get('reference', 'hg38')

    # ذخیره فایل‌ها
    r1 = request.files.get('fastq_r1')
    r2 = request.files.get('fastq_r2')

    if not r1 or not r2:
        return "Error: Both FASTQ files are required.", 400

    run_id = f"{sample_name}_{int(time.time())}"
    results_dir = os.path.join('results', run_id)
    os.makedirs(results_dir, exist_ok=True)

    r1_path = os.path.join(results_dir, secure_filename(r1.filename))
    r2_path = os.path.join(results_dir, secure_filename(r2.filename))
    r1.save(r1_path)
    r2.save(r2_path)

    # ساخت config.yaml از تنظیمات فرم
    config = {
        'sample_name': sample_name,
        'reference': ref_genome,
        'min_quality': int(request.form.get('min_quality', 20)),
        'min_depth': int(request.form.get('min_depth', 10)),
        'trim_adapters': request.form.get('trim_adapters') == 'on',
        'generate_plot': request.form.get('generate_plot') == 'on',
    }
    config_path = os.path.join(results_dir, 'config.yaml')
    with open(config_path, 'w') as f:
        yaml.dump(config, f)

    # اجرای pipeline در پس‌زمینه
    thread = threading.Thread(target=run_pipeline, args=(config_path, results_dir))
    thread.start()

    return render_template('results.html', run_id=run_id)

@app.route('/status/<run_id>')
def get_status(run_id):
    status_file = os.path.join('results', run_id, 'status.json')
    if os.path.exists(status_file):
        with open(status_file) as f:
            return jsonify(json.load(f))
    return jsonify({"step": "Waiting...", "progress": 0})

@app.route('/results/<run_id>')
def view_results(run_id):
    variants_file = os.path.join('results', run_id, 'variants.json')
    variants = []
    if os.path.exists(variants_file):
        with open(variants_file) as f:
            variants = json.load(f)
    return render_template('results_final.html', run_id=run_id, variants=variants)

# ---------- Run ----------
if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
