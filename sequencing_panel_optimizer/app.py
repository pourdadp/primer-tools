# -*- coding: utf-8 -*-
import os
import sys
import json
import shutil
import subprocess
import yaml
import sqlite3
from datetime import datetime
from flask import Flask, render_template, request, redirect, url_for, send_from_directory, jsonify, session, flash

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(BASE_DIR, 'templates')
SCRIPT_DIR = os.path.join(BASE_DIR, 'scripts')
UPLOAD_DIR = os.path.join(BASE_DIR, 'results')
HISTORY_FILE = os.path.join(BASE_DIR, 'history.json')

os.makedirs(TEMPLATE_DIR, exist_ok=True)
os.makedirs(SCRIPT_DIR, exist_ok=True)
os.makedirs(UPLOAD_DIR, exist_ok=True)

app = Flask(__name__, template_folder=TEMPLATE_DIR)
app.secret_key = 'spo-secret-key-change-in-production'
app.config['UPLOAD_FOLDER'] = UPLOAD_DIR

if not os.path.exists(HISTORY_FILE):
    with open(HISTORY_FILE, 'w', encoding='utf-8') as f:
        json.dump([], f)

def save_history(entry):
    with open(HISTORY_FILE, 'r', encoding='utf-8') as f:
        history = json.load(f)
    history.append(entry)
    with open(HISTORY_FILE, 'w', encoding='utf-8') as f:
        json.dump(history, f, indent=2, ensure_ascii=False)

def get_history():
    with open(HISTORY_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)

# ==================== پیدا کردن دیتابیس PDM ====================

def find_pdm_database():
    possible_paths = [
        os.path.join(os.path.dirname(BASE_DIR), 'primer_database_manager', 'primers.db'),
        os.path.join(BASE_DIR, 'primers.db'),
        os.path.join(os.path.dirname(os.path.dirname(BASE_DIR)), 'primer_database_manager', 'primers.db'),
        r'C:\Users\user\Desktop\primer-tools-main\primer_database_manager\primers.db',
        r'C:\Users\user\Desktop\primer-tools\primer_database_manager\primers.db',
    ]
    for path in possible_paths:
        if os.path.exists(path):
            return path
    return None

# ==================== Routes ====================

@app.route('/', methods=['GET', 'POST'])
def index():
    selected_primers = session.pop('selected_primers', None)
    
    if request.method == 'POST':
        dna_sequence = request.form['dna_sequence']
        max_mispairing = int(request.form['max_mispairing'])
        max_product_length = int(request.form['max_product_length'])
        min_product_length = int(request.form.get('min_product_length', 0))
        max_tm_diff = float(request.form.get('max_tm_diff', 5.0))
        
        primers_text = request.form['primers']
        primers = []
        for line in primers_text.strip().split('\n'):
            if line.strip():
                parts = line.split(',')
                if len(parts) == 2:
                    primers.append({
                        'name': parts[0].strip(),
                        'sequence': parts[1].strip()
                    })
        
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        run_dir = os.path.join(UPLOAD_DIR, f'run_{timestamp}')
        os.makedirs(run_dir, exist_ok=True)
        os.makedirs(os.path.join(run_dir, 'config'), exist_ok=True)
        os.makedirs(os.path.join(run_dir, 'results'), exist_ok=True)
        os.makedirs(os.path.join(run_dir, 'scripts'), exist_ok=True)
        
        config = {
            'dna_sequence': dna_sequence,
            'primers': primers,
            'max_mispairing': max_mispairing,
            'max_product_length': max_product_length,
            'min_product_length': min_product_length,
            'max_tm_diff': max_tm_diff
        }
        config_path = os.path.join(run_dir, 'config', 'config.yaml')
        with open(config_path, 'w', encoding='utf-8') as f:
            yaml.dump(config, f, allow_unicode=True)
        
        if os.path.exists(SCRIPT_DIR):
            for file in os.listdir(SCRIPT_DIR):
                if file.endswith('.py'):
                    src = os.path.join(SCRIPT_DIR, file)
                    dst = os.path.join(run_dir, 'scripts', file)
                    shutil.copy2(src, dst)
        
        scripts_order = [
            'parse_input.py',
            'check_binding.py',
            'check_dimers.py',
            'find_pairs.py',
            'visualize.py',
            'best_coverage.py'
        ]
        
        log_file = os.path.join(run_dir, 'run.log')
        success = True
        with open(log_file, 'w', encoding='utf-8') as log:
            for script in scripts_order:
                script_path = os.path.join(run_dir, 'scripts', script)
                if not os.path.exists(script_path):
                    log.write(f"WARNING: Script {script} not found. Skipping.\n")
                    continue
                log.write(f"\n--- Running {script} ---\n")
                log.flush()
                cmd = ['python', script_path]
                proc = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True, encoding='utf-8')
                log.write(proc.stdout)
                log.write(proc.stderr)
                if proc.returncode != 0:
                    log.write(f"ERROR: {script} failed with code {proc.returncode}\n")
                    success = False
                    break
        
        entry = {
            'run_id': timestamp,
            'date': timestamp,
            'dna_length': len(dna_sequence),
            'num_primers': len(primers),
            'max_mispairing': max_mispairing,
            'max_product_length': max_product_length,
            'min_product_length': min_product_length,
            'max_tm_diff': max_tm_diff,
            'status': 'success' if success else 'failed',
        }
        save_history(entry)
        
        return redirect(url_for('results', run_id=timestamp))
    
    return render_template('index.html', selected_primers=selected_primers)

@app.route('/results/<run_id>')
def results(run_id):
    run_dir = os.path.join(UPLOAD_DIR, f'run_{run_id}')
    try:
        with open(os.path.join(run_dir, 'results', 'binding_results.json'), encoding='utf-8') as f:
            binding = json.load(f)
        with open(os.path.join(run_dir, 'results', 'valid_pairs.json'), encoding='utf-8') as f:
            valid_pairs = json.load(f)
        with open(os.path.join(run_dir, 'results', 'best_coverage.json'), encoding='utf-8') as f:
            best = json.load(f)
        with open(os.path.join(run_dir, 'results', 'pair_info.json'), encoding='utf-8') as f:
            pair_info = json.load(f)
    except FileNotFoundError:
        return "Results not found", 404
    
    return render_template('results.html', 
                           binding=binding,
                           valid_pairs=valid_pairs,
                           best=best,
                           pair_info=pair_info,
                           run_id=run_id)

@app.route('/history')
def history():
    history_data = get_history()
    return render_template('history.html', history=history_data)

@app.route('/download/<run_id>/<filename>')
def download(run_id, filename):
    run_dir = os.path.join(UPLOAD_DIR, f'run_{run_id}', 'results')
    return send_from_directory(run_dir, filename)

@app.route('/static/run_<run_id>/<filename>')
def serve_result_image(run_id, filename):
    run_dir = os.path.join(UPLOAD_DIR, f'run_{run_id}', 'results')
    return send_from_directory(run_dir, filename)

@app.route('/sample_data')
def sample_data():
    sample = {
        "dna_sequence": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "primers": "F1,ATCGATCG\nR1,CGATCGAT\nF2,GATCGATC\nR2,TCGATCGA\nF3,ATCGATCGAT\nR3,CGATCGATCG",
        "max_mispairing": 2,
        "max_product_length": 80,
        "min_product_length": 20,
        "max_tm_diff": 5.0
    }
    return jsonify(sample)

# ==================== Load All from Database ====================

@app.route('/load_from_db')
def load_from_db():
    db_path = find_pdm_database()
    if not db_path:
        return jsonify({'error': 'Primer database not found.'}), 404
    
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''
            SELECT name, forward_sequence, reverse_sequence
            FROM primers WHERE is_active = 1
        ''')
        rows = cursor.fetchall()
        conn.close()
        
        primer_lines = []
        for row in rows:
            if row['forward_sequence']:
                primer_lines.append(f"{row['name']}_F,{row['forward_sequence']}")
            if row['reverse_sequence']:
                primer_lines.append(f"{row['name']}_R,{row['reverse_sequence']}")
        
        return jsonify({'primers': '\n'.join(primer_lines)})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# ==================== Select Primers from Database ====================

@app.route('/load_from_db_selector')
def load_from_db_selector():
    db_path = find_pdm_database()
    if not db_path:
        flash('Primer database not found. Please run Primer Database Manager first.', 'danger')
        return redirect(url_for('index'))
    
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''
            SELECT id, name, forward_sequence, reverse_sequence, pair_name, gene, organism
            FROM primers WHERE is_active = 1
            ORDER BY name
        ''')
        primers = cursor.fetchall()
        conn.close()
        
        if not primers:
            flash('No primers found in database.', 'warning')
            return redirect(url_for('index'))
        
        return render_template('primer_selector.html', primers=primers)
    except Exception as e:
        flash(f'Error: {str(e)}', 'danger')
        return redirect(url_for('index'))

@app.route('/add_selected_primers', methods=['POST'])
def add_selected_primers():
    primer_ids = request.form.getlist('primer_ids')
    if not primer_ids:
        flash('No primers selected.', 'warning')
        return redirect(url_for('index'))
    
    db_path = find_pdm_database()
    if not db_path:
        flash('Primer database not found.', 'danger')
        return redirect(url_for('index'))
    
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        placeholders = ','.join('?' * len(primer_ids))
        cursor.execute(f'''
            SELECT name, forward_sequence, reverse_sequence
            FROM primers WHERE id IN ({placeholders}) AND is_active = 1
        ''', primer_ids)
        rows = cursor.fetchall()
        conn.close()
        
        if not rows:
            flash('No valid primers found.', 'warning')
            return redirect(url_for('index'))
        
        primer_lines = []
        for row in rows:
            if row['forward_sequence']:
                primer_lines.append(f"{row['name']}_F,{row['forward_sequence']}")
            if row['reverse_sequence']:
                primer_lines.append(f"{row['name']}_R,{row['reverse_sequence']}")
        
        session['selected_primers'] = '\n'.join(primer_lines)
        flash(f'{len(rows)} primer(s) loaded successfully.', 'success')
        return redirect(url_for('index'))
    except Exception as e:
        flash(f'Error: {str(e)}', 'danger')
        return redirect(url_for('index'))

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
