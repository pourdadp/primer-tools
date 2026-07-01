
# quickngs/app.py
# ... (تمام importها و تنظیمات اولیه بدون تغییر)

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

    directions = []
    for i in range(len(files)):
        dir_key = f'direction_{i}'
        directions.append(request.form.get(dir_key, 'auto'))

    if all(d == 'forward' for d in directions):
        orientation = 'forward'
    elif all(d == 'reverse' for d in directions):
        orientation = 'reverse'
    else:
        orientation = 'auto'

    run_id = f"{project_name}_{int(time.time())}"
    results_dir = os.path.join(RESULTS_FOLDER, run_id)
    os.makedirs(results_dir, exist_ok=True)

    filepaths = []
    for file in files:
        fpath = os.path.join(results_dir, secure_filename(file.filename))
        file.save(fpath)
        filepaths.append(fpath)

    from assembly import assemble_reads_from_files
    try:
        result = assemble_reads_from_files(
            filepaths, mode=algorithm, orientation=orientation,
            results_folder=results_dir,
            trim_low_quality=trim_quality,
            quality_threshold=quality_threshold,
            window_size=window_size
        )

        # Filter contigs
        filtered_contigs = [c for c in result['all_contigs'] if c['reads_used'] > 1]
        extra_unused = [read for ctg in result['all_contigs'] if ctg['reads_used'] <= 1 for read in ctg['read_names']]
        all_unused = result.get('unused_reads', []) + extra_unused
        max_reads = max((c['reads_used'] for c in filtered_contigs), default=1)

        result_path = os.path.join(results_dir, 'result.json')
        with open(result_path, 'w') as f:
            json.dump(result, f)

        logger.info(f"Sanger assembly completed: {run_id}")
        return render_template('assemble_result.html',
                               result=result,
                               filtered_contigs=filtered_contigs,
                               all_unused=all_unused,
                               max_reads=max_reads,
                               filename=files[0].filename,
                               run_id=run_id)
    except ValueError as e:
        logger.error(f"Sanger Assembly Error: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 400
    except RuntimeError as e:
        logger.error(f"Sanger Assembly Error: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500
    except Exception as e:
        logger.error(f"Sanger Assembly Error: {str(e)}")
        return jsonify({"status": "error", "message": f"Unexpected error: {str(e)}"}), 500

@app.route('/assemble_guided', methods=['POST'])
def assemble_guided():
    files = request.files.getlist('seq_files')
    ref_file = request.files.get('ref_file')
    if not files or not ref_file:
        return jsonify({"status": "error", "message": "Both sequence and reference files required."}), 400

    project_name = request.form.get('project_name', 'Guided_Project')
    trim_quality = request.form.get('trim_quality') == 'on'
    quality_threshold = int(request.form.get('quality_threshold', 20))
    window_size = int(request.form.get('window_size', 5))

    directions = []
    for i in range(len(files)):
        dir_key = f'direction_{i}'
        directions.append(request.form.get(dir_key, 'auto'))
    orientation = 'auto'

    run_id = f"{project_name}_{int(time.time())}"
    results_dir = os.path.join(RESULTS_FOLDER, run_id)
    os.makedirs(results_dir, exist_ok=True)

    ref_path = os.path.join(results_dir, secure_filename(ref_file.filename))
    ref_file.save(ref_path)

    filepaths = []
    for file in files:
        fpath = os.path.join(results_dir, secure_filename(file.filename))
        file.save(fpath)
        filepaths.append(fpath)

    from assembly import assemble_reads_from_files
    try:
        result = assemble_reads_from_files(
            filepaths, ref_fasta=ref_path, results_folder=results_dir,
            trim_low_quality=trim_quality,
            quality_threshold=quality_threshold,
            window_size=window_size
        )

        # Filter contigs
        filtered_contigs = [c for c in result['all_contigs'] if c['reads_used'] > 1]
        extra_unused = [read for ctg in result['all_contigs'] if ctg['reads_used'] <= 1 for read in ctg['read_names']]
        all_unused = result.get('unused_reads', []) + extra_unused
        max_reads = max((c['reads_used'] for c in filtered_contigs), default=1)

        result_path = os.path.join(results_dir, 'result.json')
        with open(result_path, 'w') as f:
            json.dump(result, f)

        logger.info(f"Guided assembly completed: {run_id}")
        return render_template('assemble_result.html',
                               result=result,
                               filtered_contigs=filtered_contigs,
                               all_unused=all_unused,
                               max_reads=max_reads,
                               filename=files[0].filename,
                               run_id=run_id)
    except RuntimeError as e:
        logger.error(f"Guided Assembly Error: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500
    except Exception as e:
        logger.error(f"Guided Assembly Error: {str(e)}")
        return jsonify({"status": "error", "message": f"Unexpected error: {str(e)}"}), 500

# ---------- Standard Routes ----------
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/start_ngs', methods=['POST'])
def start_ngs():
    sample_name = request.form.get('project_name', 'NGS_Sample').strip()
    r1 = request.files.get('fastq_r1')
    r2 = request.files.get('fastq_r2')
    if not r1 or not r2:
        return jsonify({"status": "error", "message": "Both FASTQ files required."}), 400

    run_id = f"{sample_name}_{int(time.time())}"
    results_dir = os.path.join(RESULTS_FOLDER, run_id)
    os.makedirs(results_dir, exist_ok=True)

    r1.save(os.path.join(results_dir, secure_filename(r1.filename)))
    r2.save(os.path.join(results_dir, secure_filename(r2.filename)))

    ref_name = request.form.get('reference', 'hg38')
    if ref_name == 'custom':
        custom_ref = request.files.get('custom_reference')
        if custom_ref:
            ref_path = os.path.join(results_dir, secure_filename(custom_ref.filename))
            custom_ref.save(ref_path)
            ref_fasta = ref_path
            if not os.path.exists(ref_fasta + ".bwt"):
                subprocess.run(["bwa", "index", ref_fasta], capture_output=True, check=False)
        else:
            ref_fasta = get_reference_genome('custom', results_dir)
    else:
        ref_fasta = get_reference_genome(ref_name, results_dir)

    config = {
        'sample_name': sample_name,
        'reference': ref_name,
        'ref_fasta': ref_fasta,
        'min_quality': int(request.form.get('min_quality', 20)),
        'min_depth': int(request.form.get('min_depth', 10)),
        'trim_adapters': request.form.get('trim_adapters') == 'on',
    }
    config_path = os.path.join(results_dir, 'config.yaml')
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    threading.Thread(target=run_pipeline, args=(config_path, results_dir)).start()
    return jsonify({"status": "started", "run_id": run_id})

@app.route('/status/<run_id>')
def get_status(run_id):
    status_file = os.path.join(RESULTS_FOLDER, run_id, 'status.json')
    if os.path.exists(status_file):
        with open(status_file) as f:
            return jsonify(json.load(f))
    return jsonify({"step": "Waiting...", "progress": 0, "message": "Initializing..."})

@app.route('/results/<run_id>')
def view_results(run_id):
    result_json = os.path.join(RESULTS_FOLDER, run_id, 'result.json')
    report_json = os.path.join(RESULTS_FOLDER, run_id, 'report.json')

    # اولویت با Sanger (result.json)
    if os.path.exists(result_json):
        try:
            with open(result_json) as f:
                result = json.load(f)

            # محاسبهٔ filtered_contigs و سایر متغیرها
            filtered_contigs = [c for c in result['all_contigs'] if c['reads_used'] > 1]
            extra_unused = [read for ctg in result['all_contigs'] if ctg['reads_used'] <= 1 for read in ctg['read_names']]
            all_unused = result.get('unused_reads', []) + extra_unused
            max_reads = max((c['reads_used'] for c in filtered_contigs), default=1)

            return render_template('assemble_result.html',
                                   result=result,
                                   filtered_contigs=filtered_contigs,
                                   all_unused=all_unused,
                                   max_reads=max_reads,
                                   run_id=run_id)
        except Exception as e:
            logger.error(f"Error loading Sanger result from History: {str(e)}")
            return f"Error reading Sanger result: {str(e)}", 500

    # سپس NGS (report.json)
    if os.path.exists(report_json):
        try:
            with open(report_json) as f:
                report = json.load(f)
            variants_file = os.path.join(RESULTS_FOLDER, run_id, 'variants.json')
            variants = []
            if os.path.exists(variants_file):
                with open(variants_file) as f:
                    variants = json.load(f)
            return render_template('results_final.html', run_id=run_id, report=report, variants=variants)
        except Exception as e:
            logger.error(f"Error loading NGS result from History: {str(e)}")
            return f"Error reading NGS result: {str(e)}", 500

    return jsonify({"status": "error", "message": "No result data found for this run."}), 404

@app.route('/history')
def history():
    runs = []
    if os.path.exists(RESULTS_FOLDER):
        for folder in os.listdir(RESULTS_FOLDER):
            status_file = os.path.join(RESULTS_FOLDER, folder, 'status.json')
            if os.path.exists(status_file):
                with open(status_file) as f:
                    status = json.load(f)
                runs.append({'id': folder, 'step': status.get('step', 'Unknown'), 'progress': status.get('progress', 0)})
    runs.sort(key=lambda x: x['id'], reverse=True)
    return render_template('history.html', runs=runs)

# ... (بقیه endpointها بدون تغییر: /download, /delete, /translate, /blast, /msa, /report_pdf, /log, /about, /help)

if __name__ == '__main__':
    print("=" * 50)
    print("🧬 QuickNGS v2.0 – NGS + Assembly Suite")
    print("Powered by Pourdad Panahi – Built with DeepSeek AI")
    print(f"Uploads: {UPLOAD_FOLDER}")
    print(f"Results: {RESULTS_FOLDER}")
    print("=" * 50)
    app.run(host='0.0.0.0', port=5002, debug=False)