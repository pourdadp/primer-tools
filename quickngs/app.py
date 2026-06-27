# quickngs/app.py
"""
QuickNGS – From FASTQ/AB1 to clinical report or assembled contig.
Real NGS pipeline (BWA + Samtools + FreeBayes + SnpEff) + De Novo / Guided Assembly.
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

# ... (تمامی importها و توابع کمکی بدون تغییر)
# توجه: برای جلوگیری از طولانی شدن، کل فایل تکرار نشده. فقط routeهای اصلاح‌شده در زیر آمده‌اند.

# ---------- Assembly routes (اصلاح‌شده) ----------
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

    # محاسبهٔ orientation بر اساس direction انتخاب‌شده برای هر خوانش
    directions = []
    for i in range(len(files)):
        dir_key = f'direction_{i}'
        directions.append(request.form.get(dir_key, 'auto'))

    # تعیین orientation کلی
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
    result = assemble_reads_from_files(
        filepaths, mode=algorithm, orientation=orientation,
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

    # بررسی وجود bcftools (علاوه بر ابزارهای دیگر که در guided_assemble_fastq چک می‌شود)
    if not is_tool_available('bcftools'):
        return jsonify({"status": "error", "message": "bcftools is not installed. Please install it (sudo apt install bcftools)."}), 500

    project_name = request.form.get('project_name', 'Guided_Project')
    trim_quality = request.form.get('trim_quality') == 'on'
    quality_threshold = int(request.form.get('quality_threshold', 20))
    window_size = int(request.form.get('window_size', 5))

    # orientation برای guided هم می‌تواند استفاده شود (اگرچه فعلاً تأثیری ندارد)
    directions = []
    for i in range(len(files)):
        dir_key = f'direction_{i}'
        directions.append(request.form.get(dir_key, 'auto'))
    orientation = 'auto'  # guided از orientation استفاده نمی‌کند

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
    result = assemble_reads_from_files(
        filepaths, ref_fasta=ref_path, results_folder=results_dir,
        trim_low_quality=trim_quality,
        quality_threshold=quality_threshold,
        window_size=window_size
    )
    return render_template('assemble_result.html', result=result, filename=files[0].filename)

# ... (بقیه routeها بدون تغییر)
