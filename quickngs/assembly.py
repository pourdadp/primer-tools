# quickngs/assembly.py
# ... (همه توابع قبلی بدون تغییر)

# ---------- Guided Assembly (with reference validation) ----------
def guided_assemble_fastq(r1_fastq, ref_fasta, results_folder):
    for tool in ['bwa', 'samtools', 'bcftools']:
        if not is_tool_available(tool):
            raise RuntimeError(f"{tool} is not installed.")

    # --- اعتبارسنجی فایل Reference ---
    if not os.path.exists(ref_fasta):
        raise RuntimeError(f"Reference file not found: {ref_fasta}")
    if os.path.getsize(ref_fasta) == 0:
        raise RuntimeError("Reference file is empty.")

    # تلاش برای خواندن به‌عنوان FASTA—اگر خطای رمزگشایی داشت، پیام مناسب بده
    try:
        with open(ref_fasta, 'r', encoding='utf-8') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                raise RuntimeError(
                    "The reference file does not appear to be a valid FASTA file. "
                    "FASTA files must start with '>' followed by a header. "
                    "Please upload a plain‑text FASTA (.fasta, .fa) file."
                )
            # بازخوانی کامل برای SeqIO.parse
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

    # Index reference if needed
    if not os.path.exists(ref_fasta + ".bwt"):
        subprocess.run(["bwa", "index", ref_fasta], check=True, capture_output=True)

    # ... (بقیه کد هم‌مانطور که بود: هم‌ردیفی، consensus و غیره)
