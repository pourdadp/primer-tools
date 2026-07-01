// static/script.js
const seqInput = document.getElementById('seqFiles');
const refInput = document.getElementById('refFile');
const projectInput = document.getElementById('projectName');
const dynamicDiv = document.getElementById('dynamicSettings');
const startBtn = document.getElementById('startBtn');
let detectedMode = null;

// Tooltip toggle
document.addEventListener('click', function(e) {
    const icon = e.target.closest('.help-icon');
    if (icon) {
        e.stopPropagation();
        const tooltip = icon.querySelector('.help-tooltip');
        const allTips = document.querySelectorAll('.help-tooltip');
        const isVisible = tooltip.style.display === 'block';
        allTips.forEach(t => t.style.display = 'none');
        tooltip.style.display = isVisible ? 'none' : 'block';
    } else {
        document.querySelectorAll('.help-tooltip').forEach(t => t.style.display = 'none');
    }
});

seqInput.addEventListener('change', function() {
    const files = Array.from(this.files);
    if (files.length === 0) {
        resetUI();
        return;
    }

    const hasFastq = files.some(f => /\.(fastq|fq|gz)$/i.test(f.name));
    const hasSanger = files.some(f => /\.(ab1|seq|fasta|fa)$/i.test(f.name));

    if (hasFastq && hasSanger) {
        alert('Mixed file types detected. Please upload only NGS (FASTQ) or only Sanger (AB1/seq/FASTA) files.');
        resetUI();
        return;
    }

    if (hasFastq) {
        if (files.length !== 2) {
            alert('For NGS analysis, please upload exactly two paired FASTQ files (R1 and R2).');
            resetUI();
            return;
        }
        detectedMode = 'ngs';
        buildNGSForm(files);
    } else {
        detectedMode = 'sanger';
        buildSangerForm(files);
    }

    if (!projectInput.value.trim()) {
        const base = files[0].name.replace(/\.[^/.]+$/, "");
        const date = new Date().toISOString().slice(0,10).replace(/-/g,'');
        projectInput.value = base + '_' + date;
    }

    startBtn.disabled = false;
});

refInput.addEventListener('change', function() {
    if (detectedMode === 'sanger') {
        buildSangerForm(Array.from(seqInput.files));
    }
});

function resetUI() {
    dynamicDiv.style.display = 'none';
    dynamicDiv.innerHTML = '';
    startBtn.disabled = true;
    detectedMode = null;
}

function buildNGSForm(files) {
    dynamicDiv.style.display = 'block';
    dynamicDiv.innerHTML = `
        <div class="card mb-2 shadow-sm border-0">
            <div class="card-header bg-primary text-white py-1"><small>🔬 Reference‑Based NGS Settings</small></div>
            <div class="card-body py-2 px-3">
                <div class="row g-2">
                    <div class="col-md-6 mb-2">
                        <label class="form-label small mb-1">🧬 Reference Genome</label>
                        <select class="form-select form-select-sm" name="reference" id="refGenome">
                            <option value="hg38">Human (hg38)</option>
                            <option value="hg19">Human (hg19)</option>
                            <option value="mm10">Mouse (mm10)</option>
                            <option value="custom">Custom FASTA</option>
                        </select>
                    </div>
                    <div class="col-md-3 mb-2">
                        <label class="form-label small mb-1">
                            Min Quality (Q)
                            <span class="help-icon">?
                                <span class="help-tooltip">Minimum Phred quality for trimming. Q20 = 1% error rate. Lower values keep more bases.</span>
                            </span>
                        </label>
                        <input type="number" class="form-control form-control-sm" name="min_quality" value="20">
                    </div>
                    <div class="col-md-3 mb-2">
                        <label class="form-label small mb-1">
                            Min Depth
                            <span class="help-icon">?
                                <span class="help-tooltip">Minimum number of reads covering a position to call a variant. Higher values increase confidence.</span>
                            </span>
                        </label>
                        <input type="number" class="form-control form-control-sm" name="min_depth" value="10">
                    </div>
                </div>
                <div class="form-check form-switch">
                    <input class="form-check-input" type="checkbox" name="trim_adapters" id="trimAdapters" checked>
                    <label class="form-check-label small" for="trimAdapters">
                        Trim Adapters
                        <span class="help-icon">?
                            <span class="help-tooltip">Remove sequencing adapters and low‑quality bases from the ends of reads (Trimmomatic).</span>
                        </span>
                    </label>
                </div>
            </div>
        </div>
    `;
    document.getElementById('refGenome').addEventListener('change', function() {
        const customDiv = document.getElementById('customRefDiv');
        if (!customDiv) {
            const div = document.createElement('div');
            div.id = 'customRefDiv';
            div.innerHTML = '<input type="file" class="form-control form-control-sm mt-2" name="custom_reference" accept=".fasta,.fa,.fna">';
            this.parentNode.appendChild(div);
        }
        customDiv.style.display = this.value === 'custom' ? 'block' : 'none';
    });
}

function buildSangerForm(files) {
    const refProvided = refInput.files.length > 0;
    dynamicDiv.style.display = 'block';
    let fileRows = files.map((f, idx) => `
        <div class="row g-2 align-items-center mb-1">
            <div class="col-7"><small>${f.name}</small></div>
            <div class="col-5">
                <select class="form-select form-select-sm" name="direction_${idx}">
                    <option value="auto">Auto</option>
                    <option value="forward">Forward</option>
                    <option value="reverse">Reverse</option>
                </select>
            </div>
        </div>
    `).join('');

    const hasAB1 = files.some(f => f.name.toLowerCase().endsWith('.ab1'));

    dynamicDiv.innerHTML = `
        <div class="card mb-2 shadow-sm border-0">
            <div class="card-header bg-success text-white py-1"><small>🧬 Sanger Analysis</small></div>
            <div class="card-body py-2 px-3">
                <p class="small mb-2">
                    <strong>${files.length} read(s)</strong> detected.
                    ${refProvided ? 'Reference provided → <span class="badge bg-info">Guided Assembly</span>' : 'No reference → <span class="badge bg-warning text-dark">De Novo Assembly</span>'}
                </p>
                
                <div class="mb-2">
                    <label class="form-label small fw-bold mb-1">Direction (per read)</label>
                    ${fileRows}
                </div>

                ${hasAB1 ? `
                <div class="card bg-light mb-2">
                    <div class="card-body py-2 px-2">
                        <div class="form-check form-switch mb-1">
                            <input class="form-check-input" type="checkbox" name="trim_quality" id="trimQuality" checked>
                            <label class="form-check-label small" for="trimQuality">
                                ✂️ Trim low‑quality ends (Phred scores)
                                <span class="help-icon">?
                                    <span class="help-tooltip">If checked, the beginning and end of the sequence are trimmed when their average quality (over the window) falls below the threshold.</span>
                                </span>
                            </label>
                        </div>
                        <div id="qualityParams" class="row g-2">
                            <div class="col-6">
                                <label class="form-label small mb-1">
                                    Min Quality (Q)
                                    <span class="help-icon">?
                                        <span class="help-tooltip">Minimum acceptable Phred quality. Q20 = 1% error. Lower values keep more bases.</span>
                                    </span>
                                </label>
                                <input type="number" class="form-control form-control-sm" name="quality_threshold" value="20" min="5" max="40">
                            </div>
                            <div class="col-6">
                                <label class="form-label small mb-1">
                                    Window Size
                                    <span class="help-icon">?
                                        <span class="help-tooltip">Number of consecutive bases averaged to decide trimming. Larger windows smooth out random quality drops.</span>
                                    </span>
                                </label>
                                <input type="number" class="form-control form-control-sm" name="window_size" value="5" min="1" max="20">
                            </div>
                        </div>
                    </div>
                </div>
                ` : ''}

                ${!refProvided ? `
                <div class="mb-2">
                    <label class="form-label small fw-bold mb-1">🧠 Assembly Algorithm</label>
                    <select class="form-select form-select-sm" name="algorithm" id="algorithmSelect">
                        <option value="greedy">Greedy OLC (Fast)</option>
                        <option value="alignment">Semi‑Global Alignment (Precise)</option>
                        <option value="debruijn">De Bruijn Graph</option>
                    </select>
                    <div id="debruijnWarning" class="alert alert-warning py-2 px-2 small mt-2" style="display:none;">
                        <strong>⚠️ Note:</strong> De Bruijn Graph is optimized for <em>short NGS reads</em>. For Sanger data, <strong>Greedy</strong> or <strong>Semi‑Global</strong> may produce longer, more accurate contigs.
                    </div>
                </div>` : ''}
            </div>
        </div>
    `;
    
    const trimCheckbox = document.getElementById('trimQuality');
    if (trimCheckbox) {
        trimCheckbox.addEventListener('change', function() {
            document.getElementById('qualityParams').style.display = this.checked ? 'flex' : 'none';
        });
    }

    // De Bruijn warning toggle
    const algoSelect = document.getElementById('algorithmSelect');
    const warningDiv = document.getElementById('debruijnWarning');

    function toggleDebruijnWarning() {
        if (algoSelect && warningDiv) {
            warningDiv.style.display = algoSelect.value === 'debruijn' ? 'block' : 'none';
        }
    }

    if (algoSelect) {
        algoSelect.addEventListener('change', toggleDebruijnWarning);
        toggleDebruijnWarning();
    }
}

// ---------- Form submission handling ----------
document.getElementById('smartForm').addEventListener('submit', function(e) {
    e.preventDefault();
    if (!detectedMode) return;
    
    const formData = new FormData(this);
    const seqFiles = document.getElementById('seqFiles').files;
    const refFile = document.getElementById('refFile').files[0];
    if (refFile) formData.append('ref_file', refFile);
    
    let endpoint = '/start';
    if (detectedMode === 'ngs') {
        formData.append('fastq_r1', seqFiles[0]);
        formData.append('fastq_r2', seqFiles[1]);
        endpoint = '/start_ngs';
    } else if (detectedMode === 'sanger') {
        for (let f of seqFiles) formData.append('seq_files', f);
        if (refFile) endpoint = '/assemble_guided';
        else endpoint = '/sanger_assemble';
    }
    
    document.getElementById('progress_section').style.display = 'block';
    document.getElementById('progress_section').scrollIntoView({ behavior: 'smooth' });
    
    fetch(endpoint, { method: 'POST', body: formData })
    .then(response => {
        if (response.headers.get('content-type') && response.headers.get('content-type').includes('application/json')) {
            return response.json().then(data => {
                throw new Error(data.message || 'Server error');
            });
        }
        return response.text();
    })
    .then(html => {
        const blob = new Blob([html], { type: 'text/html' });
        const url = URL.createObjectURL(blob);
        window.location.href = url;
    })
    .catch(err => {
        alert('Error: ' + err.message);
        document.getElementById('progress_section').style.display = 'none';
    });
});