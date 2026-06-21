# -*- coding: utf-8 -*-
import os
import json
import shutil
import sqlite3
from datetime import datetime, timedelta
from flask import Flask, render_template, request, redirect, url_for, flash, session, jsonify, send_file
from werkzeug.security import generate_password_hash, check_password_hash
from functools import wraps

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = os.path.join(BASE_DIR, 'templates')
STATIC_DIR = os.path.join(BASE_DIR, 'static')
DATABASE_PATH = os.path.join(BASE_DIR, 'primers.db')
BACKUP_DIR = os.path.join(BASE_DIR, 'backups')

os.makedirs(TEMPLATE_DIR, exist_ok=True)
os.makedirs(STATIC_DIR, exist_ok=True)
os.makedirs(BACKUP_DIR, exist_ok=True)

app = Flask(__name__, template_folder=TEMPLATE_DIR, static_folder=STATIC_DIR)
app.secret_key = 'your-secret-key-change-this-in-production'

# ==================== Helper Functions ====================

def calculate_tm(sequence):
    """محاسبه دمای اتصال با فرمول Wallace (2*(A+T) + 4*(G+C))"""
    if not sequence:
        return None
    seq = sequence.upper()
    at = seq.count('A') + seq.count('T')
    gc = seq.count('G') + seq.count('C')
    tm = 2 * at + 4 * gc
    return round(tm, 1)

def get_db():
    conn = sqlite3.connect(DATABASE_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    conn = get_db()
    c = conn.cursor()

    c.execute('''
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            username TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL,
            role TEXT NOT NULL DEFAULT 'viewer',
            is_active BOOLEAN DEFAULT 1,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            last_login DATETIME
        )
    ''')
    # ... (سایر جداول به‌صورت کامل در نسخه قبلی موجود است)

    admin = c.execute("SELECT * FROM users WHERE username = 'admin'").fetchone()
    if not admin:
        admin_hash = generate_password_hash('admin123')
        c.execute(
            "INSERT INTO users (username, password_hash, role) VALUES (?, ?, ?)",
            ('admin', admin_hash, 'admin')
        )

    c.execute('''
        CREATE TRIGGER IF NOT EXISTS update_primer_timestamp
        AFTER UPDATE ON primers
        BEGIN
            UPDATE primers SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
        END
    ''')

    conn.commit()
    conn.close()

init_db()

# ==================== Database Helper Functions ====================

def get_user_by_username(username):
    conn = get_db()
    user = conn.execute("SELECT * FROM users WHERE username = ?", (username,)).fetchone()
    conn.close()
    return user

def get_user_by_id(user_id):
    conn = get_db()
    user = conn.execute("SELECT * FROM users WHERE id = ?", (user_id,)).fetchone()
    conn.close()
    return user

def get_all_users():
    conn = get_db()
    users = conn.execute("SELECT id, username, role, is_active, created_at, last_login FROM users").fetchall()
    conn.close()
    return users

def get_primers(filters=None, limit=20, offset=0):
    conn = get_db()
    query = "SELECT * FROM primers WHERE is_active = 1"
    params = []
    if filters:
        if filters.get('name'):
            query += " AND name LIKE ?"
            params.append(f"%{filters['name']}%")
        if filters.get('gene'):
            query += " AND gene LIKE ?"
            params.append(f"%{filters['gene']}%")
        if filters.get('organism'):
            query += " AND organism LIKE ?"
            params.append(f"%{filters['organism']}%")
    query += " ORDER BY name LIMIT ? OFFSET ?"
    params.extend([limit, offset])
    primers = conn.execute(query, params).fetchall()
    conn.close()
    return primers

def get_primer_by_id(primer_id):
    conn = get_db()
    primer = conn.execute("SELECT * FROM primers WHERE id = ?", (primer_id,)).fetchone()
    conn.close()
    return primer

def get_primer_by_name(name):
    conn = get_db()
    primer = conn.execute("SELECT * FROM primers WHERE name = ?", (name,)).fetchone()
    conn.close()
    return primer

def get_custom_fields(primer_id):
    conn = get_db()
    fields = conn.execute("SELECT * FROM custom_fields WHERE primer_id = ?", (primer_id,)).fetchall()
    conn.close()
    return fields

def get_probes_by_pair(pair_name):
    conn = get_db()
    probes = conn.execute("SELECT * FROM probes WHERE pair_name = ?", (pair_name,)).fetchall()
    conn.close()
    return probes

def get_pcr_programs_by_pair(pair_name):
    conn = get_db()
    programs = conn.execute("SELECT * FROM pcr_programs WHERE pair_name = ? ORDER BY is_default DESC, program_name", (pair_name,)).fetchall()
    conn.close()
    return programs

def get_pcr_program_by_id(program_id):
    conn = get_db()
    program = conn.execute("SELECT * FROM pcr_programs WHERE id = ?", (program_id,)).fetchone()
    conn.close()
    return program

def get_pcr_steps(program_id):
    conn = get_db()
    steps = conn.execute("SELECT * FROM pcr_steps WHERE program_id = ? ORDER BY step_order", (program_id,)).fetchall()
    conn.close()
    return steps

def get_panels():
    conn = get_db()
    panels = conn.execute("SELECT * FROM reaction_panels ORDER BY panel_name").fetchall()
    conn.close()
    return panels

def get_panel_by_id(panel_id):
    conn = get_db()
    panel = conn.execute("SELECT * FROM reaction_panels WHERE id = ?", (panel_id,)).fetchone()
    conn.close()
    return panel

def get_panel_primers(panel_id):
    conn = get_db()
    items = conn.execute('''
        SELECT pp.*, p.name as primer_name, pr.sequence as probe_sequence
        FROM panel_primers pp
        LEFT JOIN primers p ON pp.primer_id = p.id
        LEFT JOIN probes pr ON pp.probe_id = pr.id
        WHERE pp.panel_id = ?
    ''', (panel_id,)).fetchall()
    conn.close()
    return items

def get_reset_requests(status='pending'):
    conn = get_db()
    requests = conn.execute('''
        SELECT prr.*, u.username
        FROM password_reset_requests prr
        JOIN users u ON prr.user_id = u.id
        WHERE prr.status = ?
        ORDER BY prr.requested_at DESC
    ''', (status,)).fetchall()
    conn.close()
    return requests

def get_reset_request_by_id(request_id):
    conn = get_db()
    request = conn.execute('''
        SELECT prr.*, u.username
        FROM password_reset_requests prr
        JOIN users u ON prr.user_id = u.id
        WHERE prr.id = ?
    ''', (request_id,)).fetchone()
    conn.close()
    return request

def get_editing_lock(primer_id):
    conn = get_db()
    lock = conn.execute('''
        SELECT el.*, u.username
        FROM editing_locks el
        JOIN users u ON el.user_id = u.id
        WHERE el.primer_id = ? AND el.expires_at > CURRENT_TIMESTAMP
    ''', (primer_id,)).fetchone()
    conn.close()
    return lock

def log_audit(user_id, action, details=None):
    conn = get_db()
    conn.execute(
        "INSERT INTO audit_log (user_id, action, details) VALUES (?, ?, ?)",
        (user_id, action, details)
    )
    conn.commit()
    conn.close()

# ==================== Auth Functions ====================

def login_user(username, password):
    user = get_user_by_username(username)
    if not user or not check_password_hash(user['password_hash'], password):
        return False, "Invalid username or password."
    if not user['is_active']:
        return False, "Account is deactivated."
    session['user_id'] = user['id']
    session['username'] = user['username']
    session['role'] = user['role']
    log_audit(user['id'], 'login', f"User {username} logged in")
    return True, "Login successful."

def logout_user():
    if 'user_id' in session:
        log_audit(session['user_id'], 'logout', f"User {session['username']} logged out")
    session.clear()

def is_authenticated():
    return 'user_id' in session

def is_admin():
    return session.get('role') == 'admin'

def is_editor_or_admin():
    return session.get('role') in ['admin', 'editor']

def change_password(user_id, old_password, new_password):
    user = get_user_by_id(user_id)
    if not user:
        return False, "User not found."
    if not check_password_hash(user['password_hash'], old_password):
        return False, "Current password is incorrect."
    if len(new_password) < 6:
        return False, "Password must be at least 6 characters."
    new_hash = generate_password_hash(new_password)
    conn = get_db()
    conn.execute("UPDATE users SET password_hash = ? WHERE id = ?", (new_hash, user_id))
    conn.commit()
    conn.close()
    log_audit(user_id, 'password_change', "User changed password")
    return True, "Password changed successfully."

def admin_reset_password(user_id, new_password):
    if len(new_password) < 6:
        return False, "Password must be at least 6 characters."
    user = get_user_by_id(user_id)
    if not user:
        return False, "User not found."
    new_hash = generate_password_hash(new_password)
    conn = get_db()
    conn.execute("UPDATE users SET password_hash = ? WHERE id = ?", (new_hash, user_id))
    conn.commit()
    conn.close()
    log_audit(session.get('user_id'), 'admin_reset_password', f"Admin reset password for {user['username']}")
    return True, "Password reset successfully."

def create_reset_request(username):
    user = get_user_by_username(username)
    if not user:
        return False, "Username not found."
    conn = get_db()
    today = datetime.now().strftime('%Y-%m-%d')
    count = conn.execute('''
        SELECT COUNT(*) as count FROM password_reset_requests
        WHERE user_id = ? AND DATE(requested_at) = ? AND status = 'pending'
    ''', (user['id'], today)).fetchone()['count']
    if count >= 3:
        conn.close()
        return False, "Too many reset requests today."
    conn.execute("INSERT INTO password_reset_requests (user_id) VALUES (?)", (user['id'],))
    conn.commit()
    conn.close()
    log_audit(user['id'], 'reset_request', f"User {username} requested password reset")
    return True, "Request submitted. Admin will contact you."

def resolve_reset_request(request_id, new_password, admin_id):
    req = get_reset_request_by_id(request_id)
    if not req or req['status'] != 'pending':
        return False, "Invalid request."
    if len(new_password) < 6:
        return False, "Password must be at least 6 characters."
    new_hash = generate_password_hash(new_password)
    conn = get_db()
    conn.execute('''
        UPDATE password_reset_requests
        SET status = 'resolved', resolved_at = CURRENT_TIMESTAMP,
            resolved_by = ?, new_password_hash = ?
        WHERE id = ?
    ''', (admin_id, new_hash, request_id))
    conn.execute("UPDATE users SET password_hash = ? WHERE id = ?", (new_hash, req['user_id']))
    conn.commit()
    conn.close()
    log_audit(admin_id, 'reset_resolved', f"Resolved reset request for {req['username']}")
    return True, "Password reset successfully."

def cancel_reset_request(request_id, admin_id):
    req = get_reset_request_by_id(request_id)
    if not req:
        return False, "Request not found."
    if req['status'] != 'pending':
        return False, "Request already processed."
    conn = get_db()
    conn.execute(
        "UPDATE password_reset_requests SET status = 'canceled', resolved_by = ?, resolved_at = CURRENT_TIMESTAMP WHERE id = ?",
        (admin_id, request_id)
    )
    conn.commit()
    conn.close()
    log_audit(admin_id, 'reset_canceled', f"Canceled reset request for {req['username']}")
    return True, "Request canceled."

def create_user(username, password, role='viewer'):
    if len(username) < 3:
        return False, "Username must be at least 3 characters."
    if len(password) < 6:
        return False, "Password must be at least 6 characters."
    if get_user_by_username(username):
        return False, "Username already exists."
    conn = get_db()
    conn.execute(
        "INSERT INTO users (username, password_hash, role) VALUES (?, ?, ?)",
        (username, generate_password_hash(password), role)
    )
    conn.commit()
    conn.close()
    log_audit(session.get('user_id'), 'user_created', f"User {username} created")
    return True, "User created."

def update_user_role(user_id, new_role):
    user = get_user_by_id(user_id)
    if not user or user['username'] == 'admin':
        return False, "Cannot change admin role."
    conn = get_db()
    conn.execute("UPDATE users SET role = ? WHERE id = ?", (new_role, user_id))
    conn.commit()
    conn.close()
    log_audit(session.get('user_id'), 'user_role_changed', f"User {user['username']} role changed to {new_role}")
    return True, "Role updated."

def deactivate_user(user_id):
    user = get_user_by_id(user_id)
    if not user or user['username'] == 'admin':
        return False, "Cannot deactivate admin."
    conn = get_db()
    conn.execute("UPDATE users SET is_active = 0 WHERE id = ?", (user_id,))
    conn.commit()
    conn.close()
    log_audit(session.get('user_id'), 'user_deactivated', f"User {user['username']} deactivated")
    return True, "User deactivated."

def activate_user(user_id):
    user = get_user_by_id(user_id)
    if not user:
        return False, "User not found."
    conn = get_db()
    conn.execute("UPDATE users SET is_active = 1 WHERE id = ?", (user_id,))
    conn.commit()
    conn.close()
    log_audit(session.get('user_id'), 'user_activated', f"User {user['username']} activated")
    return True, "User activated."

def is_valid_sequence(seq):
    valid_chars = set('ATCGatcgRYWSMKBDHVNrywsmkbdhvnI')
    return all(c in valid_chars for c in seq)

# ==================== Decorators ====================

def login_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if not is_authenticated():
            flash('Please log in first.', 'warning')
            return redirect(url_for('login'))
        return f(*args, **kwargs)
    return decorated

def admin_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if not is_authenticated():
            flash('Please log in first.', 'warning')
            return redirect(url_for('login'))
        if not is_admin():
            flash('Admin access required.', 'danger')
            return redirect(url_for('dashboard'))
        return f(*args, **kwargs)
    return decorated

def editor_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if not is_authenticated():
            flash('Please log in first.', 'warning')
            return redirect(url_for('login'))
        if not is_editor_or_admin():
            flash('Edit access required.', 'danger')
            return redirect(url_for('dashboard'))
        return f(*args, **kwargs)
    return decorated

# ==================== Routes ====================

@app.route('/')
@login_required
def dashboard():
    user = get_user_by_id(session['user_id'])
    primers = get_primers(limit=1000)
    users = get_all_users()
    panels = get_panels()
    return render_template('dashboard.html',
                           user=user,
                           primers=primers,
                           primers_count=len(primers),
                           users_count=len(users),
                           panels=panels)

@app.route('/login', methods=['GET', 'POST'])
def login():
    if is_authenticated():
        return redirect(url_for('dashboard'))
    if request.method == 'POST':
        username = request.form.get('username')
        password = request.form.get('password')
        success, msg = login_user(username, password)
        flash(msg, 'success' if success else 'danger')
        if success:
            return redirect(url_for('dashboard'))
    return render_template('login.html')

@app.route('/logout')
@login_required
def logout():
    logout_user()
    flash('Logged out.', 'success')
    return redirect(url_for('login'))

@app.route('/profile', methods=['GET', 'POST'])
@login_required
def profile():
    user = get_user_by_id(session['user_id'])
    if request.method == 'POST':
        old = request.form.get('old_password')
        new = request.form.get('new_password')
        confirm = request.form.get('confirm_password')
        if new != confirm:
            flash('Passwords do not match.', 'danger')
        else:
            success, msg = change_password(user['id'], old, new)
            flash(msg, 'success' if success else 'danger')
    return render_template('profile.html', user=user)

@app.route('/forgot-password', methods=['GET', 'POST'])
def forgot_password():
    if request.method == 'POST':
        username = request.form.get('username')
        success, msg = create_reset_request(username)
        flash(msg, 'success' if success else 'danger')
        if success:
            return render_template('forgot_password_done.html')
    return render_template('forgot_password.html')

@app.route('/primers')
@login_required
def primer_list():
    page = request.args.get('page', 1, type=int)
    limit = 20
    offset = (page - 1) * limit
    filters = {}
    for key in ['name', 'gene', 'organism']:
        if request.args.get(key):
            filters[key] = request.args.get(key)
    primers = get_primers(filters, limit, offset)
    return render_template('primer_list.html', primers=primers, page=page, filters=filters)

@app.route('/primers/add', methods=['GET', 'POST'])
@editor_required
def primer_add():
    if request.method == 'POST':
        data = {
            'name': request.form.get('name'),
            'forward_sequence': request.form.get('forward_sequence'),
            'reverse_sequence': request.form.get('reverse_sequence'),
            'pair_name': request.form.get('pair_name'),
            'gene': request.form.get('gene'),
            'organism': request.form.get('organism'),
            'strain_or_serotype': request.form.get('strain_or_serotype'),
            'pcr_type': request.form.get('pcr_type'),
            'amplicon_length': request.form.get('amplicon_length', type=int),
            'estimated_tm': request.form.get('estimated_tm', type=float),
            'experimental_tm': request.form.get('experimental_tm', type=float),
            'reference': request.form.get('reference'),
            'diagnostic_limitations': request.form.get('diagnostic_limitations'),
            'for_sequencing': request.form.get('for_sequencing') == 'on',
            'for_diagnosis': request.form.get('for_diagnosis') == 'on',
            'binding_region': request.form.get('binding_region'),
            'binding_detail': request.form.get('binding_detail'),
            'stock_concentration': request.form.get('stock_concentration', type=float),
            'working_volume_per_reaction': request.form.get('working_volume_per_reaction', type=float),
            'general_notes': request.form.get('general_notes'),
        }
        if not data['name'] or not data['forward_sequence']:
            flash('Name and Forward Sequence are required.', 'danger')
            return render_template('primer_add.html')
        if not is_valid_sequence(data['forward_sequence']):
            flash('Forward Sequence contains invalid characters.', 'danger')
            return render_template('primer_add.html')
        if get_primer_by_name(data['name']):
            flash('Primer with this name already exists.', 'danger')
            return render_template('primer_add.html')

        # محاسبه خودکار Tm اگر کاربر خالی گذاشته بود
        if not data['estimated_tm'] and data['forward_sequence']:
            data['estimated_tm'] = calculate_tm(data['forward_sequence'])

        conn = get_db()
        c = conn.cursor()
        c.execute('''
            INSERT INTO primers (
                name, forward_sequence, reverse_sequence, pair_name, gene, organism,
                strain_or_serotype, pcr_type, amplicon_length, estimated_tm, experimental_tm,
                reference, diagnostic_limitations, for_sequencing, for_diagnosis,
                binding_region, binding_detail, stock_concentration,
                working_volume_per_reaction, general_notes, added_by
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            data['name'], data['forward_sequence'], data['reverse_sequence'],
            data['pair_name'], data['gene'], data['organism'],
            data['strain_or_serotype'], data['pcr_type'], data['amplicon_length'],
            data['estimated_tm'], data['experimental_tm'], data['reference'],
            data['diagnostic_limitations'], data['for_sequencing'], data['for_diagnosis'],
            data['binding_region'], data['binding_detail'], data['stock_concentration'],
            data['working_volume_per_reaction'], data['general_notes'],
            session['user_id']
        ))
        primer_id = c.lastrowid

        custom_names = request.form.getlist('custom_field_name')
        custom_values = request.form.getlist('custom_field_value')
        for name, value in zip(custom_names, custom_values):
            if name and value:
                c.execute(
                    "INSERT INTO custom_fields (primer_id, field_name, field_value) VALUES (?, ?, ?)",
                    (primer_id, name, value)
                )
        conn.commit()
        conn.close()

        log_audit(session['user_id'], 'primer_added', f"Added primer {data['name']}")
        flash('Primer added successfully.', 'success')
        return redirect(url_for('primer_detail', primer_id=primer_id))
    return render_template('primer_add.html')

@app.route('/primers/<int:primer_id>')
@login_required
def primer_detail(primer_id):
    primer = get_primer_by_id(primer_id)
    if not primer:
        flash('Primer not found.', 'danger')
        return redirect(url_for('primer_list'))
    custom_fields = get_custom_fields(primer_id)
    probes = get_probes_by_pair(primer['pair_name']) if primer['pair_name'] else []
    programs = get_pcr_programs_by_pair(primer['pair_name']) if primer['pair_name'] else []
    lock = get_editing_lock(primer_id)

    programs_with_steps = []
    for prog in programs:
        steps = get_pcr_steps(prog['id'])
        programs_with_steps.append({'program': prog, 'steps': steps})

    return render_template('primer_detail.html', primer=primer, custom_fields=custom_fields,
                           probes=probes, programs=programs_with_steps, lock=lock)

@app.route('/primers/<int:primer_id>/edit', methods=['GET', 'POST'])
@editor_required
def primer_edit(primer_id):
    primer = get_primer_by_id(primer_id)
    if not primer:
        flash('Primer not found.', 'danger')
        return redirect(url_for('primer_list'))

    if request.method == 'GET':
        lock = get_editing_lock(primer_id)
        if lock and lock['user_id'] != session['user_id']:
            flash(f"Being edited by {lock['username']}.", 'warning')
            return redirect(url_for('primer_detail', primer_id=primer_id))
        conn = get_db()
        expires = datetime.now() + timedelta(minutes=15)
        conn.execute(
            "INSERT OR REPLACE INTO editing_locks (primer_id, user_id, expires_at) VALUES (?, ?, ?)",
            (primer_id, session['user_id'], expires)
        )
        conn.commit()
        conn.close()

    if request.method == 'POST':
        data = {
            'name': request.form.get('name'),
            'forward_sequence': request.form.get('forward_sequence'),
            'reverse_sequence': request.form.get('reverse_sequence'),
            'pair_name': request.form.get('pair_name'),
            'gene': request.form.get('gene'),
            'organism': request.form.get('organism'),
            'strain_or_serotype': request.form.get('strain_or_serotype'),
            'pcr_type': request.form.get('pcr_type'),
            'amplicon_length': request.form.get('amplicon_length', type=int),
            'estimated_tm': request.form.get('estimated_tm', type=float),
            'experimental_tm': request.form.get('experimental_tm', type=float),
            'reference': request.form.get('reference'),
            'diagnostic_limitations': request.form.get('diagnostic_limitations'),
            'for_sequencing': request.form.get('for_sequencing') == 'on',
            'for_diagnosis': request.form.get('for_diagnosis') == 'on',
            'binding_region': request.form.get('binding_region'),
            'binding_detail': request.form.get('binding_detail'),
            'stock_concentration': request.form.get('stock_concentration', type=float),
            'working_volume_per_reaction': request.form.get('working_volume_per_reaction', type=float),
            'general_notes': request.form.get('general_notes'),
        }
        if not data['name'] or not data['forward_sequence']:
            flash('Name and Forward Sequence are required.', 'danger')
            return render_template('primer_edit.html', primer=primer)
        if not is_valid_sequence(data['forward_sequence']):
            flash('Forward Sequence contains invalid characters.', 'danger')
            return render_template('primer_edit.html', primer=primer)

        # محاسبه خودکار Tm اگر کاربر خالی گذاشته بود
        if not data['estimated_tm'] and data['forward_sequence']:
            data['estimated_tm'] = calculate_tm(data['forward_sequence'])

        conn = get_db()
        c = conn.cursor()
        c.execute('''
            UPDATE primers SET
                name = ?, forward_sequence = ?, reverse_sequence = ?, pair_name = ?,
                gene = ?, organism = ?, strain_or_serotype = ?, pcr_type = ?,
                amplicon_length = ?, estimated_tm = ?, experimental_tm = ?,
                reference = ?, diagnostic_limitations = ?,
                for_sequencing = ?, for_diagnosis = ?,
                binding_region = ?, binding_detail = ?,
                stock_concentration = ?, working_volume_per_reaction = ?,
                general_notes = ?
            WHERE id = ?
        ''', (
            data['name'], data['forward_sequence'], data['reverse_sequence'],
            data['pair_name'], data['gene'], data['organism'],
            data['strain_or_serotype'], data['pcr_type'], data['amplicon_length'],
            data['estimated_tm'], data['experimental_tm'], data['reference'],
            data['diagnostic_limitations'], data['for_sequencing'], data['for_diagnosis'],
            data['binding_region'], data['binding_detail'], data['stock_concentration'],
            data['working_volume_per_reaction'], data['general_notes'],
            primer_id
        ))
        c.execute("DELETE FROM custom_fields WHERE primer_id = ?", (primer_id,))
        custom_names = request.form.getlist('custom_field_name')
        custom_values = request.form.getlist('custom_field_value')
        for name, value in zip(custom_names, custom_values):
            if name and value:
                c.execute(
                    "INSERT INTO custom_fields (primer_id, field_name, field_value) VALUES (?, ?, ?)",
                    (primer_id, name, value)
                )
        conn.commit()
        conn.close()

        conn = get_db()
        conn.execute("DELETE FROM editing_locks WHERE primer_id = ?", (primer_id,))
        conn.commit()
        conn.close()

        log_audit(session['user_id'], 'primer_edited', f"Edited primer {data['name']}")
        flash('Primer updated.', 'success')
        return redirect(url_for('primer_detail', primer_id=primer_id))

    custom_fields = get_custom_fields(primer_id)
    return render_template('primer_edit.html', primer=primer, custom_fields=custom_fields)

@app.route('/primers/<int:primer_id>/unlock')
@admin_required
def primer_unlock(primer_id):
    conn = get_db()
    conn.execute("DELETE FROM editing_locks WHERE primer_id = ?", (primer_id,))
    conn.commit()
    conn.close()
    flash('Lock removed.', 'success')
    return redirect(url_for('primer_detail', primer_id=primer_id))

@app.route('/primers/<int:primer_id>/delete')
@admin_required
def primer_delete(primer_id):
    primer = get_primer_by_id(primer_id)
    if not primer:
        flash('Primer not found.', 'danger')
        return redirect(url_for('primer_list'))
    conn = get_db()
    conn.execute("UPDATE primers SET is_active = 0 WHERE id = ?", (primer_id,))
    conn.commit()
    conn.close()
    log_audit(session['user_id'], 'primer_deleted', f"Deleted primer {primer['name']}")
    flash('Primer deleted.', 'success')
    return redirect(url_for('primer_list'))

@app.route('/programs/<int:program_id>/edit', methods=['GET', 'POST'])
@editor_required
def pcr_program_edit(program_id):
    program = get_pcr_program_by_id(program_id)
    if not program:
        flash('Program not found.', 'danger')
        return redirect(url_for('primer_list'))

    primer = get_primer_by_name(program['pair_name'])
    if not primer:
        flash('Associated primer not found.', 'danger')
        return redirect(url_for('primer_list'))
    primer_id = primer['id']

    if request.method == 'POST':
        program_name = request.form.get('program_name')
        program_type = request.form.get('program_type')
        reaction_volume = request.form.get('reaction_volume', type=float)
        master_mix = request.form.get('master_mix')
        notes = request.form.get('notes')

        conn = get_db()
        conn.execute('''
            UPDATE pcr_programs SET
                program_name = ?, program_type = ?, reaction_volume = ?,
                master_mix = ?, notes = ?, updated_at = CURRENT_TIMESTAMP
            WHERE id = ?
        ''', (program_name, program_type, reaction_volume, master_mix, notes, program_id))
        conn.execute("DELETE FROM pcr_steps WHERE program_id = ?", (program_id,))

        step_types = request.form.getlist('step_type[]')
        temperatures = request.form.getlist('temperature[]')
        durations = request.form.getlist('duration_sec[]')
        cycle_repeats = request.form.getlist('cycle_repeat[]')
        is_reads = request.form.getlist('is_read_step[]')

        for i, (stype, temp, dur, cycle, is_read) in enumerate(zip(step_types, temperatures, durations, cycle_repeats, is_reads)):
            if not stype:
                continue
            conn.execute('''
                INSERT INTO pcr_steps (
                    program_id, step_order, step_type, temperature,
                    duration_sec, cycle_repeat, is_read_step
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', (program_id, i+1, stype, float(temp) if temp else None,
                  int(dur) if dur else None, int(cycle) if cycle else None,
                  1 if is_read else 0))
        conn.commit()
        conn.close()

        flash('Program updated.', 'success')
        return redirect(url_for('primer_detail', primer_id=primer_id))

    steps = get_pcr_steps(program_id)
    return render_template('pcr_program_edit.html', program=program, steps=steps, primer_id=primer_id)

@app.route('/programs/add', methods=['POST'])
@editor_required
def pcr_program_add():
    primer_id = request.form.get('primer_id')
    program_name = request.form.get('program_name')
    is_default = request.form.get('is_default') == 'on'

    if not primer_id or not program_name:
        flash('Primer ID and Program Name are required.', 'danger')
        return redirect(url_for('primer_list'))

    primer = get_primer_by_id(primer_id)
    if not primer:
        flash('Primer not found.', 'danger')
        return redirect(url_for('primer_list'))

    conn = get_db()
    if is_default:
        conn.execute("UPDATE pcr_programs SET is_default = 0 WHERE pair_name = ?", (primer['pair_name'],))
    conn.execute(
        "INSERT INTO pcr_programs (pair_name, program_name, is_default) VALUES (?, ?, ?)",
        (primer['pair_name'], program_name, is_default)
    )
    conn.commit()
    conn.close()

    log_audit(session['user_id'], 'pcr_program_added', f"Added program {program_name} for {primer['name']}")
    flash('Program created.', 'success')
    return redirect(url_for('primer_detail', primer_id=primer_id))

@app.route('/programs/<int:program_id>/set-default')
@editor_required
def pcr_program_set_default(program_id):
    program = get_pcr_program_by_id(program_id)
    if not program:
        flash('Program not found.', 'danger')
        return redirect(url_for('primer_list'))

    primer = get_primer_by_name(program['pair_name'])
    if not primer:
        flash('Associated primer not found.', 'danger')
        return redirect(url_for('primer_list'))
    primer_id = primer['id']

    conn = get_db()
    conn.execute("UPDATE pcr_programs SET is_default = 0 WHERE pair_name = ?", (program['pair_name'],))
    conn.execute("UPDATE pcr_programs SET is_default = 1 WHERE id = ?", (program_id,))
    conn.commit()
    conn.close()

    flash('Default program updated.', 'success')
    return redirect(url_for('primer_detail', primer_id=primer_id))

@app.route('/programs/<int:program_id>/delete')
@editor_required
def pcr_program_delete(program_id):
    program = get_pcr_program_by_id(program_id)
    if not program:
        flash('Program not found.', 'danger')
        return redirect(url_for('primer_list'))

    primer = get_primer_by_name(program['pair_name'])
    if not primer:
        flash('Associated primer not found.', 'danger')
        return redirect(url_for('primer_list'))
    primer_id = primer['id']

    conn = get_db()
    conn.execute("DELETE FROM pcr_programs WHERE id = ?", (program_id,))
    conn.commit()
    conn.close()

    flash('Program deleted.', 'success')
    return redirect(url_for('primer_detail', primer_id=primer_id))

@app.route('/probes/add', methods=['POST'])
@editor_required
def probe_add():
    pair_name = request.form.get('pair_name')
    sequence = request.form.get('sequence')
    probe_type = request.form.get('probe_type')
    reporter = request.form.get('reporter')
    quencher = request.form.get('quencher')
    modifications = request.form.get('modifications')
    notes = request.form.get('probe_notes')

    if not sequence or not is_valid_sequence(sequence):
        flash('Invalid probe sequence.', 'danger')
        return redirect(url_for('primer_detail', primer_id=pair_name))

    conn = get_db()
    conn.execute('''
        INSERT INTO probes (pair_name, sequence, probe_type, reporter, quencher, modifications, notes)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', (pair_name, sequence, probe_type, reporter, quencher, modifications, notes))
    conn.commit()
    conn.close()

    flash('Probe added.', 'success')
    return redirect(url_for('primer_detail', primer_id=pair_name))

@app.route('/probes/<int:probe_id>/delete')
@editor_required
def probe_delete(probe_id):
    conn = get_db()
    probe = conn.execute("SELECT pair_name FROM probes WHERE id = ?", (probe_id,)).fetchone()
    if probe:
        conn.execute("DELETE FROM probes WHERE id = ?", (probe_id,))
        conn.commit()
        flash('Probe deleted.', 'success')
    conn.close()
    return redirect(url_for('primer_detail', primer_id=probe['pair_name']))

@app.route('/panels')
@login_required
def panel_list():
    panels = get_panels()
    return render_template('panel_list.html', panels=panels)

@app.route('/panels/add', methods=['GET', 'POST'])
@editor_required
def panel_add():
    if request.method == 'POST':
        panel_name = request.form.get('panel_name')
        description = request.form.get('description')
        organism = request.form.get('organism')

        conn = get_db()
        c = conn.cursor()
        c.execute(
            "INSERT INTO reaction_panels (panel_name, description, organism, created_by) VALUES (?, ?, ?, ?)",
            (panel_name, description, organism, session['user_id'])
        )
        panel_id = c.lastrowid

        primer_ids = request.form.getlist('primer_ids')
        for pid in primer_ids:
            if pid:
                c.execute(
                    "INSERT INTO panel_primers (panel_id, primer_id) VALUES (?, ?)",
                    (panel_id, pid)
                )
        conn.commit()
        conn.close()

        flash('Panel created.', 'success')
        return redirect(url_for('panel_detail', panel_id=panel_id))

    primers = get_primers(limit=1000)
    return render_template('panel_add.html', primers=primers)

@app.route('/panels/<int:panel_id>')
@login_required
def panel_detail(panel_id):
    panel = get_panel_by_id(panel_id)
    if not panel:
        flash('Panel not found.', 'danger')
        return redirect(url_for('panel_list'))
    items = get_panel_primers(panel_id)
    return render_template('panel_detail.html', panel=panel, items=items)

@app.route('/panels/<int:panel_id>/delete')
@admin_required
def panel_delete(panel_id):
    conn = get_db()
    conn.execute("DELETE FROM reaction_panels WHERE id = ?", (panel_id,))
    conn.commit()
    conn.close()
    flash('Panel deleted.', 'success')
    return redirect(url_for('panel_list'))

@app.route('/admin/users')
@admin_required
def admin_users():
    users = get_all_users()
    return render_template('admin_users.html', users=users)

@app.route('/admin/users/add', methods=['POST'])
@admin_required
def admin_user_add():
    username = request.form.get('username')
    password = request.form.get('password')
    role = request.form.get('role', 'viewer')
    success, msg = create_user(username, password, role)
    flash(msg, 'success' if success else 'danger')
    return redirect(url_for('admin_users'))

@app.route('/admin/users/<int:user_id>/role', methods=['POST'])
@admin_required
def admin_user_role(user_id):
    role = request.form.get('role')
    success, msg = update_user_role(user_id, role)
    flash(msg, 'success' if success else 'danger')
    return redirect(url_for('admin_users'))

@app.route('/admin/users/<int:user_id>/toggle')
@admin_required
def admin_user_toggle(user_id):
    user = get_user_by_id(user_id)
    if user and user['is_active']:
        success, msg = deactivate_user(user_id)
    else:
        success, msg = activate_user(user_id)
    flash(msg, 'success' if success else 'danger')
    return redirect(url_for('admin_users'))

@app.route('/admin/users/<int:user_id>/reset-password', methods=['POST'])
@admin_required
def admin_user_reset_password(user_id):
    new_password = request.form.get('new_password')
    success, msg = admin_reset_password(user_id, new_password)
    flash(msg, 'success' if success else 'danger')
    return redirect(url_for('admin_users'))

@app.route('/admin/reset-requests')
@admin_required
def admin_reset_requests():
    requests = get_reset_requests('pending')
    return render_template('admin_reset_requests.html', requests=requests)

@app.route('/admin/reset-requests/<int:request_id>')
@admin_required
def admin_reset_detail(request_id):
    req = get_reset_request_by_id(request_id)
    if not req:
        flash('Request not found.', 'danger')
        return redirect(url_for('admin_reset_requests'))
    return render_template('admin_reset_detail.html', request=req)

@app.route('/admin/reset-requests/<int:request_id>/resolve', methods=['POST'])
@admin_required
def admin_reset_resolve(request_id):
    new_password = request.form.get('new_password')
    success, msg = resolve_reset_request(request_id, new_password, session['user_id'])
    flash(msg, 'success' if success else 'danger')
    return redirect(url_for('admin_reset_requests'))

@app.route('/admin/reset-requests/<int:request_id>/cancel')
@admin_required
def admin_reset_cancel(request_id):
    conn = get_db()
    conn.execute(
        "UPDATE password_reset_requests SET status = 'canceled', resolved_by = ?, resolved_at = CURRENT_TIMESTAMP WHERE id = ?",
        (session['user_id'], request_id)
    )
    conn.commit()
    conn.close()
    flash('Request canceled.', 'success')
    return redirect(url_for('admin_reset_requests'))

@app.route('/backup')
@admin_required
def backup_page():
    backups = []
    if os.path.exists(BACKUP_DIR):
        for f in sorted(os.listdir(BACKUP_DIR), reverse=True):
            if f.endswith('.db'):
                path = os.path.join(BACKUP_DIR, f)
                size = os.path.getsize(path)
                modified = datetime.fromtimestamp(os.path.getmtime(path))
                backups.append({'name': f, 'size': size, 'modified': modified})
    return render_template('backup.html', backups=backups)

@app.route('/backup/create')
@admin_required
def backup_create():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_name = f"backup_{timestamp}.db"
    backup_path = os.path.join(BACKUP_DIR, backup_name)
    shutil.copy2(DATABASE_PATH, backup_path)

    cutoff = datetime.now() - timedelta(days=30)
    for f in os.listdir(BACKUP_DIR):
        if f.endswith('.db'):
            path = os.path.join(BACKUP_DIR, f)
            mtime = datetime.fromtimestamp(os.path.getmtime(path))
            if mtime < cutoff:
                os.remove(path)

    log_audit(session['user_id'], 'backup_created', f"Backup {backup_name} created")
    flash(f'Backup {backup_name} created.', 'success')
    return redirect(url_for('backup_page'))

@app.route('/backup/download/<filename>')
@admin_required
def backup_download(filename):
    path = os.path.join(BACKUP_DIR, filename)
    if not os.path.exists(path):
        flash('File not found.', 'danger')
        return redirect(url_for('backup_page'))
    return send_file(path, as_attachment=True, download_name=filename)

@app.route('/backup/restore/<filename>')
@admin_required
def backup_restore(filename):
    path = os.path.join(BACKUP_DIR, filename)
    if not os.path.exists(path):
        flash('File not found.', 'danger')
        return redirect(url_for('backup_page'))

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    pre_restore = f"pre_restore_{timestamp}.db"
    shutil.copy2(DATABASE_PATH, os.path.join(BACKUP_DIR, pre_restore))
    shutil.copy2(path, DATABASE_PATH)

    log_audit(session['user_id'], 'backup_restored', f"Restored backup {filename}")
    flash('Database restored.', 'success')
    return redirect(url_for('backup_page'))

@app.route('/search', methods=['GET', 'POST'])
@login_required
def search():
    results = []
    if request.method == 'POST':
        filters = {}
        for field in ['name', 'pair_name', 'gene', 'organism', 'strain_or_serotype',
                      'pcr_type', 'reference', 'diagnostic_limitations', 'binding_region']:
            val = request.form.get(field)
            if val:
                filters[field] = val

        if request.form.get('amplicon_length_min'):
            filters['amplicon_length_min'] = int(request.form.get('amplicon_length_min'))
        if request.form.get('amplicon_length_max'):
            filters['amplicon_length_max'] = int(request.form.get('amplicon_length_max'))
        if request.form.get('estimated_tm_min'):
            filters['estimated_tm_min'] = float(request.form.get('estimated_tm_min'))
        if request.form.get('estimated_tm_max'):
            filters['estimated_tm_max'] = float(request.form.get('estimated_tm_max'))
        if request.form.get('for_sequencing') == 'on':
            filters['for_sequencing'] = 1
        if request.form.get('for_diagnosis') == 'on':
            filters['for_diagnosis'] = 1

        query = "SELECT * FROM primers WHERE is_active = 1"
        params = []
        for key, val in filters.items():
            if key.endswith('_min'):
                query += f" AND {key[:-4]} >= ?"
                params.append(val)
            elif key.endswith('_max'):
                query += f" AND {key[:-4]} <= ?"
                params.append(val)
            elif key in ['for_sequencing', 'for_diagnosis']:
                query += f" AND {key} = ?"
                params.append(val)
            else:
                query += f" AND {key} LIKE ?"
                params.append(f"%{val}%")

        conn = get_db()
        results = conn.execute(query, params).fetchall()
        conn.close()

    return render_template('search.html', results=results)

@app.route('/api/primers/export')
@login_required
def api_export_primers():
    primer_ids = request.args.getlist('ids')
    if not primer_ids:
        return jsonify({'error': 'No primers selected'}), 400

    conn = get_db()
    placeholders = ','.join('?' * len(primer_ids))
    primers = conn.execute(f'''
        SELECT name, forward_sequence, reverse_sequence, pair_name
        FROM primers WHERE id IN ({placeholders}) AND is_active = 1
    ''', primer_ids).fetchall()
    conn.close()

    result = [{'name': p['name'], 'forward_sequence': p['forward_sequence'],
               'reverse_sequence': p['reverse_sequence'], 'pair_name': p['pair_name']} for p in primers]
    return jsonify(result)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5001)
