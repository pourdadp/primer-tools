import sqlite3
from datetime import datetime
from werkzeug.security import generate_password_hash
from config import DATABASE_PATH

def get_db():
    conn = sqlite3.connect(DATABASE_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    conn = get_db()
    c = conn.cursor()

    # Users table
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

    # Primers table
    c.execute('''
        CREATE TABLE IF NOT EXISTS primers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT UNIQUE NOT NULL,
            forward_sequence TEXT NOT NULL,
            reverse_sequence TEXT,
            pair_name TEXT,
            gene TEXT,
            organism TEXT,
            strain_or_serotype TEXT,
            pcr_type TEXT,
            amplicon_length INTEGER,
            estimated_tm REAL,
            experimental_tm REAL,
            reference TEXT,
            diagnostic_limitations TEXT,
            for_sequencing BOOLEAN DEFAULT 0,
            for_diagnosis BOOLEAN DEFAULT 0,
            binding_region TEXT,
            binding_detail TEXT,
            stock_concentration REAL,
            working_volume_per_reaction REAL,
            general_notes TEXT,
            is_active BOOLEAN DEFAULT 1,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            added_by INTEGER,
            FOREIGN KEY (added_by) REFERENCES users(id)
        )
    ''')

    # Probes table
    c.execute('''
        CREATE TABLE IF NOT EXISTS probes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pair_name TEXT NOT NULL,
            sequence TEXT NOT NULL,
            probe_type TEXT,
            reporter TEXT,
            quencher TEXT,
            modifications TEXT,
            notes TEXT,
            FOREIGN KEY (pair_name) REFERENCES primers(pair_name) ON DELETE CASCADE
        )
    ''')

    # Custom fields table
    c.execute('''
        CREATE TABLE IF NOT EXISTS custom_fields (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            primer_id INTEGER NOT NULL,
            field_name TEXT NOT NULL,
            field_value TEXT,
            FOREIGN KEY (primer_id) REFERENCES primers(id) ON DELETE CASCADE
        )
    ''')

    # PCR programs table
    c.execute('''
        CREATE TABLE IF NOT EXISTS pcr_programs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pair_name TEXT NOT NULL,
            program_name TEXT NOT NULL,
            is_default BOOLEAN DEFAULT 0,
            program_type TEXT,
            reaction_volume REAL,
            master_mix TEXT,
            notes TEXT,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (pair_name) REFERENCES primers(pair_name) ON DELETE CASCADE
        )
    ''')

    # PCR cycle groups table (NEW)
    c.execute('''
        CREATE TABLE IF NOT EXISTS pcr_cycle_groups (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            program_id INTEGER NOT NULL,
            cycle_number INTEGER NOT NULL,
            repeat_count INTEGER DEFAULT 1,
            FOREIGN KEY (program_id) REFERENCES pcr_programs(id) ON DELETE CASCADE
        )
    ''')

    # PCR steps table (with cycle_group_id)
    c.execute('''
        CREATE TABLE IF NOT EXISTS pcr_steps (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            program_id INTEGER NOT NULL,
            step_order INTEGER NOT NULL,
            step_type TEXT NOT NULL,
            temperature REAL,
            duration_sec INTEGER,
            is_read_step BOOLEAN DEFAULT 0,
            cycle_group_id INTEGER,
            FOREIGN KEY (program_id) REFERENCES pcr_programs(id) ON DELETE CASCADE,
            FOREIGN KEY (cycle_group_id) REFERENCES pcr_cycle_groups(id) ON DELETE SET NULL
        )
    ''')

    # Reaction panels table
    c.execute('''
        CREATE TABLE IF NOT EXISTS reaction_panels (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            panel_name TEXT UNIQUE NOT NULL,
            description TEXT,
            organism TEXT,
            created_by INTEGER,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (created_by) REFERENCES users(id)
        )
    ''')

    # Panel primers table
    c.execute('''
        CREATE TABLE IF NOT EXISTS panel_primers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            panel_id INTEGER NOT NULL,
            primer_id INTEGER NOT NULL,
            probe_id INTEGER,
            working_volume_per_reaction REAL,
            final_concentration REAL,
            FOREIGN KEY (panel_id) REFERENCES reaction_panels(id) ON DELETE CASCADE,
            FOREIGN KEY (primer_id) REFERENCES primers(id) ON DELETE CASCADE,
            FOREIGN KEY (probe_id) REFERENCES probes(id) ON DELETE SET NULL
        )
    ''')

    # Password reset requests table
    c.execute('''
        CREATE TABLE IF NOT EXISTS password_reset_requests (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER NOT NULL,
            requested_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            resolved_at DATETIME,
            status TEXT DEFAULT 'pending',
            resolved_by INTEGER,
            new_password_hash TEXT,
            admin_notes TEXT,
            FOREIGN KEY (user_id) REFERENCES users(id),
            FOREIGN KEY (resolved_by) REFERENCES users(id)
        )
    ''')

    # Editing locks table
    c.execute('''
        CREATE TABLE IF NOT EXISTS editing_locks (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            primer_id INTEGER NOT NULL,
            user_id INTEGER NOT NULL,
            locked_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            expires_at DATETIME,
            FOREIGN KEY (primer_id) REFERENCES primers(id) ON DELETE CASCADE,
            FOREIGN KEY (user_id) REFERENCES users(id)
        )
    ''')

    # Audit log table
    c.execute('''
        CREATE TABLE IF NOT EXISTS audit_log (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            user_id INTEGER,
            action TEXT NOT NULL,
            details TEXT,
            ip_address TEXT,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (user_id) REFERENCES users(id)
        )
    ''')

    # Primer aliases table
    c.execute('''
        CREATE TABLE IF NOT EXISTS primer_aliases (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            primer_id INTEGER NOT NULL,
            alias_name TEXT NOT NULL,
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (primer_id) REFERENCES primers(id) ON DELETE CASCADE,
            UNIQUE(primer_id, alias_name)
        )
    ''')

    # Create default admin user
    admin = c.execute("SELECT * FROM users WHERE username = 'admin'").fetchone()
    if not admin:
        admin_hash = generate_password_hash('admin123')
        c.execute(
            "INSERT INTO users (username, password_hash, role) VALUES (?, ?, ?)",
            ('admin', admin_hash, 'admin')
        )

    # Trigger for updated_at
    c.execute('''
        CREATE TRIGGER IF NOT EXISTS update_primer_timestamp
        AFTER UPDATE ON primers
        BEGIN
            UPDATE primers SET updated_at = CURRENT_TIMESTAMP WHERE id = NEW.id;
        END
    ''')

    conn.commit()
    conn.close()

# ==================== Helper Functions ====================

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

def get_primer_by_pair_name(pair_name):
    if not pair_name:
        return None
    conn = get_db()
    primer = conn.execute("SELECT * FROM primers WHERE pair_name = ? LIMIT 1", (pair_name,)).fetchone()
    conn.close()
    return primer

def get_aliases(primer_id):
    conn = get_db()
    aliases = conn.execute("SELECT * FROM primer_aliases WHERE primer_id = ?", (primer_id,)).fetchall()
    conn.close()
    return aliases

def add_alias_to_primer(primer_id, alias_name):
    conn = get_db()
    try:
        conn.execute(
            "INSERT INTO primer_aliases (primer_id, alias_name) VALUES (?, ?)",
            (primer_id, alias_name.strip())
        )
        conn.commit()
        conn.close()
        return True, "Alias added successfully."
    except sqlite3.IntegrityError:
        conn.close()
        return False, "This alias already exists for this primer."

def check_duplicate_sequence(forward_seq, reverse_seq=None, exclude_id=None):
    conn = get_db()
    query = "SELECT id, name, forward_sequence, reverse_sequence FROM primers WHERE is_active = 1"
    params = []
    if exclude_id:
        query += " AND id != ?"
        params.append(exclude_id)
    primers = conn.execute(query, params).fetchall()
    conn.close()

    for p in primers:
        if p['forward_sequence'] == forward_seq:
            return True, p['id'], p['name'], "forward"
        if reverse_seq and p['reverse_sequence'] == reverse_seq:
            return True, p['id'], p['name'], "reverse"
    return False, None, None, None

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

def get_pcr_cycle_groups(program_id):
    conn = get_db()
    groups = conn.execute("SELECT * FROM pcr_cycle_groups WHERE program_id = ? ORDER BY cycle_number", (program_id,)).fetchall()
    conn.close()
    return groups

def get_pcr_steps_with_groups(program_id):
    conn = get_db()
    steps = conn.execute('''
        SELECT s.*, g.cycle_number, g.repeat_count
        FROM pcr_steps s
        LEFT JOIN pcr_cycle_groups g ON s.cycle_group_id = g.id
        WHERE s.program_id = ?
        ORDER BY s.step_order
    ''', (program_id,)).fetchall()
    conn.close()
    return steps

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
