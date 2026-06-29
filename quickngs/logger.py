# quickngs/logger.py
"""
Centralized logging for QuickNGS.
Writes logs to app.log with rotation (max 1 MB, keep 3 backups).
Powered by Pourdad Panahi – Built with DeepSeek AI
"""

import logging
import os
from logging.handlers import RotatingFileHandler

def setup_logger(log_file='app.log'):
    logger = logging.getLogger('quickngs')
    logger.setLevel(logging.INFO)

    # جلوگیری از اضافه شدن چند handler تکراری
    if not logger.handlers:
        handler = RotatingFileHandler(log_file, maxBytes=1_000_000, backupCount=3)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger
