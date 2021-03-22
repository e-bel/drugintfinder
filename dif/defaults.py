"""Default settings."""

import os
import logging

from pathlib import Path
from ebel_rest import connect

logger = logging.getLogger(__name__)

# Database credentials
db = 'pharmacome'
db_user = 'mavo_user'
db_password = 'mavo'
db_server = 'https://graphstore.scai.fraunhofer.de'

connect(db_user, db_password, db_server, db, print_url=False)

# Default Directory Paths
home_dir = str(Path.home())
PROJECT_DIR = os.path.join(home_dir, ".dif")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
DATA_DIR = os.path.join(PROJECT_DIR, "data")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "dif.log")
logging.basicConfig(filename=LOG_FILE_PATH,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
