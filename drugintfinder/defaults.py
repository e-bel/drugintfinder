"""Default settings."""

import os
import logging
from pathlib import Path

from ebel_rest import connect
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy_utils import create_database, database_exists


logger = logging.getLogger(__name__)

# Graphstore credentials
db = 'pharmacome'
db_user = 'guest'
db_password = 'guest'
db_server = 'https://graphstore.scai.fraunhofer.de'

connect(db_user, db_password, db_server, db, print_url=False)

# Similar diseases to Alzheimer's
SIMILAR_DISEASES = ["Parkinson Disease", "Amyotrophic Lateral Sclerosis", "Neurodegenerative Diseases",
                    "Huntington Disease"]

# Default Directory Paths
HOME = str(Path.home())
PROJECT_DIR = os.path.join(HOME, ".drugintfinder")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
CACHE_DIR = os.path.join(PROJECT_DIR, "cache")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "drugintfinder.log")
logging.basicConfig(filename=LOG_FILE_PATH,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# SQLite
DEFAULT_DB_PATH = os.path.join(PROJECT_DIR, "drugintfinder.db")
CONN_ENV = os.getenv('DIF_DB_PATH')
CONN_PATH = CONN_ENV if CONN_ENV else DEFAULT_DB_PATH
CONN = f"sqlite:///{CONN_PATH}"
if not database_exists(CONN):
    create_database(CONN)

engine = create_engine(CONN, convert_unicode=True)
session = sessionmaker(bind=engine)
