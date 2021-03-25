"""Default settings."""

import os
import logging

from pathlib import Path
from ebel_rest import connect
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, MetaData
from sqlalchemy_utils import create_database, database_exists


logger = logging.getLogger(__name__)

# Graphstore credentials
db = 'pharmacome'
db_user = 'mavo_user'
db_password = 'mavo'
db_server = 'https://graphstore.scai.fraunhofer.de'

connect(db_user, db_password, db_server, db, print_url=False)

# Similar diseases to Alzheimer's
SIMILAR_DISEASES = ["Parkinson Disease", "Amyotrophic Lateral Sclerosis", "Neurodegenerative Diseases",
                    "Huntington Disease"]

# Default Directory Paths
HOME = str(Path.home())
PROJECT_DIR = os.path.join(HOME, ".dif")
LOG_DIR = os.path.join(PROJECT_DIR, "logs")
CACHE_DIR = os.path.join(PROJECT_DIR, "cache")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

# Cached information
BIOASSAY_CACHE = os.path.join(CACHE_DIR, "bioassays.json")

# Logging Configuration
LOG_FILE_PATH = os.path.join(LOG_DIR, "dif.log")
logging.basicConfig(filename=LOG_FILE_PATH,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# SQLite
DB_PATH = os.path.join(PROJECT_DIR, f"dif.db")
CONN = f"sqlite:///{DB_PATH}"
if not database_exists(CONN):
    create_database(CONN)

engine = create_engine(CONN, convert_unicode=True)
session = sessionmaker(bind=engine)
