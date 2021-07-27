"""Table definitions for SQLite DB."""

from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, VARCHAR, INTEGER, BOOLEAN, ForeignKey

from drugintfinder.defaults import engine

Base = declarative_base()


class MetaClass:
    """Meta methods for handling DB tables."""

    def dump(self):
        """Delete all tables in DB."""
        return dict([(k, v) for k, v in vars(self).items() if not k.startswith('_')])


class Druggable(MetaClass, Base):
    """Druggable interactor table definitions."""

    __tablename__ = 'druggable_interactor_relations'
    id = Column(Integer, primary_key=True)

    drug = Column(VARCHAR(1000))
    drugbank_id = Column(VARCHAR(255))
    chembl_id = Column(VARCHAR(255))
    pubchem_id = Column(VARCHAR(255))
    interactor_type = Column(VARCHAR(100))
    interactor_bel = Column(VARCHAR(1000))
    interactor_name = Column(VARCHAR(255))
    capsule_interactor_bel = Column(VARCHAR(1000))
    capsule_interactor_type = Column(VARCHAR(255))
    rel_pub_year = Column(INTEGER)
    target_bel = Column(VARCHAR(1000))
    target_symbol = Column(VARCHAR(255), index=True)
    target_type = Column(VARCHAR(255), index=True)
    relation_type = Column(VARCHAR(255))
    pmod_type = Column(VARCHAR(255))
    pmid = Column(INTEGER)
    pmc = Column(VARCHAR(255))
    rel_rid = Column(VARCHAR(255))
    drug_rel_rid = Column(VARCHAR(255))
    drug_rel_actions = Column(VARCHAR(255))


class General(MetaClass, Base):
    """General table definitions."""

    __tablename__ = 'general'
    id = Column(Integer, primary_key=True)

    pmod_type = Column(VARCHAR(255))
    target_bel = Column(VARCHAR(1000))
    target_symbol = Column(VARCHAR(255), index=True)
    target_type = Column(VARCHAR(255), index=True)
    relation_type = Column(VARCHAR(255))
    interactor_bel = Column(VARCHAR(1000))
    interactor_name = Column(VARCHAR(255))
    interactor_type = Column(VARCHAR(255))
    pmid = Column(INTEGER)
    pmc = Column(VARCHAR(255))
    target_species = Column(INTEGER)


class Patents(MetaClass, Base):
    """Patents table definitions."""

    __tablename__ = 'patents'
    id = Column(Integer, primary_key=True)

    expired = Column(BOOLEAN)
    patent_number = Column(Integer)
    drug_id = Column(Integer, ForeignKey('drugs.id'), index=True)


class Products(MetaClass, Base):
    """Products table definitions."""

    __tablename__ = 'products'
    id = Column(Integer, primary_key=True)

    product_name = Column(VARCHAR(255), index=True)
    has_generic = Column(BOOLEAN)
    is_approved = Column(BOOLEAN)
    has_approved_generic = Column(BOOLEAN)
    drug_id = Column(Integer, ForeignKey('drugs.id'), index=True)


class Trials(MetaClass, Base):
    """Clinical trials table definitions."""

    __tablename__ = 'trials'
    id = Column(Integer, primary_key=True, index=True)

    trial_id = Column(VARCHAR(255), index=True, unique=True)
    trial_status = Column(VARCHAR(255))
    conditions = Column(VARCHAR(255))
    drugs_in_trial = Column(VARCHAR(255))


class TargetStats(MetaClass, Base):
    """Target statistic table definitions."""

    __tablename__ = 'target_stats'
    id = Column(Integer, primary_key=True)

    symbol = Column(VARCHAR(255), index=True)
    out_count = Column(INTEGER)
    in_count = Column(INTEGER)
    both_count = Column(INTEGER)


class Drugs(MetaClass, Base):
    """Drugs table definitions."""

    __tablename__ = 'drugs'
    id = Column(Integer, primary_key=True, index=True)

    drug_name = Column(VARCHAR(255), index=True)
    drugbank_id = Column(VARCHAR(255), index=True, unique=True)
    num_targets = Column(INTEGER)
    targets = Column(VARCHAR(255))
    patents = relationship("Patents")
    products = relationship("Products")
    clinical_trials = Column(VARCHAR(255))


class BioAssays(MetaClass, Base):
    """BioAssay table definitions."""

    __tablename__ = 'bioassays'
    id = Column(Integer, primary_key=True)

    target = Column(VARCHAR(255), index=True, unique=True)
    num_assays = Column(INTEGER)


Base.metadata.create_all(engine)
