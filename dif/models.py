"""Table definitions for SQLite DB."""

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, VARCHAR, INTEGER, BOOLEAN

from dif.defaults import engine

Base = declarative_base()


class MetaClass:

    def dump(self):
        return dict([(k, v) for k, v in vars(self).items() if not k.startswith('_')])


class Druggable(MetaClass, Base):
    __tablename__ = 'druggable'
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
    __tablename__ = 'patents'
    id = Column(Integer, primary_key=True)

    drug_name = Column(VARCHAR(255))
    has_patent = Column(BOOLEAN)
    expired = Column(BOOLEAN)
    patent_numbers = Column(VARCHAR(255))


class Products(MetaClass, Base):
    __tablename__ = 'products'
    id = Column(Integer, primary_key=True)

    drug_name = Column(VARCHAR(255))
    has_generic = Column(BOOLEAN)
    has_approved_generic = Column(BOOLEAN)
    generic_products = Column(VARCHAR(255))


class Trials(MetaClass, Base):
    __tablename__ = 'trials'
    id = Column(Integer, primary_key=True)

    drug_name = Column(VARCHAR(255), index=True)
    drugbank_id = Column(VARCHAR(255), index=True)
    trial_id = Column(VARCHAR(255), index=True)
    trial_status = Column(VARCHAR(255))
    conditions = Column(VARCHAR(255))
    drugs_in_trial = Column(VARCHAR(255))


class Edges(MetaClass, Base):
    __tablename__ = 'edges'
    id = Column(Integer, primary_key=True)

    symbol = Column(VARCHAR(255), index=True)
    out_count = Column(INTEGER)
    in_count = Column(INTEGER)
    both_count = Column(INTEGER)


Base.metadata.create_all(engine)
