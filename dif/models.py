from sqlalchemy import Column, Integer, VARCHAR, INTEGER
from sqlalchemy.ext.declarative import declarative_base

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


Base.metadata.create_all(bind=engine)
