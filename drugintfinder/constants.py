"""Defined string constants."""

# Strings
POINTS = 'points'
PATENTS = 'patents'
IDENTIFIERS = 'identifiers'
PRODUCTS = 'products'
INTERACTORS = 'interactors'
DRUG_NAME = 'drug_name'
SYNERGY = 'synergy'
TRIALS = 'trials'
TARGET_COUNT = 'target_count'
CLINICAL_TRIALS = 'clinical_trials'
TIC = 'target_interactor_contradiction'
DAC = 'drug_action_contradiction'

# MAPPER
POSREG = 'positive_regulator'
NEGREG = 'negative_regulator'
NEUTRAL = 'neutral'

# API
PUBCHEM_BIOASSAY_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{}/aids/TXT"

###########
# Mappers #
###########

# Edge types
COMPILER = ['has_modification', 'has_product', 'reactant_in', 'acts_in', 'has_variant', 'translocates', 'includes']
CAUSAL = ['increases', 'decreases', 'directly_increases', 'directly_decreases', 'causes_no_change',
          'rate_limiting_step_of', 'regulates']
CORRELATIVE = ['association', 'no_correlation', 'positive_correlation', 'negative_correlation']
OTHER = ['has_member', 'has_members', 'has_component', 'has_components', 'equivalent_to', 'is_a', 'sub_process_of',
         'analogous_to', 'biomarker_for', 'prognostic_biomarker_for']
GENOMIC = ['transcribed_to', 'translated_to', 'orthologous']
BEL_RELATION = CAUSAL + CORRELATIVE + OTHER + GENOMIC
ALL = BEL_RELATION + COMPILER

EDGE_MAPPER = {
    'bel_relation': BEL_RELATION,
    'causal': CAUSAL,
    'correlative': CORRELATIVE,
    'other': OTHER,
    'genomic': GENOMIC,
    'compiler': COMPILER,
    'E': ALL,
}

CT_MAPPER = {
    'ongoing': [
        "Recruiting",
        "Enrolling by invitation",
        "Active, not recruiting",
        "Approved for marketing",
        "Available",
    ],
    'finished': [
        "Completed",
        "Unknown status",
        "Terminated",
        "Withdrawn",
        "Suspended",
        "No longer available",
        "Withheld",
        "Not yet recruiting",
        "Temporarily not available",
    ],
}

ACTION_MAPPER = {'activator': POSREG,
                 'aggregation inhibitor': NEGREG,
                 'agonist': POSREG,
                 'allosteric modulator': NEUTRAL,  # Can't tell if it's pos or neg
                 'antagonist': NEGREG,
                 'antibody': NEUTRAL,  # Can't tell if it's pos or neg
                 'binder': NEUTRAL,  # Can't tell if it's pos or neg
                 'binding': NEUTRAL,
                 'blocker': NEGREG,
                 'chaperone': POSREG,
                 'cofactor': POSREG,  # Cofactor = component for activity
                 'inducer': POSREG,
                 'inhibitor': NEGREG,
                 'ligand': NEUTRAL,  # Can't tell if it's agonist or anatagnoist
                 'modulator': NEUTRAL,  # Can't tell if it's pos or neg
                 'multitarget': NEGREG,  # 'multitarget' only used for Dasatinib which is an inhibitor
                 'negative modulator': NEGREG,
                 'neutralizer': NEGREG,
                 'other': NEUTRAL,
                 'other/unknown': NEUTRAL,
                 'partial agonist': POSREG,
                 'positive allosteric modulator': POSREG,
                 'potentiator': NEGREG,  # 'potentiator' only used for Pimecrolimus which is an inhibitor
                 'stabilization': POSREG,
                 'stimulator': POSREG,
                 'substrate': NEUTRAL,  # Can't tell if it's pos or neg
                 'suppressor': NEGREG,
                 'unknown': NEUTRAL,
                 'weak inhibitor': NEGREG,
                 }

###########
# Queries #
###########

IN_COUNT = "SELECT count(*) as number FROM causal WHERE in.name = '{}' AND in.pure = true AND in.@class = 'protein'"
OUT_COUNT = """SELECT count(*) as number FROM causal
WHERE out.name = '{}' AND out.pure = true AND out.@class = 'protein'"""

UNIPROT_ID = "SELECT uniprot.id FROM protein WHERE name = '{}' AND namespace = 'HGNC' and pure = true LIMIT 1"

CLINICAL_TRIALS_COUNT = "SELECT count(*) as trial_count FROM clinical_trial"
CLINICAL_TRIALS_DATA = "SELECT overall_status, trial_id, condition, mesh_conditions, drugs_in_trial FROM " \
                       "clinical_trial SKIP {} LIMIT {}"


DRUG_COUNT = "SELECT count(*) as num_drugs FROM drugbank"
DRUG_METADATA = "SELECT id as drugbank_id, name as drug_name, patents as drug_patents, " \
                "products.product as drug_products, target_symbols, clinical_trials.trial_id AS clinical_trials " \
                "FROM drugbank SKIP {} LIMIT {}"

ASSOCIATED_PATHWAYS = "SELECT count(*) FROM pathway_interaction WHERE out.name = '{}' OR in.name = '{}'"

INTERACTOR_QUERY = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
{{class:{}, as:target, WHERE:(name.toLowerCase() = '{}')}}
.inE(){{class:{} ,as:relation, where:(@class != 'causes_no_change'{})}}
.outV(){{class:bel, as:interactor}}
RETURN
pmod.type as pmod_type,
relation.@class as relation_type,
target.name as target_symbol,
target.bel as target_bel,
target.@class as target_type,
interactor.bel as interactor_bel,
interactor.name as interactor_name,
interactor.@class as interactor_type,
relation.pmid as pmid,
relation.pmc as pmc,
target.species as target_species
"""

PURE_DRUGGABLE_QUERY = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
{{class:{}, as:target, WHERE:(name.toLowerCase() = '{}')}}
.inE(){{class:causal,as:relation, where:(@class != 'causes_no_change'{})}}
.outV(){{class:bel, as:interactor}}
.inE(){{class:has_drug_target, as:drug_rel}}
.outV(){{class:drug, as:drug}}
RETURN
pmod.type as pmod_type,
relation.@class as relation_type,
relation.citation.pub_date.subString(0, 4) as rel_pub_year,
target.name as target_symbol,
target.bel as target_bel,
target.@class as target_type,
interactor.bel as interactor_bel,
interactor.name as interactor_name,
interactor.@class as interactor_type,
drug.label as drug,
drug.drugbank_id as drugbank_id,
drug.drugbank.chembl_id as chembl_id,
drug.drugbank.pubchem_cid as pubchem_id,
relation.pmid as pmid,
relation.pmc as pmc,
relation.@rid.asString() as rel_rid,
drug_rel.@rid.asString() as drug_rel_rid,
drug_rel.actions as drug_rel_actions
"""

CAPSULE_PREFIX = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
{{class:{}, as:target, WHERE:(name.toLowerCase() = '{}')}}
.inE(){{class:causal,as:relation, where:(@class != 'causes_no_change')}}
.outV(){{class:bel, as:capsule_interactor}}"""

CAPSULE_MODIFIED_PROTEIN = """
.inE('has_modified_protein', 'has_variant_protein', 'has_located_protein', 'has_fragmented_protein')
.outV(){{class:protein, as:pure_interactor, WHERE:(pure=true)}}"""

CAPSULE_HAS_PROTEIN = """.outE('has__protein')
.inV(){{class:protein, as:pure_interactor, WHERE:(pure=true)}}"""

CAPSULE_SUFFIX = """.inE(){{class:has_drug_target, as:drug_rel}}
.outV(){{class:drug, as:drug}}
RETURN
pmod.type as pmod_type,
drug.label as drug,
drug.drugbank_id as drugbank_id,
drug.drugbank.chembl_id as chembl_id,
drug.drugbank.pubchem_cid as pubchem_id,
pure_interactor.@class as interactor_type,
pure_interactor.bel as interactor_bel,
pure_interactor.name as interactor_name,
capsule_interactor.bel as capsule_interactor_bel,
capsule_interactor.@class as capsule_interactor_type,
relation.@class as relation_type,
relation.citation.pub_date.subString(0, 4) as rel_pub_year,
target.name as target_symbol,
target.bel as target_bel,
target.@class as target_type,
relation.pmid as pmid,
relation.pmc as pmc,
relation.@rid.asString() as rel_rid,
drug_rel.@rid.asString() as drug_rel_rid,
drug_rel.actions as drug_rel_actions
"""

CAPSULE_DRUGGABLE_MODIFIED = CAPSULE_PREFIX + CAPSULE_MODIFIED_PROTEIN + CAPSULE_SUFFIX
CAPSULE_DRUGGABLE_COMPLEX = CAPSULE_PREFIX + CAPSULE_HAS_PROTEIN + CAPSULE_SUFFIX
