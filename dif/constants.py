"""Defined string constants."""

# Strings
POINTS = 'points'
PATENTS = 'patents'
IDENTIFIERS = 'identifiers'
PRODUCTS = 'products'
INTERACTORS = 'interactors'
CLINICAL_TRIALS = 'clinical_trials'
TIC = 'target_interactor_contradiction'
DAC = 'drug_action_contradiction'

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
    'E': ALL
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
    ]
}

ACTION_MAPPER = {'activator': 'positive_regulator',
                 'aggregation inhibitor': 'negative_regulator',
                 'agonist': 'positive_regulator',
                 'allosteric modulator': 'neutral',  # Can't tell if it's pos or neg
                 'antagonist': 'negative_regulator',
                 'antibody': 'neutral',  # Can't tell if it's pos or neg
                 'binder': 'neutral',  # Can't tell if it's pos or neg
                 'cofactor': 'positive_regulator',  # Cofactor = component for activity
                 'inducer': 'positive_regulator',
                 'inhibitor': 'negative_regulator',
                 'ligand': 'neutral',  # Can't tell if it's agonist or anatagnoist
                 'modulator': 'neutral',  # Can't tell if it's pos or neg
                 'multitarget': 'negative_regulator',  # 'multitarget' only used for Dasatinib which is an inhibitor
                 'neutralizer': 'negative_regulator',
                 'other/unknown': 'neutral',
                 'partial agonist': 'positive_regulator',
                 'potentiator': 'negative_regulator',  # 'potentiator' only used for Pimecrolimus which is an inhibitor
                 'stabilization': 'positive_regulator',
                 'substrate': 'neutral',  # Can't tell if it's pos or neg
                 'weak inhibitor': 'negative_regulator',
                 }

###########
# Queries #
###########

IN_COUNT = "SELECT count(*) as number FROM causal WHERE in.name = '{}' AND in.pure = true AND in.@class = 'protein'"
OUT_COUNT = "SELECT count(*) as number FROM causal WHERE out.name = '{}' AND out.pure = true AND out.@class = 'protein'"

UNIPROT_ID = "SELECT uniprot.id FROM protein WHERE name = '{}' and pure = true LIMIT 1"
CLINICAL_TRIAL_FROM_DRUG = """SELECT expand(clinical_trials) FROM drugbank WHERE id = '{}'"""
ASSOCIATED_PATHWAYS = "SELECT count(*) FROM pathway_interaction WHERE out.name = '{}' OR in.name = '{}'"
PATENTS_PRODUCTS = "SELECT name, patents as drug_patents, products.product as drug_products FROM drugbank " \
                   "WHERE id in {}"

INTERACTOR_QUERY = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
{{class:{}, as:target, WHERE:(name in {})}}
.inE(){{class:{} ,as:relation, where:(@class != 'causes_no_change')}}
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
        {{class:{}, as:target, WHERE:(name in {})}}
        .inE(){{class:causal,as:relation, where:(@class != 'causes_no_change')}}
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

CAPSULE_DRUGGABLE_QUERY = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
        {{class:{}, as:target, WHERE:(name in {})}}
        .inE(){{class:causal,as:relation, where:(@class != 'causes_no_change')}}
        .outV(){{class:bel, as:capsule_interactor}}
        .bothE('has__protein', 'has_modified_protein', 'has_variant_protein', 'has_located_protein', 'has_fragmented_protein')
        .bothV(){{class:protein, as:pure_interactor, WHERE:(pure=true)}}
        .inE(){{class:has_drug_target, as:drug_rel}}
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
