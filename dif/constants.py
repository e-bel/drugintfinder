"""Defined string constants."""

PUBCHEM_BIOASSAY_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{}/aids/TXT"


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
                 'potentiator': 'negative_regulator',
                 # 'potentiator' only used for Pimecrolimus which is an inhibitor
                 'stabilization': 'positive_regulator',
                 'substrate': 'neutral',  # Can't tell if it's pos or neg
                 'weak inhibitor': 'negative_regulator',
                 }

###########
# Queries #
###########

UNIPROT_ID = "SELECT uniprot.id FROM protein WHERE name = '{}'"
CLINICAL_TRIAL_FROM_DRUG = """SELECT expand(clinical_trials) FROM drugbank WHERE id = '{}'"""
ASSOCIATED_PATHWAYS = "SELECT count(*) FROM pathway_interaction WHERE out.name = '{}' OR in.name = '{}'"

INTERACTOR_QUERY = """MATCH {{class:pmod, as:pmod{}}}<-has__pmod-
{{class:{}, as:target, WHERE:(name in {})}}
.inE(){{class:{} ,as:relation, where:(@class != 'causes_no_change')}}
.outV(){{class:bel, as:interactor}}
RETURN
pmod.type as pmod_type,
relation.@class as relation_type,
target.bel as target_bel,
interactor.bel as interactor_bel,
interactor.name as interactor_name,
interactor.@class as interactor_type,
interactor.involved_genes as interactor_involved_genes,
interactor.involved_other as interactor_involved_other,
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
        relation.@class as relation_type,
        target.bel as target_bel,
        interactor.bel as interactor_bel,
        interactor.name as interactor_name,
        interactor.@class as interactor_type,
        drug.label as drug,
        relation.pmid as pmid,
        relation.pmc as pmc,
        relation.@rid.asString() as rel_rid,
        drug_rel.@rid.asString() as drug_rel_rid
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
        drug.label as drug,
        pure_interactor.@class as interactor_type,
        pure_interactor.bel as interactor_bel,
        pure_interactor.name as interactor_name,
        capsule_interactor.bel as capsule_interactor_bel,
        capsule_interactor.@class as capsule_interactor_type,
        relation.@class as relation_type,
        target.bel as target_bel,
        relation.pmid as pmid,
        relation.pmc as pmc,
        relation.@rid.asString() as rel_rid,
        drug_rel.@rid.asString() as drug_rel_rid
        """
