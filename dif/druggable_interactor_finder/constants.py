"""Defined string constants."""

PUBCHEM_BIOASSAY_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/accession/{}/aids/TXT"

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
