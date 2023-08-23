from biocypher import BioCypher
from decider_genetics.adapters.all_variants_adapter import (
    AllVariantsAdapter,
    AllVariantsAdapterNodeType,
    AllVariantsAdapterEdgeType,
    AllVariantsAdapterPatientField,
    AllVariantsAdapterSampleField,
    AllVariantsAdapterVariantField,
)
from decider_genetics.adapters.cn_genes_adapter import (
    CnGenesAdapter,
    CnGenesAdapterNodeType,
    CnGenesAdapterEdgeType,
    CnGenesAdapterSampleField,
    CnGenesAdapterGeneField,
    CnGenesAdapterEdgeField,
)

bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
)

# VARIANTS from all_variants.csv
variant_node_types = [
    AllVariantsAdapterNodeType.PATIENT,
    AllVariantsAdapterNodeType.SAMPLE,
    AllVariantsAdapterNodeType.VARIANT,
]

variant_node_fields = [
    # Patients
    AllVariantsAdapterPatientField.ID,
    AllVariantsAdapterPatientField.SAMPLES,
    # Samples
    AllVariantsAdapterSampleField.ID,
    AllVariantsAdapterSampleField.READ_COUNTS,
    # Variants
    AllVariantsAdapterVariantField.ID,
    AllVariantsAdapterVariantField.CHROMOSOME,
    AllVariantsAdapterVariantField.POSITION,
    AllVariantsAdapterVariantField.REF,
    AllVariantsAdapterVariantField.ALT,
    AllVariantsAdapterVariantField.GENE,
    AllVariantsAdapterVariantField.CADD_PHRED,
    AllVariantsAdapterVariantField.FUNCTION,
    AllVariantsAdapterVariantField.EXONIC_FUNCTION,
    AllVariantsAdapterVariantField.AA_CHANGE,
    AllVariantsAdapterVariantField.COSMIC_TOTAL_OCCURENCE,
    AllVariantsAdapterVariantField.CLNSIG,
    AllVariantsAdapterVariantField.CLNREVSTAT,
    AllVariantsAdapterVariantField.GNOMAD_GENOME_MAX,
]

variant_edge_types = [
    AllVariantsAdapterEdgeType.PATIENT_SAMPLE_ASSOCIATION,
    AllVariantsAdapterEdgeType.SAMPLE_VARIANT_ASSOCIATION,
    AllVariantsAdapterEdgeType.VARIANT_GENE_ASSOCIATION,
]

variant_adapter = AllVariantsAdapter(
    node_types=variant_node_types,
    node_fields=variant_node_fields,
    edge_types=variant_edge_types,
)

# COPY NUMBERS from CnCombinedGenes.csv
cn_node_types = [
    CnGenesAdapterNodeType.SAMPLE,
    CnGenesAdapterNodeType.GENE,
]

cn_node_fields = [
    # Samples
    CnGenesAdapterSampleField.ID,
    # Genes
    CnGenesAdapterGeneField.ENSEMBL_ID,
    CnGenesAdapterGeneField.NAME,
    CnGenesAdapterGeneField.CHR,
    CnGenesAdapterGeneField.START,
    CnGenesAdapterGeneField.END,
    CnGenesAdapterGeneField.STRAND,
    CnGenesAdapterGeneField.BAND,
    CnGenesAdapterGeneField.TYPE,
]

cn_edge_types = [
    CnGenesAdapterEdgeType.SAMPLE_GENE_ASSOCIATION,
]

cn_edge_fields = [
    CnGenesAdapterEdgeField.BREAKS_IN_GENE,
    CnGenesAdapterEdgeField.N_MAJOR,
    CnGenesAdapterEdgeField.N_MINOR,
    CnGenesAdapterEdgeField.PURIFIED_LOG_R,
    CnGenesAdapterEdgeField.MIN_PURIFIED_LOG_R,
    CnGenesAdapterEdgeField.MAX_PURIFIED_LOG_R,
    CnGenesAdapterEdgeField.PURIFIED_BAF,
    CnGenesAdapterEdgeField.PURIFIED_LOH,
]

cn_adapter = CnGenesAdapter(
    node_types=cn_node_types,
    node_fields=cn_node_fields,
    edge_types=cn_edge_types,
    edge_fields=cn_edge_fields,
)

# Create a knowledge graph from the adapters
bc.write_nodes(variant_adapter.get_nodes())
bc.write_nodes(cn_adapter.get_nodes())

bc.write_edges(variant_adapter.get_edges())
bc.write_edges(cn_adapter.get_edges())

# Write admin import statement
data = bc.write_import_call()

# Print summary
bc.summary()
