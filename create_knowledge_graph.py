from biocypher import BioCypher
from decider_genetics.adapters.all_variants_adapter import (
    AllVariantsAdapter,
    AllVariantsAdapterNodeType,
    AllVariantsAdapterEdgeType,
    AllVariantsAdapterPatientField,
    AllVariantsAdapterVariantField,
)

# Instantiate the BioCypher interface
# You can use `config/biocypher_config.yaml` to configure the framework or
# supply settings via parameters below
bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
)

# Choose node types to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_types = [
    AllVariantsAdapterNodeType.PATIENT,
    AllVariantsAdapterNodeType.VARIANT,
]

# Choose protein adapter fields to include in the knowledge graph.
# These are defined in the adapter (`adapter.py`).
node_fields = [
    # Patients
    AllVariantsAdapterPatientField.ID,
    # Variants
    AllVariantsAdapterVariantField.ID,
    AllVariantsAdapterVariantField.CHROMOSOME,
    AllVariantsAdapterVariantField.POSITION,
    AllVariantsAdapterVariantField.REF,
    AllVariantsAdapterVariantField.ALT,
    AllVariantsAdapterVariantField.GENE,
    AllVariantsAdapterVariantField.COSMIC_ID,
    AllVariantsAdapterVariantField.TRUNCAL,
]

edge_types = [
    AllVariantsAdapterEdgeType.PATIENT_VARIANT_ASSOCIATION,
]

# Create a protein adapter instance
adapter = AllVariantsAdapter(
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
)


# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Write admin import statement
data = bc.write_import_call()

# Print summary
bc.summary()
