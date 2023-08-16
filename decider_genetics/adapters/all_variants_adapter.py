import hashlib
import random
import string
import pandas as pd
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class AllVariantsAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PATIENT = auto()
    VARIANT = auto()


class AllVariantsAdapterPatientField(Enum):
    """
    Define possible fields the adapter can provide for patients.
    """

    ID = "patient"


class AllVariantsAdapterVariantField(Enum):
    """
    Define possible fields the adapter can provide for variants.
    """

    ID = "ID"
    CHROMOSOME = "CHROM"
    POSITION = "POS"
    REF = "REF"
    ALT = "ALT"
    FILTER = "FILTER"
    CYTOBAND = "cytoBand"
    FUNCTION = "Func.MANE"
    GENE = "Gene.MANE"
    GENE_DETAIL = "GeneDetail.MANE"
    EXONIC_FUNCTION = "ExonicFunc.MANE"
    AA_CHANGE = "AAChange.MANE"
    FUNCTION_REF = "Func.refGene"
    GENE_REF = "Gene.refGene"
    GENE_DETAIL_REF = "GeneDetail.refGene"
    EXONIC_FUNCTION_REF = "ExonicFunc.refGene"
    AA_CHANGE_REF = "AAChange.refGene"
    GENOMIC_SUPER_DUPS = "genomicSuperDups"
    DBSC_SNV_ADA_SCORE = "dbscSNV_ADA_SCORE"
    DBSC_SNV_RF_SCORE = "dbscSNV_RF_SCORE"
    COSMIC_ID = "COSMIC_ID"
    COSMIC_OCCURENCE = "COSMIC_OCCURRENCE"
    COSMIC_TOTAL_OCCURENCE = "COSMIC_TOTAL_OCC"
    COSMIC_CONF_SOMATIC = "COSMIC_CONF_SOMA"
    CLNSIG = "CLNSIG"
    CLNSIGCONF = "CLNSIGCONF"
    CLNDN = "CLNDN"
    CLNREVSTAT = "CLNREVSTAT"
    CLNALLELEID = "CLNALLELEID"
    CLNDISDB = "CLNDISDB"
    INTERPRO_DOMAIN = "Interpro_domain"
    REGULOME_DB = "regulomeDB"
    CADD_RAW = "CADD_raw"
    CADD_PHRED = "CADD_phred"
    THOUSAND_GENOMES_ALL = "1000G_ALL"
    THOUSAND_GENOMES_EUR = "1000G_EUR"
    GNOMAD_GENOME_ALL = "gnomAD_genome_ALL"
    GNOMAD_GENOME_NFE = "gnomAD_genome_NFE"
    GNOMAD_GENOME_FIN = "gnomAD_genome_FIN"
    GNOMAD_GENOME_MAX = "gnomAD_genome_max"
    GNOMAD_EXOME_NC_ALL = "gnomAD_exome_nc_ALL"
    GNOMAD_EXOME_NC_NFE = "gnomAD_exome_nc_NFE"
    GNOMAD_EXOME_NC_NFE_SWE = "gnomAD_exome_nc_NFE_SWE"
    GNOMAD_EXOME_NC_NC_FIN = "gnomAD_exome_nc_FIN"
    GNOMAD_EXOME_NC_MAX = "gnomAD_exome_nc_max"
    TRUNCAL = "Truncal"
    READ_COUNTS = "readCounts"
    SAMPLES = "samples"


class AllVariantsAdapterEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """

    PATIENT_VARIANT_ASSOCIATION = auto()


class AllVariantsAdapter:
    """
    Generates patient and variant nodes and edges between them.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self._set_types_and_fields(
            node_types, node_fields, edge_types, edge_fields
        )
        self._load_data()

    def _load_data(self):
        """
        Read CSV and create dataframes for variants and patients according to
        the selected fields.
        """
        logger.info("Loading data.")

        # read from csv 'data/all_variants.sampled.csv'
        self.data = pd.read_csv(
            "data/all_variants.sampled.csv", sep="\t", header=0
        )

        # one row is one variant node, first select the columns given by the
        # node_fields parameter
        self.variants = self.data[
            [
                field.value
                for field in self.node_fields
                if field.value in self.data.columns
            ]
        ]

        # PATIENTS: select the PATIENT.ID column and drop duplicates
        if AllVariantsAdapterNodeType.PATIENT in self.node_types:
            self.patients = self.data[[AllVariantsAdapterPatientField.ID.value]].drop_duplicates()

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        for _, patient in self.patients.iterrows():
            yield (
                patient[AllVariantsAdapterPatientField.ID.value],
                "patient",
                {},
            )

        # VARIANTS: for each node (row), yield a 3-tuple of node id (the 'ID'
        # column), node label (hardcode to 'variant' for now), and node
        # properties (dict of column names and values, except the 'ID')

        for _, node in self.data.iterrows():
            # if ID is '.', generate md5 hash from other columns
            if node["ID"] == ".":
                node["ID"] = hashlib.md5(
                    "".join(
                        [
                            str(node[column])
                            for column in self.data.columns
                            if column != "ID"
                        ]
                    ).encode("utf-8")
                ).hexdigest()

            yield (
                node["ID"],
                "variant",
                node.drop("ID").to_dict(),
            )

    def get_edges(self, probability: float = 0.3):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two nodes.
        """

        logger.info("Generating edges.")

        for node in self.variants:
            if random.random() < probability:
                other_node = random.choice(self.variants)

                # generate random relationship id by choosing upper or lower letters and integers, length 10, and joining them
                relationship_id = "".join(
                    random.choice(string.ascii_letters + string.digits)
                    for _ in range(10)
                )

                # determine type of edge from other_node type
                if (
                    isinstance(other_node, Protein)
                    and AllVariantsAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION
                    in self.edge_types
                ):
                    edge_type = (
                        AllVariantsAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value
                    )
                elif (
                    isinstance(other_node, Disease)
                    and AllVariantsAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION
                    in self.edge_types
                ):
                    edge_type = (
                        AllVariantsAdapterEdgeType.PROTEIN_DISEASE_ASSOCIATION.value
                    )
                else:
                    continue

                yield (
                    relationship_id,
                    node.get_id(),
                    other_node.get_id(),
                    edge_type,
                    {"example_proptery": "example_value"},
                )

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

    def _set_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in AllVariantsAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    AllVariantsAdapterPatientField,
                    AllVariantsAdapterVariantField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in AllVariantsAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain()]


class Node:
    """
    Base class for nodes.
    """

    def __init__(self):
        self.id = None
        self.label = None
        self.properties = {}

    def get_id(self):
        """
        Returns the node id.
        """
        return self.id

    def get_label(self):
        """
        Returns the node label.
        """
        return self.label

    def get_properties(self):
        """
        Returns the node properties.
        """
        return self.properties


class Protein(Node):
    """
    Generates instances of proteins.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "uniprot_protein"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random UniProt-style id.
        """
        lets = [random.choice(string.ascii_uppercase) for _ in range(3)]
        nums = [random.choice(string.digits) for _ in range(3)]

        # join alternating between lets and nums
        return "".join([x for y in zip(lets, nums) for x in y])

    def _generate_properties(self):
        properties = {}

        ## random amino acid sequence
        if (
            self.fields is not None
            and AllVariantsAdapterPatientField.SEQUENCE in self.fields
        ):

            # random int between 50 and 250
            l = random.randint(50, 250)

            properties["sequence"] = "".join(
                [random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(l)],
            )

        ## random description
        if (
            self.fields is not None
            and AllVariantsAdapterPatientField.DESCRIPTION in self.fields
        ):
            properties["description"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        ## taxon
        if (
            self.fields is not None
            and AllVariantsAdapterPatientField.TAXON in self.fields
        ):
            properties["taxon"] = "9606"

        return properties


class Disease(Node):
    """
    Generates instances of diseases.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "do_disease"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random disease id.
        """
        nums = [random.choice(string.digits) for _ in range(8)]

        return f"DOID:{''.join(nums)}"

    def _generate_properties(self):
        properties = {}

        ## random name
        if (
            self.fields is not None
            and AllVariantsAdapterVariantField.NAME in self.fields
        ):
            properties["name"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        ## random description
        if (
            self.fields is not None
            and AllVariantsAdapterVariantField.DESCRIPTION in self.fields
        ):
            properties["description"] = " ".join(
                [random.choice(string.ascii_lowercase) for _ in range(10)],
            )

        return properties
