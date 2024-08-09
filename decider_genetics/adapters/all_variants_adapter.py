import hashlib
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
    SAMPLE = auto()
    VARIANT = auto()


class AllVariantsAdapterPatientField(Enum):
    """
    Define possible fields the adapter can provide for patients.
    """

    ID = "patient"
    SAMPLES = "samples"


class AllVariantsAdapterSampleField(Enum):
    """
    Define possible fields the adapter can provide for samples.
    """

    ID = "sample"
    READ_COUNTS = "readCounts"


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


class AllVariantsAdapterEdgeType(Enum):
    """
    Enum for the edge types of the adapter.
    """

    PATIENT_SAMPLE_ASSOCIATION = auto()
    SAMPLE_VARIANT_ASSOCIATION = auto()
    VARIANT_GENE_ASSOCIATION = auto()


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

        # read from csv
        data = pd.read_csv("data/synthetic_variants.csv", sep="\t", header=0)

        # one row is one variant node, first select the columns given by the
        # node_fields parameter
        self.variants = data[
            [
                field.value
                for field in self.node_fields
                if field.value in data.columns
            ]
        ]

        # break up the 'samples' column into one row per sample, rename the
        # column to 'sample'; at the same time, the readCounts and Gene.MANE
        # columns need to be split into one row per sample as well
        self.variants = self.variants.assign(
            samples=self.variants["samples"].str.split(";")
        ).explode("samples")
        self.variants = self.variants.rename(
            columns={"samples": "sample"}
        ).reset_index(drop=True)
        self.variants = self.variants.assign(
            readCounts=self.variants["readCounts"].str.split(";")
        ).explode("readCounts")
        self.variants = self.variants.assign(
            Gene=self.variants["Gene.MANE"].str.split(";")
        ).explode("Gene")

        # remove duplicate rows
        self.variants = self.variants.drop_duplicates()

        # define dropped columns that should not contribute to hash
        self._drop_columns = [
            AllVariantsAdapterPatientField.ID.value,
            AllVariantsAdapterSampleField.ID.value,
            AllVariantsAdapterSampleField.READ_COUNTS.value,
            "Gene",
        ]

        # if ID is '.', generate md5 hash from other columns
        self.variants["ID"] = self.variants.apply(
            lambda row: (
                hashlib.md5(
                    "".join(
                        [
                            str(row[column])
                            for column in self.variants.columns
                            if column not in self._drop_columns
                        ]
                    ).encode("utf-8")
                ).hexdigest()
                if row["ID"] == "."
                else row["ID"]
            ),
            axis=1,
        )

        # PATIENTS and SAMPLES: select the PATIENT.ID and SAMPLE.ID column and
        # drop duplicates
        if AllVariantsAdapterNodeType.PATIENT in self.node_types:
            self.patients = self.variants[
                [
                    AllVariantsAdapterPatientField.ID.value,
                    AllVariantsAdapterSampleField.ID.value,
                ]
            ].drop_duplicates()

        # GENES: should already be created by the copy number adapter

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        # # get unique patients as list
        # patients = (
        #     self.patients[[AllVariantsAdapterPatientField.ID.value]]
        #     .drop_duplicates()[AllVariantsAdapterPatientField.ID.value]
        #     .tolist()
        # )

        # # get unique samples
        # samples = (
        #     self.variants[[AllVariantsAdapterSampleField.ID.value]]
        #     .drop_duplicates()[AllVariantsAdapterSampleField.ID.value]
        #     .tolist()
        # )

        # # yield tuples for patients
        # for patient in patients:
        #     yield (
        #         patient,
        #         "patient",
        #         {"name": patient},
        #     )

        # # yield tuples for samples
        # for sample in samples:
        #     yield (
        #         sample,
        #         "sample",
        #         {},
        #     )

        # VARIANTS: for each node (row), yield a 3-tuple of node id (the 'ID'
        # column), node label (hardcode to 'variant' for now), and node
        # properties (dict of column names and values, except the 'ID')

        # first remove columns in _drop_columns and drop duplicates
        unique_variants = self.variants.drop(
            self._drop_columns, axis=1
        ).drop_duplicates()

        for _, node in unique_variants.iterrows():
            yield (
                node["ID"],
                "variant",
                node.drop("ID").to_dict(),
            )

    def get_edges(self, variant_via_sample: bool = False):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """

        logger.info("Generating edges.")

        if variant_via_sample:
            # PATIENT - SAMPLE
            # yield 5-tuple of edge id (hash of patient and sample ids), source
            # node id, target node id, edge label (hardcode to 'patient_has_sample'
            # for now), and edge properties (empty dict for now)
            for _, row in self.patients.iterrows():
                p_id = row[AllVariantsAdapterPatientField.ID.value]
                s_id = row[AllVariantsAdapterSampleField.ID.value]
                _id = hashlib.md5((p_id + s_id).encode("utf-8")).hexdigest()
                yield (
                    _id,
                    p_id,
                    s_id,
                    "patient_has_sample",
                    {},
                )

            # SAMPLE - VARIANT
            # yield 5-tuple of edge id (hash of sample and variant ids), source node
            # id, target node id, edge label (hardcode to 'sample_has_variant' for
            # now), and edge properties (empty dict for now)
            for _, row in self.variants.iterrows():
                s_id = row[AllVariantsAdapterSampleField.ID.value]
                v_id = row[AllVariantsAdapterVariantField.ID.value]
                _id = hashlib.md5((s_id + v_id).encode("utf-8")).hexdigest()
                yield (
                    _id,
                    s_id,
                    v_id,
                    "sample_has_variant",
                    {},
                )

        else:
            # PATIENT - VARIANT
            # yield 5-tuple of edge id (hash of patient and variant ids), source node
            # id, target node id, edge label (hardcode to 'patient_has_variant' for
            # now), and edge properties (empty dict for now)
            for _, row in self.variants.iterrows():
                p_id = row[AllVariantsAdapterPatientField.ID.value]
                v_id = row[AllVariantsAdapterVariantField.ID.value]
                _id = hashlib.md5((p_id + v_id).encode("utf-8")).hexdigest()
                yield (
                    _id,
                    p_id,
                    v_id,
                    "patient_has_variant",
                    {},
                )

        # VARIANT - GENE
        # yield 5-tuple of edge id, source node id, target node id, edge label
        # (hardcode to 'variant_in_gene' for now), and edge properties (empty
        # dict for now)

        # first remove all columns except the variant id and gene and
        # deduplicate
        unique_variant_gene = self.variants[
            [
                AllVariantsAdapterVariantField.ID.value,
                "Gene",  # expanded above
            ]
        ].drop_duplicates()

        # remove all 'NONE' genes
        unique_variant_gene = unique_variant_gene[
            unique_variant_gene["Gene"] != "NONE"
        ]

        for _, row in unique_variant_gene.iterrows():
            v_id = row[AllVariantsAdapterVariantField.ID.value]
            g_id = f"{row['Gene']}"
            _id = hashlib.md5((v_id + g_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                v_id,
                g_id,
                "variant_in_gene",
                {},
            )

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
