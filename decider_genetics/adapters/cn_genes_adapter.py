import hashlib
import math
import random
import string
import pandas as pd
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class CnGenesAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    SAMPLE = auto()
    GENE = auto()


class CnGenesAdapterSampleField(Enum):
    """
    Define possible fields the adapter can provide for samples.
    """

    ID = "sample"


class CnGenesAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes.
    """

    ENSEMBL_ID = "ID"
    NAME = "Gene"
    CHR = "chr"
    START = "start"
    END = "end"
    STRAND = "strand"
    BAND = "band"
    TYPE = "type"  # should be prefiltered to only include 'protein_coding'


class CnGenesAdapterEdgeType(Enum):
    """
    Enum for the edge types of the adapter.
    """

    SAMPLE_GENE_ASSOCIATION = auto()


class CnGenesAdapterEdgeField(Enum):
    """
    Enum for the edge fields of the adapter.
    """

    N_PROBES_CR = "nProbesCr"
    N_PROBES_AF = "nProbesAf"
    LOG_R = "logR"
    BAF = "baf"
    N_ARAW = "nAraw"
    N_BRAW = "nBraw"
    N_MAJOR = "nMajor"
    N_MINOR = "nMinor"
    PURIFIED_LOG_R = "purifiedLogR"
    PURIFIED_BAF = "purifiedBaf"
    PURIFIED_LOH = "purifiedLoh"
    CN_STATUS = "CNstatus"
    LOH_STATUS = "LOHstatus"
    MIN_PURIFIED_LOG_R = "minPurifiedLogR"
    MAX_PURIFIED_LOG_R = "maxPurifiedLogR"
    BREAKS_IN_GENE = "breaksInGene"


class CnGenesAdapter:
    """
    Generates sample and gene nodes and edges between them.

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
        self.data = pd.read_csv("data/synthetic_cns.csv", sep="\t", header=0)

        # each sample is connected to each gene by copy number; first, filter
        # to keep only the specified node fields and edge fields
        self.data = self.data[
            [
                field.value
                for field in chain(
                    self.node_fields,
                    self.edge_fields,
                )
                if field.value in self.data.columns
            ]
        ]

        # GENES: remove all columns except the ones in CnGenesAdapterGeneField
        # and deduplicate
        self.genes = self.data[
            [
                field.value
                for field in CnGenesAdapterGeneField
                if field.value in self.data.columns
            ]
        ].drop_duplicates()

        # SAMPLES: should already be created by the all_variants adapter

        # VARIANTS: remove all columns except the ones in
        # CnGenesAdapterEdgeField, plus the sample id and gene NAME columns, and
        # deduplicate
        self.variants = self.data[
            [
                field.value
                for field in CnGenesAdapterEdgeField
                if field.value in self.data.columns
            ]
            + [
                CnGenesAdapterSampleField.ID.value,
                CnGenesAdapterGeneField.NAME.value,
            ]
        ].drop_duplicates()

        # generate an id for each variant using the md5 hash of all columns
        self.variants["VARIANT_ID"] = self.variants.apply(
            lambda row: hashlib.md5(
                "".join(
                    [str(row[column]) for column in self.variants.columns]
                ).encode("utf-8")
            ).hexdigest(),
            axis=1,
        )

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        # GENES: for each node (row), yield a 3-tuple of node id (the 'NAME'
        # column), node label (hardcode to 'gene' for now),
        # and node properties (dict of column names and values, except the
        # 'NAME')

        for _, node in self.genes.iterrows():
            node["name"] = node[CnGenesAdapterGeneField.NAME.value]
            yield (
                f"{node[CnGenesAdapterGeneField.NAME.value]}",
                "gene",
                node.drop(CnGenesAdapterGeneField.NAME.value).to_dict(),
            )

        # VARIANTS: for each node (row), yield a 3-tuple of node id (the 'VARIANT_ID'
        # column), node label (hardcode to 'copy_number_variant' for now), and
        # node properties

        for _, row in self.variants.iterrows():
            _props = row.drop(
                [
                    "VARIANT_ID",
                    CnGenesAdapterSampleField.ID.value,
                    CnGenesAdapterGeneField.NAME.value,
                ]
            ).to_dict()

            # replace 'nan' with 'NaN' in N_MAJOR and N_MINOR; otherwise, Neo4j
            # will throw an error (can't deal with 'nan')
            if math.isnan(_props[CnGenesAdapterEdgeField.N_MAJOR.value]):
                _props[CnGenesAdapterEdgeField.N_MAJOR.value] = "NaN"
            if math.isnan(_props[CnGenesAdapterEdgeField.N_MINOR.value]):
                _props[CnGenesAdapterEdgeField.N_MINOR.value] = "NaN"

            yield (
                row["VARIANT_ID"],
                "copy_number_variant",
                _props,
            )

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """

        logger.info("Generating edges.")

        # yield 5-tuple of edge id, source node id, target node id, edge label
        # (hardcode to 'copy_number_alteration' for now), and edge properties
        # (all columns except the source and target node ids)

        for _, row in self.variants.iterrows():

            # patient to variant
            yield (
                None,
                row[CnGenesAdapterSampleField.ID.value],
                row["VARIANT_ID"],
                "patient_has_copy_number_variant",
                {},
            )

            # variant to gene
            yield (
                None,
                row["VARIANT_ID"],
                f"{row[CnGenesAdapterGeneField.NAME.value]}",
                "copy_number_variant_in_gene",
                {},
            )

    def _set_types_and_fields(
        self, node_types, node_fields, edge_types, edge_fields
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in CnGenesAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    CnGenesAdapterSampleField,
                    CnGenesAdapterGeneField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in CnGenesAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain()]
