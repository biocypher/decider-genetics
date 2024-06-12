import hashlib
import pandas as pd
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class PandasAdapter:
    """
    Transforms custom (demo) data from oncodashkb.
    """

    def __init__(
        self,
    ):
        self._load_data()

    def _load_data(self):
        logger.info("Loading data.")

        self.nodes = pd.DataFrame()
        self.edges = pd.DataFrame()

        # read from csv
        bio_process = pd.read_csv(
            "data/oncodash files/BiologicalProcess-part000.csv",
            sep=";",
            names=[
                "id",
                "name",
                "preferred_id",
                "label",
            ],
        )

        bio_process.loc[:, "id"] = bio_process["id"].str.replace(
            ":biological_process", ""
        )
        bio_process["name"] = bio_process.loc[:, "id"]
        bio_process.loc[:, "label"] = "biological_process"
        bio_process = bio_process.drop(["preferred_id"], axis=1)
        bio_process = bio_process[bio_process["id"] != "None"]

        self.nodes = pd.concat([self.nodes, bio_process])

        # read from csv
        gene_to_process = pd.read_csv(
            "data/oncodash files/GeneToBiologicalProcess-part000.csv",
            sep=";",
            names=[
                "Gene",
                "BiologicalProcess",
                "Label",
            ],
        )
        gene_to_process.loc[:, "Gene"] = gene_to_process["Gene"].str.replace(
            ":gene_hugo", ""
        )
        gene_to_process.loc[:, "BiologicalProcess"] = gene_to_process[
            "BiologicalProcess"
        ].str.replace(":biological_process", "")
        gene_to_process.loc[:, "Label"] = "gene_to_process"
        gene_to_process = gene_to_process[
            gene_to_process["BiologicalProcess"] != "None"
        ]

        self.edges = pd.concat([self.edges, gene_to_process])

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """

        logger.info("Generating nodes.")

        for _, row in self.nodes.iterrows():
            properties = row.drop(["id", "label"]).to_dict()
            yield (
                row["id"],
                row["label"],
                properties,
            )

    def get_edges(self):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.
        """

        logger.info("Generating edges.")

        for _, row in self.edges.iterrows():
            p_id = row["Gene"]
            s_id = row["BiologicalProcess"]
            _id = hashlib.md5((p_id + s_id).encode("utf-8")).hexdigest()
            yield (
                _id,
                p_id,
                s_id,
                "gene_to_process",
                {},
            )
