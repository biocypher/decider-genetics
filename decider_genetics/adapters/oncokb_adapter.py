import hashlib
import pandas as pd
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class OncoKBAdapter:
    """
    Load OncoKB druggability data.
    """

    def __init__(self) -> None:
        self._load_data()

    def _load_data(self) -> None:
        logger.info("Loading data.")

        # read from csv
        raw_df = pd.read_csv(
            "data/oncokb_biomarker_drug_associations.tsv",
            sep="\t",
            header=0,
        )

        # explode the "Drugs (for therapeutic implications only)" column
        raw_df = raw_df.assign(
            **{
                "Drugs (for therapeutic implications only)": raw_df[
                    "Drugs (for therapeutic implications only)"
                ].str.split(", ")
            }
        ).explode("Drugs (for therapeutic implications only)")

        # explode the "Alterations" column
        raw_df = raw_df.assign(
            **{"Alterations": raw_df["Alterations"].str.split(", ")}
        ).explode("Alterations")

        # explode the "Cancer Types" column
        raw_df = raw_df.assign(
            **{"Cancer Types": raw_df["Cancer Types"].str.split(", ")}
        ).explode("Cancer Types")

        # remove drugs that are not strings
        raw_df = raw_df[
            raw_df["Drugs (for therapeutic implications only)"].apply(
                lambda x: isinstance(x, str)
            )
        ]

        self._data = raw_df

        logger.info("Data loaded.")

    def get_nodes(self):
        # drugs
        drugs = (
            self._data["Drugs (for therapeutic implications only)"]
            .drop_duplicates()
            .tolist()
        )
        for drug in drugs:
            yield (
                drug,
                "drug",
                {"id": drug, "name": drug},
            )

    def get_edges(self):
        # gene druggability
        for _, row in self._data.iterrows():
            row_id = hashlib.sha256(str(row).encode("utf-8")).hexdigest()
            yield (
                row_id,
                row["Gene"],
                row["Drugs (for therapeutic implications only)"],
                "potentially_druggable",
                {
                    "level": row["Level"],
                    "alteration": row["Alterations"],
                    "cancer_type": row["Cancer Types"],
                },
            )
