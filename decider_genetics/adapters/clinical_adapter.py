import hashlib
import random
import pandas as pd
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

# columns:
# Attention
# Age at Diagnosis
# BMI at Dg
# Histology
# Histology re_evaluated in DECICER
# Histology re_evaluated notes
# Stage_FIGO2014
# Treatment strategy
# Residual tumor PDS
# Residual tumor IDS
# Oper2_IDS cancelled
# Oper2_laparoscopic evaluation
# Oper2_explorative laparotomy
# Oper2_debulking surgery
# Oper1_Omental disease largest nodule_NEW
# Oper2_Omental disease largest nodule_NEW
# Primary chemotherapy cycles
# NACT cycles
# Post IDS chemotherapy cycles
# Maintenance therary after 1st line
# PARPi treatment
# Patient card::Participation in clinical trials
# Patient card::DrugTrial_name
# Patient card::Drugtrial unblinded
# CRS Omental
# CRS Omental range
# CRS3 residual tumor
# RECIST 1.1 Response to NACT
# Response to NACT radiologic evaluation made
# Primary therapy outcome
# Current phase of treatment
# Progression Yes_No_ND
# Survival
# Cause of death
# OS_KaplanM_allHGSC
# PFS_KaplanM_allHGSC
# PFI_KaplanM_allHGSC
# Time from End of 1st line maintenance to 1st prog_Days TFI
# Time from 1st prog to Death_Days Post progression survival
# BRCA mutation any
# BRCA mutation status in Clinical test
# BRCA any in Decider
# HR signature SBS3 pretreatment WGS
# HR signature SBS3 per patient
# HRD Clinical test result
# Chronic illnesses at Dg
# Chronic illnesses type
# Previous cancer yes no
# Previous cancer dg
# Previous cancer_year
# Excluded from PFI_calculations
# Excluded from PFI calculations_ reason


class ClinicalAdapter:
    """
    Load clinical patient data.
    """

    def __init__(self) -> None:
        self._load_data()

    def _load_data(self) -> None:
        logger.info("Loading data.")

        # read from csv
        self.nodes = pd.read_csv(
            "data/synthetic_clinical.csv",
            sep=";",
            header=0,
        )

    def get_nodes(self):
        """
        Iterate through the rows and create a node for each patient, yielding a
        tuple of name, "patient", and dictionary of properties (other columns).
        """

        for _, row in self.nodes.iterrows():
            # create a unique identifier for the patient
            patient_id = row["Patient"]
            _props = row.drop("Patient").to_dict()
            # lowercase the keys, replace space with underscore
            _props = {k.lower().replace(" ", "_"): v for k, v in _props.items()}
            _props["name"] = row["Patient"]
            _props["bmi"] = float(_props["bmi"].replace(",", "."))
            _props["chemotherapy_cycles"] = int(_props["chemotherapy_cycles"])

            # convert parpi, brca_mutation, and hr_deficient to boolean
            _props["parpi"] = str(_props["parpi"]).lower() == "yes"
            _props["brca_mutation"] = (
                str(_props["brca_mutation"]).lower() == "yes"
            )
            _props["hr_deficient"] = (
                str(_props["hr_deficient"]).lower() == "hrd positive"
            )
            # convert bool to lowercase string
            _props["parpi"] = str(_props["parpi"]).lower()
            _props["brca_mutation"] = str(_props["brca_mutation"]).lower()
            _props["hr_deficient"] = str(_props["hr_deficient"]).lower()

            # add fake severe_adverse_reaction randomly
            drugs = [
                "cisplatin",
                "bevacizumab",
                "abeciclib",
                "olaparib",
                "paclitaxel",
            ]
            # with 20% probability, sample one of the drugs
            if random.random() < 0.2:
                _props["severe_adverse_reaction_to"] = random.choice(drugs)

            yield (
                patient_id,
                "patient",
                _props,
            )
