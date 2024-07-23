import random
import pandas as pd

bp_genes = pd.read_csv(
    "data/oncodash files/GeneToBiologicalProcess-part000.csv",
    sep=";",
    names=["Gene.MANE", "BiologicalProcess", "Label"],
)
gene_list = bp_genes["Gene.MANE"].unique().tolist()
gene_list = [gene.replace(":gene_hugo", "") for gene in gene_list]

# add oncokb genes
oncokb_genes = pd.read_csv(
    "data/oncokb_biomarker_drug_associations.tsv", sep="\t"
)
oncokb_genes = oncokb_genes["Gene"].drop_duplicates().tolist()
gene_list.extend(oncokb_genes)
gene_list = list(set(gene_list))

tissues = ["Per", "Ova", "Ome", "Asc", "Lum"]
data_variants = pd.read_csv("data/filtered_variants.csv", sep="\t")

# create list of patients (strings) from patient1 to patient20
patients = ["patient" + str(i) for i in range(1, 21)]

# create synthetic sequence variant data
synthetic_data_variants = []
for patient in patients:
    for gene in gene_list:
        gene_data = data_variants[data_variants["Gene.MANE"] == gene]
        if gene_data.empty:
            continue
        # randomly continue with probability prob
        prob = 0.1
        if prob < random.random():
            continue
        random_patient = gene_data.sample(n=1, random_state=42)[
            "patient"
        ].values[0]
        selection = gene_data[gene_data["patient"] == random_patient]
        for index, row in selection.iterrows():
            random_tissues = random.sample(tissues, random.randint(1, 5))
            patient_with_tissues = [
                patient + "_" + tissue for tissue in random_tissues
            ]
            selection.at[index, "samples"] = ";".join(patient_with_tissues)
        selection.loc[:, "patient"] = patient
        # randomly select 10% of selection, but at least one row
        selection = selection.sample(
            n=max(1, int(len(selection) * 0.1)), random_state=42
        )
        synthetic_data_variants.append(selection)
synthetic_data_variants = pd.concat(synthetic_data_variants)
synthetic_data_variants.to_csv(
    "data/synthetic_variants.csv", index=False, sep="\t"
)

# create synthetic copy number variant data
data_cns = pd.read_csv("data/filtered_cns.csv", sep="\t")
synthetic_data_cns = []
for patient in patients:
    for gene in gene_list:
        gene_data_cns = data_cns[data_cns["Gene"] == gene]
        if gene_data_cns.empty:
            continue
        # randomly continue with probability prob
        prob = 0.1
        if prob < random.random():
            continue
        alleles = gene_data_cns["nMajor"] + gene_data_cns["nMinor"]
        gene_data_cns = gene_data_cns[alleles > 4]
        random_patient_cns = gene_data_cns.sample(n=1, random_state=42)[
            "sample"
        ].values[0]
        selection_cns = gene_data_cns[
            gene_data_cns["sample"] == random_patient_cns
        ]
        selection_cns.loc[:, "sample"] = patient
        synthetic_data_cns.append(selection_cns)
synthetic_data_cns = pd.concat(synthetic_data_cns)
synthetic_data_cns.to_csv("data/synthetic_cns.csv", index=False, sep="\t")

# create synthetic clinical data
data_clinical = pd.read_csv(
    "data/clinical_raw.tsv", sep="\t", encoding="utf-16"
)
synthetic_data_clinical = []
for patient in patients:
    # construct a synthetic row by selecting random rows from the clinical data
    # for each individual column of the original data
    synthetic_row = {}
    synthetic_row["patient"] = patient
    for column in data_clinical.columns:
        # select a random row from the original data
        if column in [
            "Patient ID",
            "Patient card::Patient cohort code_Patient Card",
            "Patient card::Publication code",
            "Patient card::Patient with sequenced samples_sent_ proceed",
        ]:
            continue
        random_row = data_clinical.sample(n=1)
        # set the value of the synthetic row to the value of the random row
        data = random_row[column].values[0]
        if isinstance(data, str):
            data = data.replace("\n", ", ")
        synthetic_row[column] = data

    # set the patient ID to the current patient
    synthetic_data_clinical.append(synthetic_row)
synthetic_data_clinical = pd.DataFrame(synthetic_data_clinical)
synthetic_data_clinical.to_csv(
    "data/synthetic_clinical.csv", index=False, sep=";"
)
