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

synthetic_data_variants = []
for patient in [
    "patient1",
    "patient2",
    "patient3",
    "patient4",
    "patient5",
]:
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
        synthetic_data_variants.append(selection)
synthetic_data_variants = pd.concat(synthetic_data_variants)
synthetic_data_variants.to_csv(
    "data/synthetic_variants.csv", index=False, sep="\t"
)

data_cns = pd.read_csv("data/filtered_cns.csv", sep="\t")
synthetic_data_cns = []
for patient in [
    "patient1",
    "patient2",
    "patient3",
    "patient4",
    "patient5",
]:
    for gene in gene_list:
        gene_data_cns = data_cns[data_cns["Gene"] == gene]
        if gene_data_cns.empty:
            continue
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
