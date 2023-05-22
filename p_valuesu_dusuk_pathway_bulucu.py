import pandas as pd 
import time
import numpy as np
from sklearn.cluster import SpectralCoclustering
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from get_gene_set import GeneSet, get_next_pathway

start = time.time()

"""
pathway_count = 0
for line in open("pathways_with_p-values.txt").readlines(): pathway_count += 1
while True:
    pathway = get_next_pathway(pathway_count)
    #pathway = "GO_BRCA1_A_COMPLEX"
    gs = GeneSet(pathway)

    mRNA_table = pd.read_csv("z-score_normalized_mRNA_table.txt", index_col=[0], header=0, usecols=["Unnamed: 0"] +  gs.gene_set)
    cols_with_all_zeros = mRNA_table.loc[:, (mRNA_table == 0).all()].columns.tolist()
    print(mRNA_table.shape)
    mRNA_table = mRNA_table.drop(labels=cols_with_all_zeros, axis=1).dropna(axis=1)
    print(mRNA_table.shape)
    model = SpectralCoclustering(n_clusters=2)
    if not mRNA_table.empty:
        model.fit(mRNA_table)
        mRNA_table["patient_group"] = model.row_labels_

        group_1 = mRNA_table[mRNA_table["patient_group"] == 0]
        group_2 = mRNA_table[mRNA_table["patient_group"] == 1]


        df = pd.read_csv("clinical_table.txt", index_col=[0])
        df["Duration"] = np.where(df["Vital Status"] == 'Dead', df["Days to death"], df["Days to last followup"])
        df = df[df["Duration"] != "None"]
        df["Duration"] = df["Duration"].astype(float)
        df["Event_Observed"] = np.where(df["Vital Status"] == 'Dead', True, False)
        group_1_df = df.merge(group_1, how='inner', left_index=True, right_index=True)
        group_2_df = df.merge(group_2, how='inner', left_index=True, right_index=True)

        logrank = logrank_test(group_1_df["Duration"], group_2_df["Duration"], event_observed_A=group_1_df["Event_Observed"], event_observed_B=group_2_df["Event_Observed"])

        print(pathway, logrank.p_value)
        with open("pathways_with_p-values.txt", "a") as output:
            output.write(f"{pathway} - {logrank.p_value}\n")
            pathway_count += 1
        cur_time = time.time()
        #if cur_time - start > 300:
        #    break
    else:
        pathway_count += 1
        continue
"""
df = pd.read_csv("pathways_with_p-values.txt", sep=" - ", names=["pathway", "p-value"])
df = df.sort_values("p-value")
pathways = df[df["p-value"] < 0.001].iloc[:,0].values
print(pathways)
end = time.time()
print(end - start, "secs", sep=" ")