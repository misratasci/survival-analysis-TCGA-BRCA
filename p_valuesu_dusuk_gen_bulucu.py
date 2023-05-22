import pandas as pd 
import time
import numpy as np
from sklearn.cluster import SpectralCoclustering
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from get_gene_set import GeneSet, get_rand_pathway

start = time.time()

patient_group_count = 2
mRNAs = pd.read_csv("mRNA_table_with_gene_names.txt")
clinicals = pd.read_csv("clinical_table.txt")
df = mRNAs.merge(clinicals, how='inner', on="Unnamed: 0")
df["Duration"] = np.where(df["Vital Status"] == 'Dead', df["Days to death"], df["Days to last followup"])
df = df[df["Duration"] != "None"]
df["Duration"] = df["Duration"].astype(float)
df["Event_Observed"] = np.where(df["Vital Status"] == 'Dead', True, False)
#df = df.drop(["Days to last followup", "Days to death", "Vital Status"], axis=1)
print(df)
"""
df = df.sort_values(by=mRNA_id)
a = int(len(df)/patient_group_count)
df_list = [df[a*i:a*(i+1)] for i in range(patient_group_count)]
df_list[-1] = pd.concat([df_list[-1], df[a*patient_group_count:]], ignore_index=True)

logrank = logrank_test(df_list[0]["Duration"], df_list[1]["Duration"], event_observed_A=df_list[0]["Event_Observed"], event_observed_B=df_list[1]["Event_Observed"])





while True:
    try:

        model = SpectralCoclustering(n_clusters=2)
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
    except:
        continue
    if logrank.p_value <= 0.01:
        print(pathway, logrank.p_value)
        break
"""
end = time.time()
print(end - start, "secs", sep=" ")