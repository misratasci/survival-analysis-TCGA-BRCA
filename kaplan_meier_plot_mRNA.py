from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

id_to_gene_name = pd.read_csv("protein_coding_mRNA_ids_names.txt", index_col="gene_name")

#gene_name = input("Gene name: ")
gene_name = "RB1" #TP53, RB1, BRCA1/2..
mRNA_id = id_to_gene_name.loc[gene_name].values[0]

#patient_group_count = int(input("Hastaları kaça bölelim? "))
patient_group_count = 2


mRNAs = pd.read_csv("mRNA_table.txt", usecols=["Unnamed: 0", mRNA_id])
clinicals = pd.read_csv("clinical_table.txt")
df = mRNAs.merge(clinicals, how='inner', on="Unnamed: 0")
df["Duration"] = np.where(df["Vital Status"] == 'Dead', df["Days to death"], df["Days to last followup"])
df = df[df["Duration"] != "None"]
df["Duration"] = df["Duration"].astype(float)
df["Event_Observed"] = np.where(df["Vital Status"] == 'Dead', True, False)

df = df.sort_values(by=mRNA_id)
a = int(len(df)/patient_group_count)
df_list = [df[a*i:a*(i+1)] for i in range(patient_group_count)]
df_list[-1] = pd.concat([df_list[-1], df[a*patient_group_count:]], ignore_index=True)

kmf = KaplanMeierFitter()
ax = plt.subplot()
print(df_list)
ax.set_ylabel('survival ratio')
if len(df_list) == 2:
    low_label = f"{gene_name} low ({df_list[0].shape[0]} patients)"
    kmf.fit(df_list[0]["Duration"], event_observed=df_list[0]["Event_Observed"], label = low_label)
    ax = kmf.survival_function_.plot(ax=ax)
    high_label = f"{gene_name} high ({df_list[1].shape[0]} patients)"
    kmf.fit(df_list[1]["Duration"], event_observed=df_list[1]["Event_Observed"], label = high_label)
    ax = kmf.survival_function_.plot(ax=ax)
    ax.set_xlabel('days')
    ax.set_ylim(ymin=0)
    logrank = logrank_test(df_list[0]["Duration"], df_list[1]["Duration"], event_observed_A=df_list[0]["Event_Observed"], event_observed_B=df_list[1]["Event_Observed"])
    p_value_text = '{:.5f}'.format(logrank.p_value)
    if float(p_value_text) < 10**(-5):
        p_value_text = "< 10^-5"
    plt.title(f"p-value {p_value_text}")

elif len(df_list) == 3:
    low_label = f"{gene_name} low ({df_list[0].shape[0]} patients)"
    kmf.fit(df_list[0]["Duration"], event_observed=df_list[0]["Event_Observed"], label = low_label)
    ax = kmf.survival_function_.plot(ax=ax)
    mid_label = f"{gene_name} mid ({df_list[1].shape[0]} patients)"
    kmf.fit(df_list[1]["Duration"], event_observed=df_list[1]["Event_Observed"], label = mid_label)
    ax = kmf.survival_function_.plot(ax=ax)
    low_label = f"{gene_name} high ({df_list[2].shape[0]} patients)"
    kmf.fit(df_list[2]["Duration"], event_observed=df_list[2]["Event_Observed"], label = low_label)
    ax = kmf.survival_function_.plot(ax=ax)
    ax.set_xlabel('days')
    ax.set_ylim(ymin=0)

else:
    for d in df_list:
        label = f"{gene_name}, mRNA exp: {round(float(d[mRNA_id].values[0]), 2)} - {round(float(d[mRNA_id].values[-1]), 2)}, ({d.shape[0]} patients) "
        kmf.fit(d["Duration"], event_observed=d["Event_Observed"], label = label)
        ax = kmf.survival_function_.plot(ax=ax)
    ax.set_xlabel('days')
    ax.set_ylim(ymin=0)

plt.show()
