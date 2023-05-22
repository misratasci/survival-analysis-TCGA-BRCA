from lifelines import KaplanMeierFitter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#miRNA_id = input("miRNA ID: ")
miRNA_id = "hsa-let-7a-1"
#patient_group_count = int(input("Hastaları kaça bölelim? "))
patient_group_count = 2

miRNAs = pd.read_csv("miRNA_table.txt", usecols=["Unnamed: 0", miRNA_id])
clinicals = pd.read_csv("clinical_table.txt")
df = miRNAs.merge(clinicals, how='inner', on="Unnamed: 0")
df["Duration"] = np.where(df["Vital Status"] == 'Dead', df["Days to death"], df["Days to last followup"])
df = df[df["Duration"] != "None"]
df["Duration"] = df["Duration"].astype(float)
df["Event_Observed"] = np.where(df["Vital Status"] == 'Dead', True, False)
df = df.sort_values(by=miRNA_id)

a = int(len(df)/patient_group_count) + 1
df_list = [df[a*i:a*(i+1)] for i in range(patient_group_count)]

kmf = KaplanMeierFitter()
ax = plt.subplot()
for d in df_list:

    label = f"miRNA exp: {round(float(d[miRNA_id].values[0]), 2)} - {round(float(d[miRNA_id].values[-1]), 2)} "
    kmf.fit(d["Duration"], event_observed=d["Event_Observed"], label = label)
    ax = kmf.survival_function_.plot(ax=ax)

plt.show()