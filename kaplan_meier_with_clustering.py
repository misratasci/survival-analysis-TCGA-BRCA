import pandas as pd 
import time
import numpy as np
from sklearn.cluster import SpectralCoclustering
from matplotlib import pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from get_gene_set import GeneSet, get_rand_pathway

start = time.time()
#pathway = get_rand_pathway()
pathway = 'GO_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY'
"""
these pathways have p value < 0.001 :
['GO_ACETYL_COA_BIOSYNTHETIC_PROCESS'
 'GO_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY'
 'GO_RESPONSE_TO_PHORBOL_13_ACETATE_12_MYRISTATE'
 'GO_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE'
 'GO_POSITIVE_REGULATION_OF_DNA_LIGATION'
 'GO_PEPTIDYL_PROLINE_HYDROXYLATION_TO_4_HYDROXY_L_PROLINE'
 'GO_NEGATIVE_REGULATION_OF_KERATINOCYTE_DIFFERENTIATION'
 'GO_MEIOSIS_I_CELL_CYCLE_PROCESS' 'GO_ACETYL_COA_METABOLIC_PROCESS'
 'GO_REGULATION_OF_DNA_LIGATION']
"""
gs = GeneSet(pathway)

mRNA_table = pd.read_csv("z-score_normalized_mRNA_table.txt", index_col=[0], header=0, usecols=["Unnamed: 0"] +  gs.gene_set)


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


kmf = KaplanMeierFitter()
ax = plt.subplot()

ax.set_ylabel('survival ratio')

low_label = f"Group 1: ({group_1_df.shape[0]} patients)"
kmf.fit(group_1_df["Duration"], event_observed=group_1_df["Event_Observed"], label = low_label)
ax = kmf.survival_function_.plot(ax=ax)
high_label = f"Group 2: ({group_2_df.shape[0]} patients)"
kmf.fit(group_2_df["Duration"], event_observed=group_2_df["Event_Observed"], label = high_label)
ax = kmf.survival_function_.plot(ax=ax)
ax.set_xlabel('days')
ax.set_ylim(ymin=0)
logrank = logrank_test(group_1_df["Duration"], group_2_df["Duration"], event_observed_A=group_1_df["Event_Observed"], event_observed_B=group_2_df["Event_Observed"])
pathway_title = pathway.replace("_", " ").title()
p_value_text = '{:.5f}'.format(logrank.p_value)
if float(p_value_text) < 10**(-5):
    p_value_text = "< 10^-5"
plt.title(f"{pathway_title}\np-value {p_value_text}")
plt.show()

end = time.time()
print(end - start, "secs", sep=" ")