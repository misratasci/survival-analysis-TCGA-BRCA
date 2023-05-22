import pandas as pd 
import time
import numpy as np
from sklearn.cluster import SpectralCoclustering
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl
from get_gene_set import GeneSet, get_rand_pathway

import random
import matplotlib as mpl

start = time.time()
pathway = 'GO_POSITIVE_REGULATION_OF_DNA_LIGATION'
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
#pathway = get_rand_pathway()
gs = GeneSet(pathway)
mRNA_table = pd.read_csv("z-score_normalized_mRNA_table.txt", usecols=["Unnamed: 0"] +  gs.gene_set, index_col=[0], header=0)

demographics = pd.read_csv("age_gender_race_table.txt", index_col=[0])
merged = mRNA_table.merge(demographics, how='inner', left_index=True, right_index=True)

table_for_heatmap = merged.iloc[:,:-3].transpose()

model = SpectralCoclustering(n_clusters=2)
model.fit(table_for_heatmap)
merged["Group"] = model.column_labels_

fit_data = table_for_heatmap.iloc[np.argsort(model.row_labels_)]
fit_data = fit_data.iloc[:, np.argsort(model.column_labels_)]
merged = merged.iloc[np.argsort(model.column_labels_)]
annotation_df = merged.iloc[:,-4:]

col_colors = pd.DataFrame()
color_palettes = ["coolwarm_r", "Greens", "Spectral", "coolwarm"]
demographic_values = []
for i in range(len(annotation_df.columns)):
    col = annotation_df.iloc[:, i]
    cmap = sns.color_palette(color_palettes[i], as_cmap=True)
    color_list = [cmap(a/(len(col.unique())-1)) for a in range(len(col.unique()))]
    lut = dict(zip(np.sort(col.unique()), color_list))
    demographic_values.append(np.sort(col.unique()))
    colors = col.map(lut)
    col_colors = pd.concat([col_colors,colors],axis=1)

heatmap = sns.clustermap(fit_data, vmin=-3.5, vmax=3.5, cmap="RdBu_r", xticklabels=False, 
row_cluster=False, col_cluster=False, z_score=None, 
col_colors=col_colors, linewidths=0)
plt.xlabel('')
plt.ylabel('')
heatmap.fig.suptitle(gs.pathway.title().replace("_", " "))

cbar_axes = []
for i in range(3):
    cb = heatmap.fig.add_axes((-.155,.5-i*.25,.2,.33))
    cbar_axes.append(cb)

gender_cbar_im = cbar_axes[0].imshow(col_colors.apply(lambda x: x[0]).to_numpy(dtype="f"), vmin=0, vmax=len(col_colors.iloc[:,0].unique()), cmap=plt.cm.get_cmap(color_palettes[0], len(col_colors.iloc[:,0].unique())))
cbar_axes[0].set_axis_off()
cbar_axes[0].set_visible(False)

age_cbar_im = cbar_axes[1].imshow(col_colors.apply(lambda x: x[0]).to_numpy(dtype="f"), vmin=0, vmax=max(demographic_values[1]), cmap=plt.cm.get_cmap(color_palettes[1], len(col_colors.iloc[:,1].unique())))
cbar_axes[1].set_axis_off()
cbar_axes[1].set_visible(False)

race_cbar_im = cbar_axes[2].imshow(col_colors.apply(lambda x: x[0]).to_numpy(dtype="f"), vmin=0, vmax=len(col_colors.iloc[:,2].unique()), cmap=plt.cm.get_cmap(color_palettes[2], len(col_colors.iloc[:,2].unique())))
cbar_axes[2].set_axis_off()
cbar_axes[2].set_visible(False)


gender_cbar = plt.colorbar(gender_cbar_im, aspect=5, ax=cbar_axes[0], shrink=.4, label=col_colors.columns[0])
gender_cbar.set_ticks([i+0.5 for i in range(len(annotation_df.iloc[:,0].unique().tolist()))])
gender_cbar.set_ticklabels(demographic_values[0].tolist())

age_cbar = plt.colorbar(age_cbar_im, aspect=5, ax=cbar_axes[1], shrink=.4, label=col_colors.columns[1])
age_cbar.set_ticks([int(90/5)*i for i in range(6)])
age_cbar.set_ticklabels(demographic_values[1].tolist()[0:len(demographic_values[1].tolist()):int(len(demographic_values[1].tolist())/5)] + [max(demographic_values[1])])

race_cbar = plt.colorbar(race_cbar_im, aspect=5, ax=cbar_axes[2], shrink=.4, label=col_colors.columns[2])
race_cbar.set_ticks([i+0.5 for i in range(len(col_colors.iloc[:,2].unique()))])
racelist = []
for race in demographic_values[2].tolist():
    race = race.replace("OR ", "OR\n")
    race = race.replace("None", "UNKNOWN")
    racelist.append(race)
print(racelist)
race_cbar.set_ticklabels(racelist)

race_cbar.ax.tick_params(labelsize=7)
#plt.savefig("clustermap.png")
plt.show()

end = time.time()
print(end - start, "secs", sep=" ")
