import pandas as pd 
import numpy as np
from sklearn.cluster import SpectralCoclustering
import seaborn as sb
from matplotlib import pyplot as plt
from get_gene_set import GeneSet, get_rand_pathway

pathway = "KEGG_PATHWAYS_IN_CANCER"
gs = GeneSet(pathway)
mRNA_table = pd.read_csv("z-score_normalized_mRNA_table.txt", usecols=["Unnamed: 0"] +  gs.gene_set, index_col=[0], header=0)
mRNA_table = mRNA_table.transpose()

model = SpectralCoclustering(n_clusters=2)
model.fit(mRNA_table)

fit_data = mRNA_table.iloc[np.argsort(model.row_labels_)]
fit_data = fit_data.iloc[:, np.argsort(model.column_labels_)]

#plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
#plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

#fig, ax = plt.subplots()
heatmap = sb.heatmap(fit_data, vmin=-3.5, vmax=3.5, cmap="RdBu_r", xticklabels=False)
plt.xlabel('')
plt.ylabel('')
plt.title(gs.pathway.title().replace("_", " "), loc='left')

plt.show()
