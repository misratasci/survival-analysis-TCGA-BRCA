import pandas as pd 
import os
import numpy as np
from sklearn.cluster import SpectralCoclustering
import math

def read_gene_set_file(filename):
    df = pd.read_csv(filename, sep="\n", header=None)
    df = df[0].str.split("\t", expand=True).drop(columns=[1]).set_index(0)
    return df

def get_rand_pathway():
    df = pd.read_csv("all_pathways.txt", header=None).squeeze()
    return df.sample().values[0]

def get_next_pathway(count=0):
    df = pd.read_csv("all_pathways.txt", header=None).squeeze()
    return df.values[count]

df_list = []
for file in os.listdir("gene_sets"):
    df_list.append(read_gene_set_file(f"gene_sets/{file}"))
master_df = pd.concat(df_list, axis=0).dropna(axis=1, how='all')


class GeneSet():
    def __init__(self, pathway):
        self.pathway = pathway
        print(self.pathway)
        self.gene_set = master_df.loc[pathway].dropna().tolist()
        genes_in_table = pd.read_csv("mRNA_table_with_gene_names.txt", nrows=1, header=None).squeeze().dropna()
        self.gene_set = list(set(self.gene_set).intersection(genes_in_table))

        mRNA_table = pd.read_csv("mRNA_table_with_gene_names.txt", index_col=[0], header=0, usecols=["Unnamed: 0"] +  self.gene_set)
        cols_with_all_zeros = mRNA_table.loc[:, (mRNA_table == 0).all()].columns.tolist()
        mRNA_table = mRNA_table.drop(labels=cols_with_all_zeros, axis=1)
        if (not mRNA_table.empty) and len(mRNA_table.columns) > 1:
            #print(mRNA_table, mRNA_table.shape)
            mRNA_table = mRNA_table.replace(0, 0.000001)
            mRNA_table = mRNA_table.replace([np.inf, -np.inf], np.nan).dropna(how="all")
            #print(mRNA_table, mRNA_table.shape)
            model = SpectralCoclustering(n_clusters=2)
            model.fit(mRNA_table)
            mRNA_table["patient_group"] = model.row_labels_
            
            group_1 = mRNA_table[mRNA_table["patient_group"] == 0].drop(labels="patient_group", axis=1)
            group_2 = mRNA_table[mRNA_table["patient_group"] == 1].drop(labels="patient_group", axis=1)

            mean_1 = group_1.mean(axis=0)
            mean_2 = group_2.mean(axis=0)
            self.log2_fold_change = np.log2(mean_1.divide(mean_2)).abs().sort_values(ascending=False)
            
            self.max_genes = 100
            if len(self.gene_set) > self.max_genes:
                self.gene_set = list(self.log2_fold_change[:self.max_genes].index.values)
        else:
            self.gene_set = []    


#pathway = "GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY"
#print(len(master_df.loc[pathway].dropna().tolist()), master_df.loc[pathway].dropna())
#g = GeneSet(get_rand_pathway())
#print(g.pathway, g.gene_set)
