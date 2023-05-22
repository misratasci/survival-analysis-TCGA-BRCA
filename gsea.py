import pandas as pd 
import time
import numpy as np
from scipy.stats import ttest_ind, hypergeom
import os

start = time.time()

mRNA_table = pd.read_csv("mRNA_table_with_gene_names.txt", header=0)
demographics = pd.read_csv("age_gender_race_table.txt")

df = demographics.merge(mRNA_table, on=["Unnamed: 0"]).set_index(["Unnamed: 0"])
white_patients = df[df["Race"] == "WHITE"].iloc[:, 3:]
#print(white_patients)
other_patients = df[df["Race"] != "WHITE"].iloc[:, 3:]

pvalues = pd.Series(ttest_ind(white_patients, other_patients, axis=0).pvalue)
genes = pd.Series(mRNA_table.columns[1:])
gene_no = int(len(genes)*.05)
genes_with_low_pvalues = pd.concat([genes, pvalues], axis=1).sort_values(by=1, ignore_index=True).iloc[:gene_no, 0]

def read_gene_set_file(filename):
    df = pd.read_csv(filename, sep="\n", header=None)
    df = df[0].str.split("\t", expand=False)
    ind = df.apply(lambda x: x[0])
    genes = df.apply(lambda y: y[2:])
    df = pd.concat([genes], axis=1).set_index(ind)
    return df

df_list = []
for file in os.listdir("gene_sets"):
    df_list.append(read_gene_set_file(f"gene_sets/{file}"))
pathways = pd.concat(df_list, axis=0)

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

pathways["Interesting Gene Count"] = pathways.apply(lambda x: len(intersection(x[0], genes_with_low_pvalues.tolist())), axis=1)
pathways["Hypergeom pvalue"] = pathways.apply(lambda x: hypergeom.cdf(x["Interesting Gene Count"], len(genes), gene_no, len(x[0])), axis=1)
pathways["Gene Count"] = pathways.apply(lambda y: len(y[0]), axis=1)
pathways = pathways.sort_values(by=["Hypergeom pvalue"])
output = pathways[["Gene Count", "Interesting Gene Count", "Hypergeom pvalue"]]
output.to_csv("white_vs_others_gsea.txt")
end = time.time()
print(end - start, "secs", sep=" ")