import pandas as pd
import time

start = time.time()
"""
genome = pd.read_csv("gencode.v22.annotation.gtf", comment='#', sep='\t', header=None, usecols=[8], names=["col"])
pr = " gene_type \"protein_coding\";"
genome = genome[genome["col"].str.contains(pr)]
genome["gene_id"] = genome["col"].str[9:26]
genome = genome.drop_duplicates(subset="gene_id")
genome["col"].to_csv("protein_coding_mRNAs.txt", index=False, header=False)
print(genome)
"""
genome = pd.read_csv("protein_coding_mRNAs.txt",  header=None, sep="gene", usecols=[1, 4], names=["gene_id", "gene_name"])
for col in genome.columns:
    genome[col] = genome[col].str.split('"')
    genome[col] = genome[col].apply(lambda x: x[2])
genome = genome.drop_duplicates(subset="gene_name")
print(genome)
genome.to_csv("protein_coding_mRNA_ids_names.txt", index=False)
end = time.time()
print(end - start, "secs")
