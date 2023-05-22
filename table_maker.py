import numpy as np
import pandas as pd 
import os
import time
#588 saniyede bütün dosyaları yazdı. 
start = time.time()

#mRNA_ids = pd.read_csv('protein_coding_mRNA_ids.txt', header=None).sort_values(0, ignore_index=True)
mRNA_ids = np.loadtxt("protein_coding_mRNA_ids.txt", dtype = "U")
filename_to_id = pd.read_csv("filename_submitter_id.txt", header=None, sep='\t').set_index(0)

f = os.listdir("fpkm_downloads")[0]
fol_dir = f"fpkm_downloads/{f}"
first_txt_dir = f"{fol_dir}/{os.listdir(fol_dir)[0]}"
all_mRNAs = np.loadtxt(first_txt_dir, dtype="U", usecols=0)
mRNA_mask = np.isin(all_mRNAs, mRNA_ids)
df_cols = all_mRNAs[mRNA_mask]
df = pd.DataFrame(columns=df_cols)
#bütün dosyaların ilk sütunundaki mrna id lerinin sırası aynı diye varsaydım. (umarım öyledir)
for folder in os.listdir("fpkm_downloads")[:-1]:
    txtfile = os.listdir(f"fpkm_downloads/{folder}")[0]
    txtdir = f"fpkm_downloads/{folder}/{txtfile}"
    #arr = pd.read_csv(txtdir, sep="\t", header=None).sort_values(0, ignore_index=True)
    arr = np.loadtxt(txtdir, usecols=1)
    arr = arr[mRNA_mask]
    sub_id = filename_to_id.loc[txtfile].values[0]
    df.loc[sub_id] = arr

df.to_csv("mRNA_table.txt")


end = time.time()
print(end - start, "secs", sep=" ")