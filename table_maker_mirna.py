import numpy as np
import pandas as pd 
import os
import time
#30 saniyede çalıştı.
start = time.time()

filename_to_id = pd.read_csv("filename_submitter_id_mirna.txt", header=None, sep='\t').set_index(0)

f = os.listdir("mirna_downloads")[0]
fol_dir = f"mirna_downloads/{f}"
first_txt_dir = f"{fol_dir}/{os.listdir(fol_dir)[0]}"
mirna_ids = np.loadtxt(first_txt_dir, dtype="U", usecols=0, skiprows=1)
df = pd.DataFrame(columns=mirna_ids)

#bütün dosyaların ilk sütunundaki mirna id lerinin sırası aynı diye varsaydım. (umarım öyledir)

for folder in os.listdir("mirna_downloads")[:-1]:
    txtfile = os.listdir(f"mirna_downloads/{folder}")[0]
    txtdir = f"mirna_downloads/{folder}/{txtfile}"
    arr = np.loadtxt(txtdir, usecols=2, skiprows=1)
    sub_id = filename_to_id.loc[txtfile].values[0]
    df.loc[sub_id] = arr

df.to_csv("miRNA_table.txt")
print(df)

end = time.time()
print(end - start, "secs", sep=" ")