import time
import os
import pandas as pd
import xml.etree.cElementTree as et

start = time.time()

df = pd.DataFrame(columns=["Gender", "Age", "Race"])

for folder in os.listdir("clinical_downloads")[:-1]:
    xmlfile = os.listdir(f"clinical_downloads/{folder}")[0]
    submitter_id = xmlfile.split(".")[2]
    tree = et.parse(f"clinical_downloads/{folder}/{xmlfile}")
    patient = tree.getroot()[1]
    gender = patient.find("{http://tcga.nci/bcr/xml/shared/2.7}gender").text
    age = patient.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}age_at_initial_pathologic_diagnosis").text
    race_list =patient.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}race_list")
    race =race_list.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}race").text
    df.loc[submitter_id] = [gender, age, race]
df = df.fillna("None")
#print(df.isna().sum())
df.to_csv("age_gender_race_table.txt")
print(df)



end = time.time()
print(end - start, "secs", sep=" ")