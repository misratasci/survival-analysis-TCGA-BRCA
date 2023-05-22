import time
import os
import pandas as pd
import xml.etree.cElementTree as et

start = time.time()

df = pd.DataFrame(columns=["Vital Status", "Days to last followup", "Days to death"])

for folder in os.listdir("clinical_downloads")[:-1]:
    xmlfile = os.listdir(f"clinical_downloads/{folder}")[0]
    submitter_id = xmlfile.split(".")[2]
    tree = et.parse(f"clinical_downloads/{folder}/{xmlfile}")
    patient = tree.getroot()[1]
    follow_ups = patient.find("{http://tcga.nci/bcr/xml/clinical/brca/2.7}follow_ups")
    if len(follow_ups) != 0:
        vital_status = follow_ups[-1].find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}vital_status").text
        max_days = 0
        for f in follow_ups:
            try: #bazı followuplarda days to last followup verisi yok, nonetype çıkıyor diye buraya try dedim. 
                if max_days < int(f.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_last_followup").text):
                    max_days = int(f.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_last_followup").text)
            except:
                pass
        if max_days != 0:
            days_to_last_followup = max_days
        else:
            days_to_last_followup = f.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_last_followup").text
        
        days_to_death = follow_ups[-1].find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_death").text
    else:
        vital_status = patient.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}vital_status").text
        days_to_last_followup = patient.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_last_followup").text
        days_to_death = patient.find("{http://tcga.nci/bcr/xml/clinical/shared/2.7}days_to_death").text
    df.loc[submitter_id] = [vital_status, days_to_last_followup, days_to_death]
df = df.fillna("None")
df.to_csv("clinical_table.txt")
            



end = time.time()
print(end - start, "secs", sep=" ")