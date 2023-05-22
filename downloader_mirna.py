import json
import requests
import re
import time
start = time.time()
fields = [   
    "file_name",
    "cases.submitter_id",
    #"cases.samples.sample_type",
    #"cases.disease_type",
    #"cases.project.project_id",
    "data_type",
    #"experimental_strategy"
    #"cases.case_id"

    ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/files"

#TCGA-BRCA projesi, FPKM ve primer tümör olsun:
filters = {
    "op": "and", 
    "content":[
        {
            "op": "=",
            "content":{
            "field": "cases.project.project_id",
            "value": ["TCGA-BRCA"]
            }
        },
        {
            "op": "=",
            "content":{
            "field": "data_type",
            "value": ["miRNA Expression Quantification"]
            }
        },
        {
            "op": "=",
            "content":{
            "field": "analysis.workflow_type",
            "value": ["BCGSC miRNA Profiling"]
            }
        },
        {
            "op": "=",
            "content":{
            "field": "cases.samples.sample_type",
            "value": ["Primary Tumor"]
            }
        }
    ]
    }

file_no = "2000" #hepsini indirmek için sayıyı 1102den büyük yap


params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "JSON",
    "size": file_no
    }

response = requests.get(files_endpt, params = params)
filename_to_submitter_id_dict = {}
file_ids = []
entry_count = 0
for entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    print(entry)
    file_ids.append(entry["id"])
    entry_count += 1
    filename_to_submitter_id_dict[entry["file_name"]] = entry["cases"][0]["submitter_id"]
print(entry_count)

#burası indirdiğin dosya isimlerini case submitter id ile eşleştiren bi txt dosyası yazıyor.
with open("filename_submitter_id_mirna.txt", "w") as text_file:
    for key in filename_to_submitter_id_dict:
        text_file.write(f"{key}\t{filename_to_submitter_id_dict[key]}\n")


#burası indiriyor
data_endpt = "https://api.gdc.cancer.gov/data"

params = {"ids": file_ids}
response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})
response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)


end = time.time()
print(end - start, "secs", sep=" ")
