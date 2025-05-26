import requests
import json
import sys
import pandas as pd
from io import StringIO

sys.stderr = open(snakemake.log[0], "w") 

primary_site = snakemake.params[0]
sample_type = snakemake.params[1]

if (sample_type == "cancer"):
    sample_type = "primary tumor"
elif (sample_type == "normal"):
    sample_type = "solid tissue normal"
else:
    sys.exit("Incorrect sample type")

print("Getting data for: " + primary_site + ", " + sample_type)

# Fields for the query
fields = [
    "cases.case_id",
    "file_name"
    ]

fields = ",".join(fields)

# Endpoints used
files_endpt = "https://api.gdc.cancer.gov/files"
manifest_endpt = "https://api.gdc.cancer.gov/manifest"

# Filters to get ASCAT cnvs data
filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.primary_site",
            "value": [primary_site]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.samples.sample_type",
            "value": [sample_type]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.analysis.workflow_type",
            "value": ["ASCAT2"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_type",
            "value": ["Gene Level Copy Number"]
            }
        }
    ]
}

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "tsv",
    "size": "2000"
    }

# Getting ASCAT2 data
response = requests.get(files_endpt, params = params)
data = response.content.decode("utf-8")

print("ASCAT2 query done")

if not data.strip():
    open(snakemake.output[0], 'a').close()
else:
  df_ascat = pd.read_csv(StringIO(data), sep ="\t")
  df_ascat.drop_duplicates(subset=["cases.0.case_id"], inplace=True)
  df_ascat.to_csv(snakemake.output[0], sep="\t", index=False)  

  # Getting ASCAT2 files manifest for future download
  params = { "ids" : df_ascat["id"].tolist() }
  response = requests.post(manifest_endpt, data= json.dumps(params), headers = {"Content-Type":
      "application/json"})
  
print("Got ASCAT2 manifest")

f = open(snakemake.output[1], "w")
f.write(response.content.decode("utf-8"))
f.close()

print("ASCAT2 manifest written")
