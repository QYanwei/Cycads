import json
import pprint
data = [{"studentid": 1, "name": "ABC",
         "subjects": ["Python", "Data Structures"]},
        {"studentid": 2, "name": "PQR",
         "subjects": ["Java", "Operating System"]}]

with open("../test/filename.json", "w") as write_file:
    pprint.pprint(data, write_file, indent=2, width=50, compact=True)
    
fastq_json_file_path = '../test/ecoli.seq.json'

with open(fastq_json_file_path, 'r', encoding='utf-8') as jsonfile:
    s = jsonfile.read()
    js = json.loads(json.dumps(eval(s)))