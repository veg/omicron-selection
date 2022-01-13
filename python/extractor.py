import json, sys
from collections import Counter

with open (sys.argv[1], "r") as fh:
    slac = json.load (fh)
    records = []
    for i in range (slac["input"]["number of sites"]):
        records.append ({"I" : Counter(), "L" : Counter()})
    for branch, codons in slac["branch attributes"]["0"].items():
        if branch [0:4] == "Node":
            key = "I"
        else:
            key = "L"
        for s,aa in enumerate (codons["amino-acid"][0]):
            records[s][key][aa] += 1
        
        
json.dump (records, sys.stdout)

for i,muts in enumerate(records):
    z = [[v[0],v[1]] for v in muts['L'].items() if len (v[0]) == 1 and v[0] != '-' and v[1] >= 2]
    if len (z) > 1:
        print (i+1, z, file = sys.stderr)