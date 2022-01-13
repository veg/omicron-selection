import json, sys
from collections import Counter

sites = [int (k) for k in sys.argv[2].split (",")]

with open (sys.argv[1], "r") as fh:
    slac = json.load (fh)
    records = []
    
    counter = Counter()
    counter_codon = Counter()
    
    for branch, codons in slac["branch attributes"]["0"].items():
        if branch [0:4] != "Node":
            aa = "".join ([k for k in [codons["amino-acid"][0][s-1] for s in sites] if k!='-' and len(k) == 1])
            if len (aa) == len (sites):
                counter[aa] += 1
            codon = "".join ([k for k in [codons["codon"][0][s-1] for s in sites] if k!='---' and len(k) == 3])
            if len (codon) == 3*len (sites):
                counter_codon[codon] += 1
            
for k, c in counter.items():
    print (k, c)
        
        
print ("\nCODONS\n")        
for k, c in counter_codon.items():
    print (k, c)
