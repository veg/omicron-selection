import sys, csv, json, datetime
from collections import Counter

from Bio import SeqIO


variants = {}

byPosition = []
N          = 0
dates      = Counter()

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        N += 1
        try:
            D = datetime.datetime.strptime (record.id.split ("|")[-1], "%Y-%m-%d").strftime ("%Y%m%d")
        except Exception as e:
            D = None
            
        dates [D] += 1

        S = str (record.seq)
        for i in range (0,len(record.seq),3):
            if N == 1:
                byPosition.append (Counter())
            codon = S[i:i+3]
            byPosition[i//3][codon] += 1
           
variants ['N']      = N            
variants ['counts'] = byPosition
variants ['dates'] = dates

hbp        = []
H          = 0

mapping = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3,
    'N' : 4,
    '-' : 5
}

dotplot = []

with open(sys.argv[2]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        H += 1
        S = str (record.seq)
        
        dotplot.append ([mapping[k] if k in mapping else 6 for k in S])
        
        for i in range (0,len(record.seq),3):
            if H == 1:
                hbp.append (Counter())
            codon = S[i:i+3]
            hbp[i//3][codon] += 1
           
variants ['H']      = H            
variants ['haplotypes'] = hbp

hbpa        = []
HA          = 0

with open(sys.argv[3]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        HA += 1
        S = str (record.seq)
               
        for i in range (0,len(record.seq),3):
            if HA == 1:
                hbpa.append (Counter())
            codon = S[i:i+3]
            hbpa[i//3][codon] += 1
           
variants ['HA']             = HA            
variants ['all-haplotypes'] = hbpa
    
json.dump (variants, sys.stdout)
json.dump (dotplot, sys.stderr)

        
                
