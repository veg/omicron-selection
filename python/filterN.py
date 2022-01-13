import sys, csv, json
from collections import Counter

from Bio import SeqIO

max_N = 1000

if len (sys.argv) > 2:
    max_N = int (sys.argv[2])

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len ([k for k in record.seq if k == 'N']) >= max_N : continue
        print (">%s\n%s\n" % (record.id, record.seq))
                
