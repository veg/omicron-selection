import sys, csv, json
from collections import Counter

from Bio import SeqIO

annotation = {}
with open (sys.argv[2]) as fh:
    pangolin = csv.reader (fh)
    next (pangolin)
    for line in pangolin:
        ID = line[0].split ('.')[0]
        annotation[ID] = [line [4], line[1]]
        

max_N = 1000

if len (sys.argv) > 4:
    max_N = int (sys.argv[4])

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        tag = annotation[record.id]
        if (tag [1] == sys.argv[3]):
            if len ([k for k in record.seq if k == 'N']) >= max_N : continue
            print (">%s\n%s\n" % (record.id, record.seq))
                
