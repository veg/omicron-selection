import sys, csv, json
from   operator import itemgetter, attrgetter

from Bio import SeqIO

sequences = {}

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        S = str (record.seq)
        if not S in sequences:
            sequences[S] = set ()
            
        sequences[S].add (record.id)
        
output = []

for seq, ids in sequences.items():
    m = [k for k in ids]
    output.append ({
        'size' : len (ids),
        'members'  : m,
        'centroid' : ">%s||%d\n%s" % (m[0], len (ids), seq)
    })       

json.dump (output, sys.stdout)