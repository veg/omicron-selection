import sys, csv
from collections import Counter

from Bio import SeqIO


ids = set ()
seqs = set ()

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if not record.id in ids:
            ids.add (record.id)
            s = str (record.seq)
            if not s in seqs:
                seqs.add (s)
                print (">%s\n%s\n" % (record.id, s))
           


        
                
