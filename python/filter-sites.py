import sys, csv, json
from collections import Counter

from Bio import SeqIO

ranges = [int(k) for k in sys.argv[2].split (',')]

with open(sys.argv[1]) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        s = str (record.seq)
        if len (s) > 28000:
            print ('>%s\n%s\n' % (record.description.replace (' ', '_').replace ('__','_'), s[ranges[0]:ranges[1]]))
        
                
